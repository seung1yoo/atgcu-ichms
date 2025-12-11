import logging
from pathlib import Path
from typing import Tuple, Dict, List, Optional

import pandas as pd


class VCFParser:
    """
    VCF 파일 파싱 클래스 (apt-vcf, imputed-vcf 지원)
    """

    NO_CALL = "./."

    def __init__(self, file_path: str, logger: logging.Logger = None):
        self.file_path = Path(file_path)
        self.logger = logger or logging.getLogger(__name__)
        self.meta_lines: List[str] = []
        self.header: List[str] = []
        self.samples: List[str] = []
        self.data: pd.DataFrame = pd.DataFrame()
        self.position_index: Dict[Tuple[str, str], Dict[str, Dict[str, str]]] = {}
        self.rsid_index: Dict[str, Dict[str, Dict[str, str]]] = {}
        self.indexed_variant_count: int = 0
        self.tabix_hits: set = set()
        self.stream_hits: set = set()

        self.logger.info(f"VCFParser initialized with file: {self.file_path}")

    def parse_header_only(self) -> None:
        if self.header:
            return

        meta_lines = []
        header_line = None

        import gzip

        open_func = gzip.open if str(self.file_path).endswith(".gz") else open
        mode = "rt" if str(self.file_path).endswith(".gz") else "r"

        with open_func(self.file_path, mode, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("##"):
                    meta_lines.append(line)
                elif line.startswith("#CHROM"):
                    header_line = line
                    break

        self.meta_lines = meta_lines

        if header_line:
            self.header = header_line.lstrip("#").split("\t")
            if "FORMAT" in self.header:
                format_idx = self.header.index("FORMAT")
                self.samples = self.header[format_idx + 1 :]

        self.logger.info(f"Parsed header only: {len(self.samples)} samples detected")

    def parse(self) -> pd.DataFrame:
        self.logger.info(f"Parsing VCF file: {self.file_path}")

        if not self.file_path.exists():
            raise FileNotFoundError(f"VCF file not found: {self.file_path}")

        meta_lines = []
        header_line = None
        data_lines = []

        import gzip

        open_func = gzip.open if str(self.file_path).endswith(".gz") else open
        mode = "rt" if str(self.file_path).endswith(".gz") else "r"

        with open_func(self.file_path, mode, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("##"):
                    meta_lines.append(line)
                elif line.startswith("#CHROM"):
                    header_line = line
                elif line:
                    data_lines.append(line)

        self.meta_lines = meta_lines

        if header_line:
            self.header = header_line.lstrip("#").split("\t")
            if "FORMAT" in self.header:
                format_idx = self.header.index("FORMAT")
                self.samples = self.header[format_idx + 1 :]

        self.logger.info(f"Found {len(self.samples)} samples in VCF")
        self.logger.info(f"Parsing {len(data_lines)} variant lines...")

        if data_lines:
            data_rows = [line.split("\t") for line in data_lines]
            self.data = pd.DataFrame(data_rows, columns=self.header)
            if "CHROM" in self.data.columns:
                self.data["CHROM"] = self.data["CHROM"].str.replace("chr", "", regex=False)

        self.logger.info(f"Parsed {len(self.data)} variants")
        return self.data

    def get_samples(self) -> List[str]:
        return self.samples

    def get_variant_by_rsid(self, rsid: str) -> pd.Series:
        if rsid in self.rsid_index:
            entry = self.rsid_index[rsid]
            row = {
                "CHROM": entry["CHROM"],
                "POS": entry["POS"],
                "REF": entry["REF"],
                "ALT": entry["ALT"],
                "ID": rsid,
            }
            row.update(entry["samples"])
            return pd.Series(row)

        if "ID" in self.data.columns:
            matches = self.data[self.data["ID"] == rsid]
            if len(matches) > 0:
                return matches.iloc[0]
        return None

    def get_variant_by_position(self, chrom: str, pos: str) -> pd.Series:
        chrom = str(chrom).replace("chr", "")
        key = (chrom, str(pos))
        if key in self.position_index:
            entry = self.position_index[key]
            row = {
                "CHROM": chrom,
                "POS": str(pos),
                "REF": entry["REF"],
                "ALT": entry["ALT"],
            }
            row.update(entry["samples"])
            return pd.Series(row)

        if "CHROM" in self.data.columns and "POS" in self.data.columns:
            matches = self.data[(self.data["CHROM"] == chrom) & (self.data["POS"] == str(pos))]
            if len(matches) > 0:
                return matches.iloc[0]
        return None

    def decode_genotype(self, gt: str, ref: str, alt: str) -> Tuple[str, str]:
        if gt in ["./.", ".", "./."]:
            return (None, None)

        gt_core = gt.split(":", 1)[0] if gt else gt
        separator = "|" if "|" in gt_core else "/"
        parts = gt_core.split(separator)

        if len(parts) != 2:
            return (None, None)

        alleles = [ref, alt]
        try:
            a1_idx = int(parts[0])
            a2_idx = int(parts[1])
            a1 = alleles[a1_idx] if a1_idx < len(alleles) else None
            a2 = alleles[a2_idx] if a2_idx < len(alleles) else None
            return (a1, a2)
        except (ValueError, IndexError):
            return (None, None)

    def index_positions(
        self,
        positions: set,
        target_samples: Optional[set] = None,
        missing_rs_ids: Optional[set] = None,
        mode: str = "auto",
    ) -> None:
        normalized_positions = {(str(c).replace("chr", ""), str(p)) for c, p in positions}
        missing_rs_ids = set(missing_rs_ids) if missing_rs_ids else set()
        if not normalized_positions and not missing_rs_ids:
            self.logger.info("No positions or RS IDs requested; skipping imputed VCF indexing.")
            return

        mode = (mode or "auto").lower()
        allow_tabix = mode in ("auto", "tabix")
        allow_stream = mode in ("auto", "stream")

        self.parse_header_only()
        self.position_index = {}
        self.rsid_index = {}
        self.indexed_variant_count = 0
        self.tabix_hits = set()
        self.stream_hits = set()

        target_samples = set(target_samples) if target_samples else set(self.samples)
        sample_to_idx = {s: i for i, s in enumerate(self.samples)}

        tbi_path = Path(str(self.file_path) + ".tbi")
        csi_path = Path(str(self.file_path) + ".csi")
        has_index = tbi_path.exists() or csi_path.exists()

        used_tabix = False
        if allow_tabix and has_index and normalized_positions:
            try:
                import pysam  # type: ignore

                self.logger.info("Using tabix index for targeted fetching")
                vcf = pysam.VariantFile(str(self.file_path))
                for chrom, pos in normalized_positions:
                    for query_chrom in (chrom, f"chr{chrom}"):
                        try:
                            for rec in vcf.fetch(query_chrom, int(pos) - 1, int(pos)):
                                if str(rec.pos) != pos:
                                    continue
                                ref = rec.ref
                                alt = rec.alts[0] if rec.alts else "."
                                sample_gts = {}
                                for sample in target_samples:
                                    if sample in rec.samples:
                                        gt_tuple = rec.samples[sample].get("GT")
                                        if gt_tuple:
                                            gt_str = "/".join("." if a is None else str(a) for a in gt_tuple)
                                        else:
                                            gt_str = "./."
                                        sample_gts[sample] = gt_str
                                self.position_index[(chrom, pos)] = {
                                    "REF": ref,
                                    "ALT": alt,
                                    "samples": sample_gts,
                                }
                                self.indexed_variant_count += 1
                                self.tabix_hits.add((chrom, pos))
                                break
                            else:
                                continue
                            break
                        except ValueError:
                            continue
                used_tabix = True
            except ImportError:
                self.logger.warning("pysam not installed; falling back to streaming read.")
            except Exception as e:
                self.logger.warning(f"Tabix fetch failed ({e}); falling back to streaming read.")

        import gzip

        open_func = gzip.open if str(self.file_path).endswith(".gz") else open
        mode = "rt" if str(self.file_path).endswith(".gz") else "r"

        need_stream = allow_stream
        if allow_stream:
            if used_tabix and not missing_rs_ids and len(self.position_index) >= len(normalized_positions):
                need_stream = False
        else:
            need_stream = False

        if not need_stream:
            return

        self.logger.info(
            f"Streaming VCF to find positions ({len(normalized_positions)}) "
            f"and RS IDs ({len(missing_rs_ids)})..."
        )
        with open_func(self.file_path, mode, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                raw_chrom = parts[0]
                chrom = raw_chrom.replace("chr", "")
                pos = parts[1]
                rsid = parts[2] if len(parts) > 2 else ""
                key = (chrom, pos)
                hit_position = key in normalized_positions
                hit_rsid = rsid in missing_rs_ids
                if not hit_position and not hit_rsid:
                    continue

                ref = parts[3]
                alt = parts[4].split(",")[0] if parts[4] else "."
                sample_fields = parts[9:]

                sample_gts = {}
                for sample in target_samples:
                    idx = sample_to_idx.get(sample)
                    if idx is None or idx >= len(sample_fields):
                        continue
                    sample_value = sample_fields[idx]
                    gt_value = sample_value.split(":", 1)[0] if sample_value else "./."
                    sample_gts[sample] = gt_value

                if hit_position:
                    self.position_index[key] = {"REF": ref, "ALT": alt, "samples": sample_gts}
                    self.indexed_variant_count += 1
                    self.stream_hits.add(key)

                if hit_rsid:
                    self.rsid_index[rsid] = {
                        "CHROM": chrom,
                        "POS": pos,
                        "REF": ref,
                        "ALT": alt,
                        "samples": sample_gts,
                    }
                    self.indexed_variant_count += 1

    def extract_markers(
        self,
        marker_df: pd.DataFrame,
        sample_list: List[str],
    ) -> Tuple[Dict[str, Dict[str, Dict]], Dict[str, Tuple]]:
        """
        마커 리스트와 샘플 리스트에 해당하는 genotype 추출
        Returns:
            (results, marker_positions)
            results: {rs_id: {sample_id: {...}}}
            marker_positions: {rs_id: (chrom, pos, ref, alt)}
        """
        results: Dict[str, Dict[str, Dict]] = {}
        marker_positions: Dict[str, Tuple] = {}

        for _, row in marker_df.iterrows():
            rs_id = row["RS_ID"]
            variant = self.get_variant_by_rsid(rs_id)
            if variant is None:
                self.logger.warning(f"RS ID not found in VCF: {rs_id}")
                continue

            chrom = variant["CHROM"]
            pos = variant["POS"]
            ref = variant["REF"]
            alt = variant["ALT"]

            marker_positions[rs_id] = (chrom, pos, ref, alt)
            results[rs_id] = {}

            for sample in sample_list:
                if sample in self.samples:
                    gt = variant[sample]
                    allele1, allele2 = self.decode_genotype(gt, ref, alt)
                    results[rs_id][sample] = {
                        "allele1": allele1,
                        "allele2": allele2,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "gt_raw": gt,
                        "is_no_call": allele1 is None or allele2 is None,
                    }
                else:
                    self.logger.warning(f"Sample not found in VCF: {sample}")

        self.logger.info(f"Extracted {len(results)} markers for {len(sample_list)} samples")
        return results, marker_positions

