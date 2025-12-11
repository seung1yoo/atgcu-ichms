import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import pandas as pd

from .io_utils import sanitize_allele, load_oa_table


class QCReporter:
    def __init__(self, logger: logging.Logger = None):
        self.logger = logger or logging.getLogger(__name__)
        self.report_lines: List[str] = []

    def add_section(self, title: str) -> None:
        self.report_lines.append("")
        self.report_lines.append("=" * 70)
        self.report_lines.append(f"  {title}")
        self.report_lines.append("=" * 70)

    def add_line(self, text: str) -> None:
        self.report_lines.append(text)

    def add_key_value(self, key: str, value) -> None:
        self.report_lines.append(f"  - {key}: {value}")

    def analyze_list_file(self, items: List[str], list_name: str) -> Dict[str, any]:
        total_count = len(items)
        unique_items = list(dict.fromkeys(items))
        unique_count = len(unique_items)
        duplicates = [item for item in set(items) if items.count(item) > 1]

        self.add_section(f"{list_name} 분석")
        self.add_key_value("총 항목 수", total_count)
        self.add_key_value("중복 제거 후 항목 수", unique_count)
        self.add_key_value("중복 항목 수", total_count - unique_count)

        if duplicates:
            self.add_line(f"  - 중복된 항목들: {duplicates}")

        return {
            "total": total_count,
            "unique": unique_count,
            "duplicates": duplicates,
            "unique_items": unique_items,
        }

    def get_report(self) -> str:
        return "\n".join(self.report_lines)

    def write_report(self, output_path: str) -> None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "a", encoding="utf-8") as f:
            f.write(self.get_report())
            f.write("\n")
        self.logger.info(f"QC report written to: {output_path}")

    def print_report(self) -> None:
        print(self.get_report())



#---- helper functions ----

def add_run_args(qc: QCReporter, product: str, args) -> None:
    qc.add_section(f"{product} Markers Parser QC Report")
    qc.add_line(f"  실행 시간: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    qc.add_line("")
    qc.add_line("  [실행 인자]")
    qc.add_line(f"    --product       : {getattr(args, 'product', '')}")
    qc.add_line(f"    --sample-sheet  : {getattr(args, 'sample_sheet', '')}")
    qc.add_line(f"    --path-anl-for-all: {getattr(args, 'path_anl_for_all', '')}")
    qc.add_line(f"    --path-wkdir    : {getattr(args, 'path_wkdir', '')}")
    qc.add_line(f"    --imputed-mode  : {getattr(args, 'imputed_mode', 'auto')}")
    qc.add_line("")


def add_marker_info(qc: QCReporter, marker_df: pd.DataFrame) -> None:
    qc.add_section("RS Marker TSV 정보")
    qc.add_key_value("총 마커 수", len(marker_df))
    qc.add_line("")
    qc.add_line("  [마커별 Allele 정보 (처음 10개)]")
    for _, row in marker_df.head(10).iterrows():
        qc.add_line(f"    - {row['RS_ID']}: NORMAL={row['NORMAL_ALLELE']}, RISK={row['RISK_ALLELE']}")
    if len(marker_df) > 10:
        qc.add_line(f"    ... 외 {len(marker_df) - 10}개")


def init_rs_rows(marker_list: List[str]) -> Dict[str, Dict]:
    return {rs_id: {"RS_ID": rs_id} for rs_id in marker_list}


def update_rs_with_apt(qc_rs_rows: Dict[str, Dict], marker_list: List[str], apt_results, marker_positions):
    for rs_id in marker_list:
        rs_entry = qc_rs_rows[rs_id]
        rs_entry["apt_found"] = rs_id in apt_results
        if rs_id in marker_positions:
            chrom, pos, ref, alt = marker_positions[rs_id]
            rs_entry["apt_chrom"] = chrom
            rs_entry["apt_pos"] = pos
            rs_entry["apt_ref"] = ref
            rs_entry["apt_alt"] = alt
        else:
            rs_entry["apt_chrom"] = "."
            rs_entry["apt_pos"] = "."
            rs_entry["apt_ref"] = "."
            rs_entry["apt_alt"] = "."


def add_imputed_section(qc: QCReporter, imputed_vcf: str, imputed_parser, total_positions_len: int, imputed_mode: str):
    qc.add_section("Imputed VCF 파싱")
    qc.add_key_value("파일", imputed_vcf)
    qc.add_key_value("Imputed 모드", imputed_mode)
    qc.add_key_value("VCF 내 샘플 수", len(imputed_parser.get_samples()))
    qc.add_key_value("Imputed 조회 위치 수", total_positions_len)


def add_imputed_hits_detail(qc: QCReporter, imputed_parser) -> None:
    qc.add_key_value("조회된 variant 수", imputed_parser.indexed_variant_count)
    qc.add_key_value("tabix 조회 성공 위치 수", len(imputed_parser.tabix_hits))
    qc.add_key_value("스트리밍 조회 성공 위치 수", len(imputed_parser.stream_hits))


def write_rs_qc(output_dir: Path, product: str, qc_rs_rows: Dict[str, Dict], logger: logging.Logger):
    qc_rs_df = pd.DataFrame(list(qc_rs_rows.values()))
    qc_rs_path = output_dir / f"{product}.qc_rs.tsv"
    qc_rs_df.to_csv(qc_rs_path, sep="\t", index=False)
    logger.info(f"Written RS-level QC to {qc_rs_path}")
    return qc_rs_df, qc_rs_path


def compare_with_genocare(
    apt_results,
    sample_list: List[str],
    marker_list: List[str],
    sample_id_map: Dict[str, str],
    genocare_path: Path,
    logger: logging.Logger,
):
    geno_df = load_oa_table(str(genocare_path), logger)

    def _canon(a1, a2):
        a1 = sanitize_allele(a1)
        a2 = sanitize_allele(a2)
        return "".join(sorted([a1, a2]))

    genocare_map = {}
    required_cols = {"Sample ID", "NCBI SNP Reference", "Allele 1 Call", "Allele 2 Call"}
    if required_cols.issubset(set(geno_df.columns)):
        for _, row in geno_df.iterrows():
            sid = str(row["Sample ID"])
            rs = str(row["NCBI SNP Reference"])
            genocare_map[(sid, rs)] = _canon(row["Allele 1 Call"], row["Allele 2 Call"])

    compare_rows = []
    total = match = mismatch = missing_genocare = 0
    for rs_id in marker_list:
        for sample_name in sample_list:
            ours_sid = sample_name  # SAMPLE_NAME 기준
            if rs_id not in apt_results or sample_name not in apt_results[rs_id]:
                continue
            ours_info = apt_results[rs_id][sample_name]
            ours_gt = _canon(ours_info["allele1"], ours_info["allele2"])
            geno_gt = genocare_map.get((ours_sid, rs_id))
            status = ""
            if geno_gt is None:
                missing_genocare += 1
                status = "missing_genocare"
            else:
                total += 1
                if geno_gt == ours_gt:
                    match += 1
                    status = "match"
                else:
                    mismatch += 1
                    status = "mismatch"
            compare_rows.append(
                {
                    "Sample ID": ours_sid,
                    "RS_ID": rs_id,
                    "ours_gt": ours_gt,
                    "genocare_gt": geno_gt if geno_gt is not None else ".",
                    "status": status,
                }
            )

    compare_df = pd.DataFrame(compare_rows)
    summary = {
        "total": total,
        "match": match,
        "mismatch": mismatch,
        "missing_genocare": missing_genocare,
    }
    return compare_df, summary

