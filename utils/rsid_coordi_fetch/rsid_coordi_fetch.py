#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RS 좌표 보충 스크립트
 - 입력: RS_ID, NORMAL_ALLELE, RISK_ALLELE, CHROM, POS 컬럼을 가진 TSV
 - CHROM 또는 POS가 비어 있는 RS에 대해 dbSNP API에서 좌표를 조회해 채움
 - 결과는 별도 출력 파일로 저장
 - 기준 assembly는 --assembly 옵션으로 선택 (default: hg38)
"""

import argparse
import sys
import time
from typing import Optional, Tuple

import pandas as pd
import requests


def _normalize_assembly(assembly: str) -> str:
    """
    Normalize assembly keyword to GRCh37 / GRCh38.
    Accepts hg19/hg38 or GRCh37/GRCh38 (case-insensitive).
    """
    a = assembly.strip().lower()
    if a in ["hg38", "grch38"]:
        return "GRCh38"
    if a in ["hg19", "grch37"]:
        return "GRCh37"
    # fallback: return upper for visibility
    return assembly.upper()


def _normalize_chrom(chrom: str, add_chr: bool = True) -> Optional[str]:
    """
    Normalize chromosome: strip chr/NC_/leading zeros; optionally add chr prefix.
    """
    if chrom is None:
        return None
    c = str(chrom)
    c = c.replace("chr", "").replace("CHR", "")
    c = c.replace("NC_", "")
    c = c.split(".")[0]
    c = c.lstrip("0")
    c = c if c else "0"
    if add_chr and not c.startswith("chr"):
        c = "chr" + c
    return c


def fetch_coords(
    rsid: str,
    assembly: str = "GRCh38",
    timeout: float = 10.0,
    retries: int = 1,
    verbose: bool = False,
) -> Tuple[Optional[str], Optional[int]]:
    """
    NCBI dbSNP API에서 rsID 좌표 조회 (GRCh37 또는 GRCh38)
    Returns: (chrom, pos) 또는 (None, None)
    """
    target_assembly = _normalize_assembly(assembly)
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid.lstrip('rs')}"

    last_err = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(
                url,
                timeout=timeout,
                headers={
                    "User-Agent": "rsid-coord-fetch/0.1 (contact: anonymous@example.com)",
                    "Accept": "application/json",
                },
            )
            if r.status_code != 200:
                if verbose:
                    print(f"[WARN] {rsid} status={r.status_code} body={r.text[:200]}", file=sys.stderr)
                last_err = f"status {r.status_code}"
                time.sleep(0.1)
                continue
            try:
                data = r.json()
            except Exception as e:
                last_err = f"json decode error: {e}"
                if verbose:
                    print(f"[WARN] {rsid} json decode error: {e}", file=sys.stderr)
                time.sleep(0.1)
                continue
            break
        except Exception as e:
            last_err = str(e)
            if verbose:
                print(f"[WARN] {rsid} request error: {e}", file=sys.stderr)
            time.sleep(0.1)
            continue
    else:
        if verbose and last_err:
            print(f"[WARN] {rsid} failed after {retries} retries ({last_err})", file=sys.stderr)
        return None, None

    for placement in data.get("primary_snapshot_data", {}).get("placements_with_allele", []):
        # assembly_name가 placement_annot.assembly_name에 없고,
        # placement_annot.seq_id_traits_by_assembly[*].assembly_name에 들어있는 경우가 많다.
        pa = placement.get("placement_annot", {}) or {}
        assembly_name = pa.get("assembly_name", "")
        trait_names = [
            trait.get("assembly_name", "")
            for trait in pa.get("seq_id_traits_by_assembly", []) or []
        ]
        names_to_check = [assembly_name] + trait_names
        names_to_check = [n for n in names_to_check if n]

        if not any(n.startswith(target_assembly) for n in names_to_check):
            continue

        # placement["alleles"] 안의 첫 번째 allele 위치 사용
        alleles = placement.get("alleles", [])
        if not alleles:
            continue
        loc = alleles[0].get("allele", {}).get("spdi", {})
        chrom = placement.get("seq_id", None)
        pos0 = loc.get("position", None)
        if chrom is None or pos0 is None:
            continue
        norm_chrom = _normalize_chrom(chrom, add_chr=True)
        return norm_chrom, pos0 + 1  # 1-based로 변환

    if verbose:
        placements = data.get("primary_snapshot_data", {}).get("placements_with_allele", [])
        names = []
        for placement in placements:
            pa = placement.get("placement_annot", {}) or {}
            names.append(pa.get("assembly_name", ""))
            names.extend(
                trait.get("assembly_name", "")
                for trait in pa.get("seq_id_traits_by_assembly", []) or []
            )
        names = [n for n in names if n]
        names_preview = ", ".join(names[:5]) + ("..." if len(names) > 5 else "")
        print(
            f"[WARN] {rsid} no placement matched target assembly={target_assembly} "
            f"(placements={len(placements)}, assembly_names={names_preview})",
            file=sys.stderr,
        )
    return None, None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fill missing CHROM/POS in marker TSV using dbSNP API",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", "-i", required=True, help="Input marker TSV (RS_ID, NORMAL_ALLELE, RISK_ALLELE, CHROM, POS)")
    parser.add_argument("--output", "-o", required=True, help="Output TSV path (will be created)")
    parser.add_argument("--assembly", choices=["GRCh37", "GRCh38", "hg19", "hg38"], default="GRCh38", help="Reference assembly for coordinates")
    parser.add_argument("--sleep", type=float, default=0.2, help="Sleep seconds between API calls (rate-limit)")
    parser.add_argument("--timeout", type=float, default=10.0, help="HTTP request timeout (seconds)")
    parser.add_argument("--retries", type=int, default=3, help="HTTP request retries per rsID")
    parser.add_argument("--verbose", action="store_true", help="Print HTTP status/body when fetch fails")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    required_cols = ["RS_ID", "NORMAL_ALLELE", "RISK_ALLELE", "CHROM", "POS"]
    for col in required_cols:
        if col not in df.columns:
            print(f"[ERROR] Missing column: {col}", file=sys.stderr)
            return 1

    # 비어있는 대상 선택
    def _is_empty(val: str) -> bool:
        return val is None or val == "" or str(val).lower() == "nan"

    missing_mask = df["CHROM"].apply(_is_empty) | df["POS"].apply(_is_empty)
    missing_rs = df.loc[missing_mask, "RS_ID"].unique().tolist()

    print(f"총 {len(df)}개 RS 중 좌표가 비어있는 RS: {len(missing_rs)}개")

    fetched = 0
    for rsid in missing_rs:
        chrom, pos = fetch_coords(
            rsid,
            assembly=args.assembly,
            timeout=args.timeout,
            retries=args.retries,
            verbose=args.verbose,
        )
        if chrom is not None and pos is not None:
            df.loc[df["RS_ID"] == rsid, "CHROM"] = str(chrom).replace("chr", "")
            df.loc[df["RS_ID"] == rsid, "POS"] = str(pos)
            fetched += 1
            print(f"Filled {rsid}: {chrom}:{pos}")
        else:
            print(f"[WARN] Failed to fetch {rsid}")
        time.sleep(max(args.sleep, 0.0))

    df.to_csv(args.output, sep="\t", index=False)
    print(f"\n완료: {fetched}개 좌표 채움, 출력 → {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
