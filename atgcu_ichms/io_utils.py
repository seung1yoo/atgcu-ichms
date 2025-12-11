import json
import logging
from pathlib import Path
from typing import Optional, Dict

import pandas as pd


def load_sample_sheet(file_path: str, logger: logging.Logger = None) -> Dict[str, any]:
    """
    Load sample sheet TSV with columns:
    - SAMPLE_NAME
    - PLATE_BARCODE
    - SAMPLE_ID_IN_A_PLATE
    Returns dict with sample_names (list), plate_barcode (str), and sample_id_map (dict).
    """
    logger = logger or logging.getLogger(__name__)
    df = pd.read_csv(file_path, sep="\t", dtype=str)
    required_cols = ["SAMPLE_NAME", "PLATE_BARCODE", "SAMPLE_ID_IN_A_PLATE"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in sample sheet: {missing_cols}")

    df = (
        df[required_cols]
        .fillna("")
        .map(lambda x: str(x).strip())
    )

    sample_names = df["SAMPLE_NAME"].tolist()
    plate_barcodes = df["PLATE_BARCODE"].unique().tolist()
    sample_id_map = dict(zip(df["SAMPLE_NAME"], df["SAMPLE_ID_IN_A_PLATE"]))

    # Plate barcode should be single value to match previous CLI semantics
    plate_barcodes = [p for p in plate_barcodes if p != ""]
    if not plate_barcodes:
        raise ValueError("PLATE_BARCODE is empty in sample sheet.")
    if len(set(plate_barcodes)) > 1:
        raise ValueError(f"Multiple PLATE_BARCODE values found: {plate_barcodes}")
    plate_barcode = plate_barcodes[0]

    if len(sample_id_map) != len(sample_names):
        dupes = [s for s in sample_names if sample_names.count(s) > 1]
        raise ValueError(f"Duplicate SAMPLE_NAME detected: {sorted(set(dupes))}")

    logger.info(f"Loaded {len(sample_names)} samples from {file_path} (plate: {plate_barcode})")
    return {
        "sample_names": sample_names,
        "plate_barcode": plate_barcode,
        "sample_id_map": sample_id_map,
    }


def load_marker_tsv(file_path: str, logger: logging.Logger = None) -> pd.DataFrame:
    logger = logger or logging.getLogger(__name__)
    df = pd.read_csv(file_path, sep="\t", encoding="utf-8")
    required_cols = ["RS_ID", "NORMAL_ALLELE", "RISK_ALLELE"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in marker TSV: {missing_cols}")
    logger.info(f"Loaded {len(df)} markers from {file_path}")
    return df


def sanitize_allele(allele: Optional[str]) -> str:
    if allele is None:
        return "."
    a = str(allele)
    if a in ["N", "*", "-"]:
        return "."
    return a


def load_oa_table(file_path: str, logger: logging.Logger = None) -> pd.DataFrame:
    logger = logger or logging.getLogger(__name__)
    header = None
    rows = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(parts[: len(header)])
    if header is None:
        raise ValueError(f"OA table header not found in {file_path}")
    df = pd.DataFrame(rows, columns=header)
    logger.info(f"Loaded OA table: {len(df)} rows, {len(df.columns)} cols from {file_path}")
    return df


def write_oa(output_path: Path, meta: dict, data: pd.DataFrame, logger: logging.Logger) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        for key, value in meta.items():
            f.write(f"# {key} : {value}\n")
        f.write("\n")
        data.to_csv(f, sep="\t", index=False)
    logger.info(f"Written {len(data)} rows to {output_path}")


def write_json(output_path: Path, rows: list, logger: logging.Logger) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as jf:
        json.dump(rows, jf, ensure_ascii=False, indent=2)
    logger.info(f"Written {len(rows)} JSON records to {output_path}")

