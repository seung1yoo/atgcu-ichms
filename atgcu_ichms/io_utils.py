import json
import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd


def load_list_file(file_path: str, logger: logging.Logger = None) -> List[str]:
    logger = logger or logging.getLogger(__name__)
    items = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                items.append(line)
    logger.info(f"Loaded {len(items)} items from {file_path}")
    return items


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

