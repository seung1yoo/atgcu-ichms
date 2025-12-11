import json
import logging
import sys
from pathlib import Path
from typing import Dict


def load_conf(conf_path: str, product: str) -> Dict[str, any]:
    conf_path = Path(conf_path)
    if not conf_path.exists():
        raise FileNotFoundError(f"Config file not found: {conf_path}")

    with open(conf_path, "r", encoding="utf-8") as f:
        conf_obj = json.load(f)

    if product not in conf_obj:
        raise KeyError(f"Product '{product}' not found in config {conf_path}")

    conf = conf_obj[product]
    required_keys = ["apt_vcf", "imputed_vcf", "rsmarker_list", "output_dir"]
    missing = [k for k in required_keys if k not in conf or not conf[k]]
    if missing:
        raise KeyError(f"Missing required keys in product '{product}': {missing}")

    return conf


def setup_logger(name: str, log_level: str = "INFO", log_file: str = None) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, log_level.upper()))

    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)-8s %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, log_level.upper()))
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    if log_file and not any(
        isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", "") == str(Path(log_file))
        for h in logger.handlers
    ):
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)  # ensure log dir exists
        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setLevel(getattr(logging, log_level.upper()))
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger

