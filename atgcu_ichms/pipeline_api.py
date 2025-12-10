import json
import logging
from pathlib import Path
from typing import Dict, Any

import requests


def run_api(conf: Dict[str, Any], args, logger: logging.Logger) -> int:
    api_conf = conf.get("api", {})
    endpoint = api_conf.get("endpoint")
    if not endpoint:
        logger.error("API endpoint not configured in conf.json (product.api.endpoint).")
        return 1

    headers = api_conf.get("headers", {}) or {}
    token = api_conf.get("auth_token")
    if token:
        headers.setdefault("Authorization", f"Bearer {token}")
    headers.setdefault("Content-Type", "application/json")

    json_path = Path(args.json)
    if not json_path.exists():
        logger.error(f"JSON file not found: {json_path}")
        return 1

    with open(json_path, "r", encoding="utf-8") as f:
        payload = json.load(f)

    timeout = args.timeout
    retries = args.retries
    for attempt in range(1, retries + 1):
        try:
            resp = requests.post(endpoint, headers=headers, json=payload, timeout=timeout)
            if resp.status_code >= 200 and resp.status_code < 300:
                logger.info(f"API POST success: status={resp.status_code}")
                return 0
            logger.warning(f"API POST failed (attempt {attempt}/{retries}): status={resp.status_code}, body={resp.text[:500]}")
        except Exception as e:
            logger.warning(f"API POST error (attempt {attempt}/{retries}): {e}")
    logger.error("API POST failed after retries.")
    return 1

