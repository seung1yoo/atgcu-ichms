import argparse
import sys

from .config import load_conf, setup_logger
from .pipeline_api import run_api
from .pipeline_extract import run_extract


def build_parser() -> argparse.ArgumentParser:
    # 공통 옵션 (subparser에도 적용)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--conf-json", type=str, default="conf.json", help="Settings JSON path (default: conf.json)")
    common.add_argument(
        "--log-level", type=str, choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO"
    )

    parser = argparse.ArgumentParser(
        description="atgcu-ichms marker parser",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[common],
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # extract
    p_extract = sub.add_parser("extract", help="Run extract pipeline (OA + JSON)", parents=[common])
    p_extract.add_argument("--product", type=str, required=True, help="Product key in conf.json")
    p_extract.add_argument("--sample-list", type=str, required=True, help="Sample list file path")
    p_extract.add_argument("--plate-barcode", type=str, required=True, help="Plate barcode")

    # api
    p_api = sub.add_parser("api-to-lis", help="Send JSON to LIS API", parents=[common])
    p_api.add_argument("--product", type=str, required=True, help="Product key in conf.json")
    p_api.add_argument("--json", type=str, required=True, help="JSON file to post")
    p_api.add_argument("--timeout", type=float, default=10.0, help="API timeout seconds")
    p_api.add_argument("--retries", type=int, default=3, help="API retries")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        conf = load_conf(args.conf_json, args.product)
    except Exception as e:
        print(f"[ERROR] Failed to load config: {e}", file=sys.stderr)
        return 1

    output_dir = conf.get("output_dir", ".")
    qc_report_path = None
    if args.command == "extract":
        from pathlib import Path
        qc_report_path = conf.get("qc_report") or str(Path(output_dir) / f"{args.product}.oa.txt.qc.report.txt")

    logger = setup_logger("atgcu-ichms", args.log_level, log_file=qc_report_path)

    if args.command == "extract":
        return run_extract(conf, args, logger)
    if args.command == "api-to-lis":
        return run_api(conf, args, logger)
    return 1


if __name__ == "__main__":
    sys.exit(main())

