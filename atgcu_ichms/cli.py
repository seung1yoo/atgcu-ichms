import argparse
import sys
from pathlib import Path

from .config import load_conf, setup_logger
from .pipeline_api import run_api
from .pipeline_extract import run_extract

PRODUCT_CHOICES = ["irs_v2", "dtc_53", "dtc_62"]

def build_parser() -> argparse.ArgumentParser:
    # 공통 옵션 (subparser에도 적용)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "--log-level", 
        type=str, 
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
        default="INFO",
    )
    common.add_argument(
        "--conf-json",
        type=str, 
        default="conf.json", 
        help="Settings JSON path (default: conf.json)",
    )

    parser = argparse.ArgumentParser(
        description="atgcu-ichms marker parser",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[common],
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # extract
    p_extract = sub.add_parser("extract", help="Run extract pipeline (OA + JSON)", parents=[common])
    p_extract.add_argument(
        "--product", 
        type=str, 
        required=True, 
        help="Product key in conf.json",
        choices=PRODUCT_CHOICES
    )
    p_extract.add_argument(
        "--imputed-mode",
        type=str,
        choices=["auto", "tabix", "stream"],
        default="auto",
        help="Imputed fetch mode: auto (tabix then stream), tabix-only, or stream-only (default: auto)",
    )
    p_extract.add_argument(
        "--path-anl-for-all",
        type=str,
        default=None,
        help="Base path for Analysis_For_All; apt_vcf/imputed_vcf will be joined to this path",
        required=True
    )
    p_extract.add_argument(
        "--path-wkdir",
        type=str,
        default=None,
        help="Base path for working dir; output_dir and qc_report will be created under this path",
        required=True
    )
    p_extract.add_argument(
        "--sample-sheet",
        type=str,
        required=True,
        help="Sample sheet TSV path"
    )

    # api
    p_api = sub.add_parser("api-to-lis", help="Send JSON to LIS API", parents=[common])
    p_api.add_argument(
        "--product",
        type=str,
        required=True,
        help="Product key in conf.json",
        choices=PRODUCT_CHOICES
    )
    p_api.add_argument(
        "--json",
        type=str,
        required=True,
        help="JSON file to post"
    )
    p_api.add_argument(
        "--timeout",
        type=float,
        default=10.0,
        help="API timeout seconds")
    p_api.add_argument(
        "--retries",
        type=int,
        default=3,
        help="API retries (default: 3)",
    )
    p_api.add_argument(
        "--insecure",
        action="store_true",
        help="Disable TLS certificate verification for API POST (not recommended)",
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        conf = load_conf(args.conf_json, args.product)
    except Exception as e:
        print(f"[ERROR] Failed to load config: {e}", file=sys.stderr)
        return 1

    def _join_if_relative(p: str, base: str) -> str:
        path_obj = Path(p)
        return str(path_obj) if path_obj.is_absolute() else str(Path(base) / path_obj)

    if args.command == "extract":
        base_path = args.path_anl_for_all
        if base_path:
            conf["apt_vcf"] = _join_if_relative(conf["apt_vcf"], base_path)
            conf["imputed_vcf"] = _join_if_relative(conf["imputed_vcf"], base_path)
            conf["genocare_all_gt"] = _join_if_relative(conf["genocare_all_gt"], base_path)

        work_base = args.path_wkdir
        if work_base:
            conf["output_dir"] = _join_if_relative(conf.get("output_dir", "."), work_base)
            if conf.get("qc_report"):
                conf["qc_report"] = _join_if_relative(conf["qc_report"], work_base)

        output_dir = conf.get("output_dir", ".")
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        qc_report_path = conf.get("qc_report") or str(Path(output_dir) / f"{args.product}.oa.txt.qc.report.txt")
        logger = setup_logger("atgcu-ichms", args.log_level, log_file=qc_report_path)
        return run_extract(conf, args, logger)
    if args.command == "api-to-lis":
        logger = setup_logger("atgcu-ichms", args.log_level)
        return run_api(conf, args, logger)
    return 1


if __name__ == "__main__":
    sys.exit(main())

