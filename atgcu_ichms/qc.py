import logging
from pathlib import Path
from typing import List, Dict
import pandas as pd


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

