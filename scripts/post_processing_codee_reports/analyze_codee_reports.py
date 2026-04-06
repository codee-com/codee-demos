#!/usr/bin/env python3
"""Analyze Codee JSON or HTML reports over time and generate time-series visualizations."""

import argparse
import json
import logging
import re
import sys
import pandas as pd
from collections import defaultdict
from datetime import datetime
from pathlib import Path


def setup_logging() -> logging.Logger:
    """Configure logging for the script."""
    logger = logging.getLogger("codee_analyzer")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(handler)
    return logger


def load_json_file(file_path: Path, logger: logging.Logger) -> dict | None:
    """Load and parse a JSON file, returning None on failure."""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        logger.error(f"Malformed JSON in {file_path.name}: {e}")
        return None
    except OSError as e:
        logger.error(f"Failed to read {file_path.name}: {e}")
        return None


def load_html_report(report_dir: Path, logger: logging.Logger) -> dict | None:
    """Load a Codee HTML report by parsing report.js file."""
    report_js_path = report_dir / "report.js"
    if not report_js_path.exists():
        logger.error(f"report.js not found in {report_dir.name}")
        return None
    try:
        # Extract JavaScript object from report.js
        content = report_js_path.read_text(encoding="utf-8")
        match = re.search(r"const\s+report\s*=\s*(\{.*\});", content, re.DOTALL)
        if not match:
            match = re.search(r"var\s+report\s*=\s*(\{.*\});", content, re.DOTALL)

        if not match:
            logger.error(f"Could not find report object in {report_js_path}")
            return None

        json_str = match.group(1)
        return json.loads(json_str)
    except Exception as e:
        logger.error(f"Failed to parse report.js in {report_dir.name}: {e}")
        return None


def extract_timestamp(report: dict, logger: logging.Logger) -> int | None:
    """Extract timestamp from report."""
    try:
        exec_info = report.get("CodeeExecutionInfo", {})
        ts_utc = exec_info.get("TimestampUTC")
        if ts_utc:
            dt = datetime.fromisoformat(ts_utc.replace("Z", "+00:00"))
            return int(dt.timestamp())
        return exec_info.get("TimestampEpochSeconds")
    except Exception as e:
        logger.error(f"Failed to extract timestamp: {e}")
        return None


def _extract_table_values(
    report: dict,
    table_key: str,
    value_key: str,
    logger: logging.Logger,
    skip_empty: bool = False,
) -> dict[str, str]:
    """Extract values from Quality and Optimization checker tables."""
    result: dict[str, str] = {}
    try:
        screening = report.get("Screening", {})
        for table_name in (
            "Ranking of Quality Checkers",
            "Ranking of Optimization Checkers",
        ):
            table = screening.get(table_name, {}).get("DataTable", [])
            for row in table:
                checker = row.get("Checker", "")
                if checker and checker != "Total":
                    value = row.get(value_key, "")
                    if skip_empty and not value:
                        continue
                    result[checker] = value
    except Exception as e:
        logger.warning(f"Failed to extract {value_key} from {table_key}: {e}")
    return result


def extract_checker_counts(report: dict, logger: logging.Logger) -> dict[str, int]:
    """Extract checker counts from Quality and Optimization rankings."""
    raw_counts = _extract_table_values(report, "checker counts", "#", logger)
    counts: dict[str, int] = {}
    for checker, count_str in raw_counts.items():
        try:
            counts[checker] = int(count_str)
        except ValueError:
            logger.error(f"Invalid count '{count_str}' for checker {checker}")
    return counts


def extract_checker_priorities(report: dict, logger: logging.Logger) -> dict[str, str]:
    """Extract priority for each checker."""
    return _extract_table_values(
        report, "checker priorities", "Priority", logger, skip_empty=True
    )


def extract_l_level(priority_str: str) -> str:
    """Extract L-level from priority string like 'P18 (L1)'."""
    if not priority_str:
        return "Unknown"
    match priority_str:
        case s if "L1" in s:
            return "L1"
        case s if "L2" in s:
            return "L2"
        case s if "L3" in s:
            return "L3"
        case s if "L4" in s:
            return "L4"
        case _:
            return "Unknown"


def detect_input_type(input_dir: Path) -> str:
    """Detect if input is JSON or HTML format by searching recursively."""
    has_json = list(input_dir.rglob("*.json"))
    has_html_dirs = list(input_dir.rglob("report.js"))

    if has_html_dirs:
        return "html"
    elif has_json:
        return "json"
    return "unknown"


def load_reports(
    input_dir: Path, logger: logging.Logger
) -> list[tuple[int, dict, Path | None]]:
    """Load all reports from directory recursively, sorted by timestamp.

    Returns list of tuples: (timestamp, report_dict, path)
    path is the path to the report file (for linking)
    """
    reports: list[tuple[int, dict, Path | None]] = []
    input_type = detect_input_type(input_dir)

    if input_type == "html":
        logger.info("Detected HTML format")
        html_dirs = sorted(input_dir.rglob("report.js"))
        for report_js in html_dirs:
            report_dir = report_js.parent
            report = load_html_report(report_dir, logger)
            if report is None:
                continue

            timestamp = extract_timestamp(report, logger)
            if timestamp is None:
                logger.warning(f"Skipping {report_dir.name}: no valid timestamp")
                continue

            reports.append((timestamp, report, report_dir))

    elif input_type == "json":
        logger.info("Detected JSON format")
        json_files = sorted(input_dir.rglob("*.json"))
        if not json_files:
            logger.warning(f"No JSON files found in {input_dir}")
            return reports

        for file_path in json_files:
            report = load_json_file(file_path, logger)
            if report is None:
                continue

            timestamp = extract_timestamp(report, logger)
            if timestamp is None:
                logger.warning(f"Skipping {file_path.name}: no valid timestamp")
                continue

            reports.append((timestamp, report, file_path))

    else:
        logger.error(f"No valid reports found in {input_dir}")
        return reports

    reports.sort(key=lambda x: x[0])
    return reports


def build_checker_dataframe(
    reports: list[tuple[int, dict, Path | None]], logger: logging.Logger
) -> pd.DataFrame:
    """Build a DataFrame with timestamps as index and checker counts as columns."""
    all_checkers: set[str] = set()

    for _, report, _ in reports:
        counts = extract_checker_counts(report, logger)
        all_checkers.update(counts.keys())

    sorted_checkers = sorted(all_checkers)
    data: dict[str, list[int]] = {checker: [] for checker in sorted_checkers}
    timestamps: list[int] = []

    for timestamp, report, _ in reports:
        timestamps.append(timestamp)
        counts = extract_checker_counts(report, logger)
        for checker in sorted_checkers:
            data[checker].append(counts.get(checker, 0))

    df = pd.DataFrame(data, index=pd.to_datetime(timestamps, unit="s"))
    df.index.name = "timestamp"
    return df


def get_checker_priorities_for_df(
    df: pd.DataFrame,
    reports: list[tuple[int, dict, Path | None]],
    logger: logging.Logger,
) -> dict[str, str]:
    """Get priority mapping for all checkers in the DataFrame."""
    priorities: dict[str, str] = {}

    for _, report, _ in reports:
        report_priorities = extract_checker_priorities(report, logger)
        priorities.update(report_priorities)

    return priorities


def generate_html_report(
    df: pd.DataFrame,
    reports: list[tuple[int, dict, Path | None]],
    input_dir: Path,
    output_dir: Path,
    logger: logging.Logger,
) -> None:
    """Generate interactive HTML report with Chart.js."""
    logger.info("Generating interactive HTML report")

    priorities = get_checker_priorities_for_df(df, reports, logger)
    labels = [d.strftime("%Y-%m-%d") for d in df.index]
    total_data = df.sum(axis=1).tolist()

    checker_data: dict[str, list[int]] = {}
    for checker in df.columns:
        checker_data[checker] = df[checker].tolist()

    priority_order = ["L1", "L2", "L3", "L4", "Unknown"]
    priority_colors = {
        "L1": "#DC143C",
        "L2": "#FF8C00",
        "L3": "#228B22",
        "L4": "#90EE90",
        "Unknown": "#D3D3D3",
    }

    # l_level_data: dict[str, list[int]] = {"L1": [5, 2, ...], "L2": [3, 1, ...], ...}
    # p_groups: dict[str, int] = {"L1": 5, "L2": 3, ...}
    # per-row counts by priority level
    l_level_data: dict[str, list[int]] = {l: [] for l in priority_order}
    for idx in df.index:
        row = df.loc[idx]
        p_groups: dict[str, int] = {l: 0 for l in priority_order}
        for checker, count in row.items():
            l_level = extract_l_level(priorities.get(checker, ""))
            p_groups[l_level] += count
        for l_level in priority_order:
            l_level_data[l_level].append(p_groups[l_level])

    # Build links to original reports (JSON or HTML)
    report_links: list[dict | None] = []
    for _, _, report_path in reports:
        if report_path:
            abs_path = report_path.resolve()
            if report_path.suffix == ".json":
                report_links.append({"link": str(abs_path), "type": "json"})
            else:
                index_html = abs_path / "index.html"
                if index_html.exists():
                    report_links.append({"link": str(index_html), "type": "html"})
                else:
                    report_links.append({"link": str(abs_path), "type": "html"})
        else:
            report_links.append(None)

    chart_data = {
        "labels": labels,
        "total": total_data,
        "checkers": checker_data,
        "priorities": l_level_data,
        "priorityColors": priority_colors,
        "reportLinks": report_links,
        "hasLinks": any(link is not None for link in report_links),
    }

    template_path = Path(__file__).parent / "templates" / "codee_report.html"
    css_path = Path(__file__).parent / "templates" / "codee_report.css"
    js_path = Path(__file__).parent / "templates" / "codee_report.js"
    html_template = template_path.read_text(encoding="utf-8")
    css_content = css_path.read_text(encoding="utf-8")
    js_template = js_path.read_text(encoding="utf-8")

    replacements = {
        "${REPORT_COUNT}": str(len(reports)),
        "${CHECKER_COUNT}": str(len(checker_data)),
        "${TOTAL_FINDINGS}": str(int(total_data[-1])),
        "${INPUT_DIR}": str(input_dir),
        "${CHART_DATA}": json.dumps(chart_data, indent=2),
    }

    for placeholder, value in replacements.items():
        html_template = html_template.replace(placeholder, value)
        js_template = js_template.replace(placeholder, value)

    html_path = output_dir / "codee_analysis.html"
    css_output_path = output_dir / "codee_report.css"
    js_output_path = output_dir / "codee_report.js"
    html_path.write_text(html_template, encoding="utf-8")
    css_output_path.write_text(css_content, encoding="utf-8")
    js_output_path.write_text(js_template, encoding="utf-8")
    logger.info(f"HTML report saved to {html_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse Codee JSON or HTML reports and generate an HTML summary visualization."
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing Codee report files (JSON or HTML with report.js)",
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory where the HTML report will be saved",
    )
    args = parser.parse_args()

    logger = setup_logging()

    if not args.input_dir.is_dir():
        logger.error(f"Input directory does not exist: {args.input_dir}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading reports from {args.input_dir}")
    reports = load_reports(args.input_dir, logger)

    if not reports:
        logger.error("No valid reports loaded. Exiting.")
        return 0

    logger.info(f"Loaded {len(reports)} reports")

    logger.info("Building checker DataFrame")
    df = build_checker_dataframe(reports, logger)

    generate_html_report(df, reports, args.input_dir, args.output_dir, logger)

    return 0


if __name__ == "__main__":
    sys.exit(main())
