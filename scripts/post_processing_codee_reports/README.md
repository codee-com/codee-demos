# Codee Report Analyzer

Analyze Codee reports over time and generate interactive HTML visualizations.

## Overview

This script analyzes multiple Codee reports, tracking how checker findings evolve over time.
It supports both JSON and HTML report formats and generates an interactive HTML dashboard with three
charts:

- **Total Findings Over Time**: Line chart showing total findings per report
- **Findings by Checker**: Line chart for each unique checker type found
- **Findings by Priority Level**: Line chart grouped by L1/L2/L3/L4 priority

## Requirements

- Python 3.10+
- pandas

## Installation

```bash
pip install pandas
```

## Usage

### Basic Usage

```bash
python analyze_codee_reports.py <input_directory> <output_directory>
```

### Input Formats

The script automatically detects the input format:

** Codee JSON Reports** (from `codee checks --json`):
- Searches recursively for `*.json` files
- Example directory structure:
  ```
  reports/
  ├── screening-1.json
  ├── screening-2.json
  └── subdir/
      └── screening-3.json
  ```

**Codee HTML Reports** (from `codee --html`):
- Searches recursively for `report.js` files
- When HTML reports are found, the generated dashboard includes links to open each original report
- Example directory structure:
  ```
  reports/
  ├── run-2024-01-01/
  │   ├── report.js
  │   └── ...
  └── run-2024-02-01/
      ├── report.js
      └── ...
  ```

## Dashboard Features

The generated HTML dashboard includes:

- **Interactive charts**: Hover for tooltips, click on points to open original reports
- **Summary cards**: Show report count, unique checkers found, and the total findings from the latest report
- **Links section**: Table with links to open each original HTML report
- **Responsive design**: Works on desktop and mobile
- **No external dependencies**: Uses Chart.js from CDN (internet required for first load)

## Understanding the Charts

### Total Findings
Shows the sum of all checker findings for each report over time. Click on any point to open the original report (if available).

### Findings by Checker
Displays each unique checker type (e.g., PWR007, PWR068) as a separate line, allowing you to see which specific checkers contribute to findings.

### Findings by Priority Level
Groups findings by their severity level:
- **L1** (Red): High priority
- **L2** (Orange): Medium priority
- **L3** (Green): Low priority
- **L4** (Light Green): Very-low priority

## Troubleshooting

### "No JSON/HTML files found"
- Ensure your input directory contains valid report files
- For JSON: look for `*.json` files
- For HTML: look for directories containing `report.js`

### "Malformed JSON/JS"
- One or more files are corrupted
- The script will skip these files and continue

### "No valid timestamp"
- Some files may be missing the timestamp field
- The script will skip these files and continue
