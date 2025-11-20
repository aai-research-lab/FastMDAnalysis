#!/usr/bin/env python3
"""Generate lines-of-code metrics and visualization for each workflow."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import benchmark_full_workflow as base_benchmark

DEFAULT_OUTPUT_DIR = Path("assets/benchmarks")
DEFAULT_JSON_NAME = "loc_metrics.json"
DEFAULT_IMAGE_NAME = "loc_comparison.png"


def write_loc_json(records: Sequence[dict], path: Path) -> None:
    payload = {
        "metadata": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "unit": "effective_loc",
            "notes": "Non-blank/non-comment LOC computed from workflow implementations.",
        },
        "records": list(records),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))
    print(f"Saved LOC metrics JSON -> {path}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate LOC metrics/visualization for workflow comparisons.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory to store LOC artefacts (default: assets/benchmarks).",
    )
    parser.add_argument(
        "--json-name",
        default=DEFAULT_JSON_NAME,
        help="Filename for LOC JSON payload (default: loc_metrics.json).",
    )
    parser.add_argument(
        "--image-name",
        default=DEFAULT_IMAGE_NAME,
        help="Filename for LOC PNG output (default: loc_comparison.png).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    loc_data = base_benchmark.compute_loc_benchmark()
    output_dir: Path = args.output_dir

    json_path = output_dir / args.json_name
    write_loc_json(loc_data, json_path)

    image_path = output_dir / args.image_name
    base_benchmark.generate_loc_slide(loc_data, image_path)

    print("LOC artefacts generated:")
    print(f"  • {json_path}")
    print(f"  • {image_path}")


if __name__ == "__main__":
    main()
