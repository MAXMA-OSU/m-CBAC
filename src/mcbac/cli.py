"""Command-line interface for mcbac."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mcbac",
        description="Run m-CBAC charge assignment on CIF files.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the m-CBAC workflow")
    run_parser.add_argument("--data-dir", type=Path, required=True, help="Directory containing input .cif files")
    run_parser.add_argument("--output-dir", type=Path, required=True, help="Directory for FINAL_*.cif outputs")
    run_parser.add_argument("--log-file", type=Path, default=None, help="Path to write pipeline log")
    run_parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="Optional working directory (default: temporary directory)",
    )
    run_parser.add_argument(
        "--keep-work-dir",
        action="store_true",
        help="Keep working directory after completion",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run":
        result = run_pipeline(
            data_dir=args.data_dir,
            output_dir=args.output_dir,
            log_path=args.log_file,
            work_dir=args.work_dir,
            keep_work_dir=args.keep_work_dir,
        )
        print(f"Processed {result.input_count} CIF file(s).")
        print(f"Generated {len(result.final_files)} FINAL file(s).")
        for final_file in result.final_files:
            print(final_file)
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    sys.exit(main())