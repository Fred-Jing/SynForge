from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


PRIMER_BUILD_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PRIMER_BUILD_ROOT / "src"))

from primer_build.project_pipeline import run_project  # noqa: E402


def positive_integer(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be an integer of at least 1")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Design mini-protein oligo/OE/Gibson constructs from a project "
            "folder containing structures/."
        )
    )
    parser.add_argument(
        "--project",
        required=True,
        help=(
            "Project directory under this checkout, or an absolute path. "
            "It must contain structures/."
        ),
    )
    parser.add_argument(
        "--chain",
        required=True,
        help="Target protein chain identifier, for example A",
    )
    parser.add_argument(
        "--map-start",
        required=True,
        help="First construct label, for example KYN-1",
    )
    parser.add_argument(
        "--primer-start",
        required=True,
        help="First outer-primer naming example, for example KYN-1-1-1-FW",
    )
    parser.add_argument(
        "--candidates",
        type=positive_integer,
        default=2,
        help="Synonymous CDS candidates evaluated per structure (default: 2)",
    )
    parser.add_argument(
        "--workers",
        type=positive_integer,
        default=min(4, os.cpu_count() or 1),
        help=(
            "Independent structure-design processes (default: up to 4). "
            "Start with 4 workers on Apple Silicon; benchmark before increasing."
        ),
    )
    parser.add_argument(
        "--rules",
        type=Path,
        help="Alternative design rules JSON (default: config/design_rules.json)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace only generated maps/workbooks with matching target names",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    project = Path(args.project)
    if not project.is_absolute():
        project = PRIMER_BUILD_ROOT / project
    try:
        run_project(
            primer_build_root=PRIMER_BUILD_ROOT,
            project_dir=project,
            chain=args.chain,
            map_start=args.map_start,
            primer_start=args.primer_start,
            max_candidates=args.candidates,
            workers=args.workers,
            overwrite=args.overwrite,
            rules_path=(args.rules if args.rules.is_absolute() else PRIMER_BUILD_ROOT / args.rules)
            if args.rules else None,
        )
    except KeyboardInterrupt:
        print(
            "\nSTOPPED by user. Completed maps, if any, remain in the "
            "project maps directory.",
            file=sys.stderr,
        )
        raise SystemExit(130)
