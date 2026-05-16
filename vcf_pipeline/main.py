#!/usr/bin/env python3
"""
main.py  —  VCF Variant Annotation Pipeline
============================================
CLI entry point. Parses a VCF file, filters variants by user-defined
conditions, and outputs an interactive HTML clinical report.

Usage
-----
    python main.py --vcf data/sample_vcfs/sample_diabetes.vcf \
                   --conditions "diabetes" "cad" \
                   --output reports/my_report.html \
                   --patient-id PATIENT_001

    # List available built-in condition aliases
    python main.py --list-conditions
"""

import argparse
import sys
from pathlib import Path

# ── make src importable whether run from project root or anywhere ──────────
sys.path.insert(0, str(Path(__file__).parent))

from src.vcf_parser      import parse_vcf
from src.condition_filter import filter_by_condition, group_by_condition, list_available_conditions
from src.report_generator import generate_report


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="vcf_pipeline",
        description="VCF Variant Annotation Pipeline — filter by condition and export HTML report",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument(
        "--vcf", "-v",
        metavar="FILE",
        help="Path to input VCF file (v4.x)",
    )
    p.add_argument(
        "--conditions", "-c",
        nargs="+",
        metavar="CONDITION",
        help='One or more condition keywords, e.g. "diabetes" "cad" "breast cancer"',
    )
    p.add_argument(
        "--output", "-o",
        metavar="FILE",
        default="reports/report.html",
        help="Output HTML report path (default: reports/report.html)",
    )
    p.add_argument(
        "--patient-id",
        metavar="ID",
        default="PATIENT_001",
        help="Patient identifier for the report header",
    )
    p.add_argument(
        "--include-benign",
        action="store_true",
        help="Include Benign / VUS variants in the report",
    )
    p.add_argument(
        "--include-non-carriers",
        action="store_true",
        help="Include homozygous-ref (non-carrier) records",
    )
    p.add_argument(
        "--list-conditions",
        action="store_true",
        help="Print all built-in condition aliases and exit",
    )
    return p


def main() -> int:
    parser = build_parser()
    args   = parser.parse_args()

    # ── List conditions mode ───────────────────────────────────────────────
    if args.list_conditions:
        print("\nAvailable built-in condition aliases:\n")
        for c in list_available_conditions():
            print(f"  • {c}")
        print("\nYou can also type any free-text condition — the pipeline will search for partial matches.\n")
        return 0

    # ── Validate required args ─────────────────────────────────────────────
    if not args.vcf:
        parser.error("--vcf is required (use --list-conditions to see available conditions)")
    if not args.conditions:
        parser.error("--conditions is required (e.g. --conditions diabetes cad)")

    vcf_path = Path(args.vcf)

    # ── Step 1: Parse VCF ──────────────────────────────────────────────────
    print(f"\n[1/4] Parsing VCF: {vcf_path} ...")
    try:
        variants = parse_vcf(vcf_path)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}")
        return 1

    print(f"      → {len(variants)} total records loaded")

    # ── Step 2: Filter by condition ────────────────────────────────────────
    print(f"[2/4] Filtering for conditions: {args.conditions} ...")
    filtered = filter_by_condition(
        variants,
        conditions=args.conditions,
        include_benign=args.include_benign,
        only_carriers=not args.include_non_carriers,
    )
    print(f"      → {len(filtered)} variant(s) matched")

    if not filtered:
        print("      ⚠  No variants matched. Try --include-benign or --include-non-carriers, or check condition spelling.")
        print("         Use --list-conditions to see available aliases.\n")
        return 0

    # ── Step 3: Group by condition ─────────────────────────────────────────
    print("[3/4] Grouping variants by condition ...")
    groups = group_by_condition(filtered, args.conditions)
    for cond, vars_ in groups.items():
        print(f"      {cond}: {len(vars_)} variant(s)")

    # ── Step 4: Generate HTML report ───────────────────────────────────────
    out_path = Path(args.output)
    print(f"[4/4] Generating HTML report → {out_path} ...")
    generate_report(
        groups=groups,
        vcf_filename=vcf_path.name,
        output_path=out_path,
        patient_id=args.patient_id,
    )

    print(f"\n✅ Report saved to: {out_path.resolve()}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
