#!/usr/bin/env python3
"""
benchmark_convert_pretty_to_vcf.py
Convert a *_variants_pretty.csv to a minimal VCF for input to vcfeval or hap.py.
Optionally merges a truth-label column from a truth TSV/CSV (by Chr/Start/Ref/Alt key).

Usage:
    python benchmark_convert_pretty_to_vcf.py \
        --input  SAMPLE_variants_pretty.csv \
        --output SAMPLE_calls.vcf \
        [--truth  truth_variants.tsv]   # optional: adds TRUTH_LABEL INFO field
        [--sample SAMPLE_ID]
"""

import argparse
import csv
import gzip
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input",  required=True, help="Input pretty CSV (or .csv.gz)")
    p.add_argument("--output", required=True, help="Output VCF file")
    p.add_argument("--truth",  default="",    help="Optional truth TSV/CSV with same key columns")
    p.add_argument("--sample", default="SAMPLE", help="Sample name for VCF header")
    return p.parse_args()


def open_file(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def load_truth(truth_path):
    """Return set of (chr, pos, ref, alt) keys from a truth file."""
    if not truth_path or not os.path.exists(truth_path):
        return None
    truth_keys = set()
    with open_file(truth_path) as fh:
        sep = "\t" if truth_path.endswith((".tsv", ".tsv.gz")) else ","
        reader = csv.DictReader(fh, delimiter=sep)
        for row in reader:
            chr_ = row.get("Chr") or row.get("CHROM") or row.get("#CHROM", "")
            pos  = row.get("Start") or row.get("POS", "")
            ref  = row.get("Ref")  or row.get("REF", "")
            alt  = row.get("Alt")  or row.get("ALT", "")
            if chr_ and pos and ref and alt:
                truth_keys.add((chr_.lstrip("chr"), pos, ref, alt))
    return truth_keys


VCF_HEADER = """\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=FREQ,Number=1,Type=Float,Description="Allele frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=CADD,Number=1,Type=Float,Description="CADD13 PHRED score">
##INFO=<ID=METHOD,Number=1,Type=String,Description="Calling method(s)">
##INFO=<ID=TRUTH_LABEL,Number=1,Type=Integer,Description="1=in truth set 0=not">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def main():
    args = parse_args()
    truth_keys = load_truth(args.truth)
    has_truth  = truth_keys is not None

    required_missing = []
    rows_written = 0

    with open_file(args.input) as fh, open(args.output, "w") as out:
        out.write(VCF_HEADER.format(sample=args.sample))

        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []

        for row in reader:
            chrom = (row.get("Chr") or row.get("CHROM") or "").strip()
            pos   = (row.get("Start") or row.get("POS") or "").strip()
            ref   = (row.get("Ref")   or row.get("REF") or "").strip()
            alt   = (row.get("Alt")   or row.get("ALT") or "").strip()

            if not chrom or not pos or not ref or not alt:
                continue

            # Normalise chromosome name (VCF style: keep as-is; both chr-prefixed and bare are valid)
            vcf_chrom = chrom

            try:
                int_pos = int(pos)
            except ValueError:
                continue

            # Build INFO
            info_parts = []
            freq = (row.get("FREQ") or "").strip()
            if freq:
                try:
                    info_parts.append(f"FREQ={float(freq):.4f}")
                except ValueError:
                    pass
            dp = (row.get("DP") or "").strip()
            if dp:
                try:
                    info_parts.append(f"DP={int(float(dp))}")
                except ValueError:
                    pass
            cadd = (row.get("CADD13_PHRED") or "").strip()
            if cadd:
                try:
                    info_parts.append(f"CADD={float(cadd):.2f}")
                except ValueError:
                    pass
            method = (row.get("method") or "").strip().replace(" ", "_")
            if method:
                info_parts.append(f"METHOD={method}")

            if has_truth:
                key    = (chrom.lstrip("chr"), pos, ref, alt)
                label  = 1 if key in truth_keys else 0
                info_parts.append(f"TRUTH_LABEL={label}")

            info_str = ";".join(info_parts) if info_parts else "."

            out.write(
                f"{vcf_chrom}\t{int_pos}\t.\t{ref}\t{alt}\t.\tPASS\t"
                f"{info_str}\tGT\t0/1\n"
            )
            rows_written += 1

    print(f"Wrote {rows_written} variants to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
