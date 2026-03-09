#!/usr/bin/env python3
"""
Merge MuTect + VarScan VCFs, filter by min_af/min_dp, build key/DP/FREQ/method, write ANNOVAR input.
HaTSPiL VariantCalling.prepare_for_annovar logic (simplified).
"""
import argparse
import re
import gzip

def parse_vcf(path, caller):
    rows = []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            if line.startswith("#"): continue
            toks = line.strip().split("\t")
            if len(toks) < 8: continue
            chrom, pos, ref, alt = toks[0], toks[1], toks[3], toks[4].split(",")[0]
            key = f"{chrom}:{pos}-{pos}_{ref}_{alt}"
            info = dict(x.split("=", 1) for x in toks[7].split(";") if "=" in x)
            dp = int(info.get("DP", 0))
            # FREQ from FORMAT if present
            if len(toks) > 9:
                fmt = toks[8].split(":")
                vals = toks[9].split(":")
                for i, k in enumerate(fmt):
                    if k == "FREQ" and i < len(vals):
                        freq = float(vals[i].rstrip("%")) / 100.0
                        break
                else:
                    freq = 0.0
            else:
                freq = 0.0
            rows.append((key, dp, freq, caller))
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample_id", required=True)
    p.add_argument("--mutect", required=True)
    p.add_argument("--varscan_snp", required=True)
    p.add_argument("--varscan_indel", required=True)
    p.add_argument("--min_af", type=float, default=0.01)
    p.add_argument("--min_dp", type=int, default=10)
    p.add_argument("--out_variants", required=True)
    p.add_argument("--out_annovar_input", required=True)
    args = p.parse_args()

    all_rows = []
    for path, caller in [(args.mutect, "Mutect1.17"), (args.varscan_snp, "VarScan2"), (args.varscan_indel, "VarScan2")]:
        try:
            all_rows.extend(parse_vcf(path, caller))
        except Exception as e:
            pass
    seen = {}
    for key, dp, freq, caller in all_rows:
        if dp < args.min_dp or freq < args.min_af: continue
        if key not in seen:
            seen[key] = [key, dp, freq, [caller]]
        else:
            if caller not in seen[key][3]:
                seen[key][3].append(caller)

    with open(args.out_variants, "w") as vf:
        vf.write("key,DP,FREQ,method\n")
        for key, (k, dp, freq, methods) in seen.items():
            vf.write(f"{k},{dp},{freq},{':'.join(methods)}\n")

    with open(args.out_annovar_input, "w") as af:
        for key in seen:
            m = re.match(r"^(chr)?(\d+|[XYxy]+):(\d+)-(\d+)_([ACGT]+)_([ACGT]+)$", key)
            if m:
                chrom, start, end, ref, alt = m.group(1) or "chr" + m.group(2), m.group(3), m.group(4), m.group(5), m.group(6)
                af.write(f"{chrom}\t{start}\t{end}\t{ref}\t{alt}\n")

if __name__ == "__main__":
    main()
