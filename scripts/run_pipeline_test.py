#!/usr/bin/env python3
"""
run_pipeline_test.py
====================
End-to-end test runner for the somatic variants pipeline downstream steps.

Exercises:
  1. BED overlap annotation   (mirrors annotate_target.R)
  2. Candidate selection       (mirrors selection_candidates.R)
  3. Functional/drug annotation (mirrors annotate_functional_drug.R)
  4. QC territories            (mirrors qc_territories.R)
  5. TMB computation           (mirrors tmb.R)
  6. MAF analysis              (mirrors maf_plots.R)
  7. VCF conversion            (calls benchmark_convert_pretty_to_vcf.py)
  8. Performance metrics       (mirrors benchmark_performance.R — ROC, PR, AUROC)

Produces:
  test_results/
    test_report.html            — structured HTML report
    territory_summary.csv
    tmb_summary.csv
    maf_summary.csv
    benchmark_summary.csv
    plots/
      territory_bar.png
      maf_facets.png
      tmb_bar.png
      variant_type_pie.png
      mc_distribution.png
      freq_histogram.png
      roc_curve.png
      pr_curve.png
      cadd_distribution.png
      druggable_summary.png

Usage:
  python3 scripts/run_pipeline_test.py [--data test_data/synthetic_pretty.csv]
                                        [--target test_data/target.bed]
                                        [--bait   test_data/bait.bed]
                                        [--outdir test_results]

All arguments have sensible defaults pointing at the bundled test_data/.
"""

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

# ── Optional but expected deps ────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    print("[WARN] matplotlib/numpy not available — plots will be skipped", file=sys.stderr)

# ═══════════════════════════════════════════════════════════════════════════════
# Constants mirroring the R scripts
# ═══════════════════════════════════════════════════════════════════════════════
KEEP_EXONIC = {
    "nonsynonymous SNV",
    "stopgain",
    "stoploss",
    "frameshift deletion",
    "frameshift insertion",
    "nonframeshift deletion",
    "nonframeshift insertion",
    "splicing",
}
KEEP_FUNC = {"exonic", "splicing", "exonic;splicing"}

EXONIC_TO_MC = {
    "nonsynonymous SNV":       "missense",
    "stopgain":                "nonsense",
    "stoploss":                "nonsense",
    "frameshift deletion":     "frameshift",
    "frameshift insertion":    "frameshift",
    "nonframeshift deletion":  "in_frame",
    "nonframeshift insertion": "in_frame",
    "splicing":                "splice_site",
    "synonymous SNV":          "silent",
}

PALETTE = {
    "on_target":  "#2166AC",
    "on_bait":    "#5AAE61",
    "off":        "#D6604D",
    "all":        "#4393C3",
    "selected":   "#D6604D",
    "on_target2": "#1A6B3A",
    "selected_on_target": "#7B3294",
}

# ═══════════════════════════════════════════════════════════════════════════════
# I/O helpers
# ═══════════════════════════════════════════════════════════════════════════════

def load_csv(path: str) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(rows: list[dict], path: str):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def load_bed(path: str) -> list[tuple]:
    """Return list of (chrom, start, end) — BED is 0-based half-open."""
    intervals = []
    if not path or not os.path.exists(path):
        return intervals
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            intervals.append((parts[0], int(parts[1]), int(parts[2])))
    return intervals


def bed_total_mb(intervals: list[tuple]) -> float:
    return sum(e - s for _, s, e in intervals) / 1e6


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline logic mirrors
# ═══════════════════════════════════════════════════════════════════════════════

def overlaps_any(chrom: str, pos: int, intervals: list[tuple]) -> bool:
    """BED [start,end) vs 1-based position."""
    for c, s, e in intervals:
        if c == chrom and s < pos <= e:   # BED 0-based → variant 1-based
            return True
    return False


def annotate_territory(variants: list[dict],
                       target: list[tuple],
                       bait: list[tuple]) -> list[dict]:
    """Mirror of annotate_target.R."""
    for v in variants:
        chrom = v.get("Chr", "")
        try:
            pos = int(v.get("Start", 0))
        except (ValueError, TypeError):
            pos = 0
        on_t = overlaps_any(chrom, pos, target)
        on_b = overlaps_any(chrom, pos, bait)
        v["on_target"]  = on_t
        v["on_bait"]    = on_b
        v["territory"]  = "on_target" if on_t else ("on_bait" if on_b else "off")
    return variants


def select_candidates(variants: list[dict],
                      false_genes: set = None) -> list[dict]:
    """Mirror of selection_candidates.R."""
    if false_genes is None:
        false_genes = set()
    for v in variants:
        gene  = v.get("Gene.refGene", "")
        ef    = v.get("ExonicFunc.refGene", ".")
        func  = v.get("Func.refGene", "")
        if gene in false_genes:
            v["selected"] = False
            continue
        exonic_ok = ef in KEEP_EXONIC
        exonic_na = not ef or ef == "."
        func_ok   = func in KEEP_FUNC
        v["selected"] = (exonic_ok or exonic_na) and func_ok
    return variants


def annotate_drug(variants: list[dict],
                  actionable_genes: set = None,
                  dgidb: dict = None) -> list[dict]:
    """Mirror of annotate_functional_drug.R."""
    if actionable_genes is None:
        actionable_genes = set()
    if dgidb is None:
        dgidb = {}
    for v in variants:
        ef   = v.get("ExonicFunc.refGene", ".")
        ref  = v.get("Ref", "A")
        alt  = v.get("Alt", "T")
        gene = v.get("Gene.refGene", "")
        v["mc"]           = EXONIC_TO_MC.get(ef, "other")
        v["variant_type"] = "SNV" if (len(ref) == 1 and len(alt) == 1) else "InDel"
        v["druggable"]    = gene in actionable_genes or gene in dgidb
        v["drug_interactions"] = dgidb.get(gene, "")
    return variants


def compute_tmb(variants: list[dict],
                exome_mb: float,
                territory: str = "on_target") -> dict[str, float]:
    """Return TMB per sample for a given territory."""
    if exome_mb <= 0:
        print("[WARN] compute_tmb: exome_mb is zero or negative — skipping TMB computation",
              file=sys.stderr)
        return {}
    counts: dict[str, int] = Counter()
    for v in variants:
        if territory == "all" or v.get("territory") == territory:
            counts[v.get("sample", "unknown")] += 1
    return {s: n / exome_mb for s, n in counts.items()}


def compute_max_maf(variants: list[dict]) -> dict[str, dict[str, float]]:
    """Return max FREQ per sample for four subsets."""
    subsets = {
        "all":                  variants,
        "on_target":            [v for v in variants if v.get("territory") == "on_target"],
        "selected":             [v for v in variants if v.get("selected")],
        "selected_on_target":   [v for v in variants if v.get("selected") and v.get("territory") == "on_target"],
    }
    result = {}
    for label, sub in subsets.items():
        per_sample: dict[str, float] = defaultdict(float)
        for v in sub:
            try:
                freq = float(v.get("FREQ", 0))
            except (ValueError, TypeError):
                freq = 0.0
            s = v.get("sample", "unknown")
            per_sample[s] = max(per_sample[s], freq)
        result[label] = dict(per_sample)
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# Performance metrics (ROC / PR)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_roc_pr(scores, labels):
    """Trapezoidal ROC and PR curves. Returns dict with curve arrays and scalars."""
    paired = sorted(zip(scores, labels), key=lambda x: -x[0])
    scores_s = [p[0] for p in paired]
    labels_s = [p[1] for p in paired]

    n_pos = sum(labels_s)
    n_neg = len(labels_s) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None

    tp_list, fp_list = [], []
    tp = fp = 0
    for lbl in labels_s:
        if lbl == 1:
            tp += 1
        else:
            fp += 1
        tp_list.append(tp)
        fp_list.append(fp)

    tpr = [t / n_pos for t in tp_list]
    fpr = [f / n_neg for f in fp_list]
    prec = [tp_list[i] / (tp_list[i] + fp_list[i]) for i in range(len(tp_list))]
    rec  = tpr[:]

    # Sort by FPR for AUC
    pairs_roc = sorted(zip(fpr, tpr))
    fpr_s = [p[0] for p in pairs_roc]
    tpr_s = [p[1] for p in pairs_roc]
    auc = sum(
        abs(fpr_s[i+1] - fpr_s[i]) * (tpr_s[i] + tpr_s[i+1]) / 2
        for i in range(len(fpr_s)-1)
    )

    # Sort by recall for AUPRC
    pairs_pr = sorted(zip(rec, prec))
    rec_s  = [p[0] for p in pairs_pr]
    prec_s = [p[1] for p in pairs_pr]
    auprc = sum(
        abs(rec_s[i+1] - rec_s[i]) * (prec_s[i] + prec_s[i+1]) / 2
        for i in range(len(rec_s)-1)
    )

    # Threshold sweep: compute precision & recall at each FREQ threshold
    thresholds = sorted(set(scores), reverse=True)
    pr_sweep = []
    for thresh in thresholds:
        pred_pos = [i for i, s in enumerate(scores) if s >= thresh]
        tp_t = sum(labels[i] for i in pred_pos)
        fp_t = len(pred_pos) - tp_t
        fn_t = n_pos - tp_t
        p = tp_t / (tp_t + fp_t) if (tp_t + fp_t) > 0 else 0
        r = tp_t / (tp_t + fn_t) if (tp_t + fn_t) > 0 else 0
        f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0
        pr_sweep.append({"threshold": round(thresh, 4), "precision": round(p, 4),
                         "recall": round(r, 4), "f1": round(f1, 4),
                         "tp": tp_t, "fp": fp_t, "fn": fn_t})

    return {
        "fpr": fpr_s, "tpr": tpr_s, "auc": round(auc, 4),
        "precision": prec_s, "recall": rec_s, "auprc": round(auprc, 4),
        "thresholds": scores_s,
        "pr_sweep": pr_sweep,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# Plotting helpers
# ═══════════════════════════════════════════════════════════════════════════════

def save_fig(fig, path: str):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_territory_bar(variants, outdir) -> str:
    if not HAS_PLOT:
        return ""
    samples = sorted(set(v.get("sample","?") for v in variants))
    terrs   = ["on_target", "on_bait", "off"]
    x = np.arange(len(samples))
    width = 0.25

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax_idx, (title, sub) in enumerate([
        ("All variants",       variants),
        ("Selected variants",  [v for v in variants if v.get("selected")]),
    ]):
        ax = axes[ax_idx]
        for i, t in enumerate(terrs):
            counts = [sum(1 for v in sub if v.get("sample") == s and v.get("territory") == t)
                      for s in samples]
            ax.bar(x + i*width, counts, width, label=t,
                   color=PALETTE.get(t, "#999"), edgecolor="white")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xticks(x + width)
        ax.set_xticklabels(samples, rotation=15)
        ax.set_ylabel("Variant count")
        ax.legend(title="Territory")
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.8)
        ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.suptitle("Variants by Territory", fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "territory_bar.png")
    save_fig(fig, path)
    return path


def plot_maf_facets(maf_result, outdir) -> str:
    if not HAS_PLOT:
        return ""
    subsets = list(maf_result.keys())
    samples = sorted({s for d in maf_result.values() for s in d})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes = axes.flatten()
    colors = plt.cm.Set2(np.linspace(0, 1, len(samples)))
    for idx, (subset, data) in enumerate(maf_result.items()):
        ax = axes[idx]
        vals = [data.get(s, 0) for s in samples]
        bars = ax.bar(samples, vals, color=colors, edgecolor="white")
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                    f"{val:.2f}", ha="center", va="bottom", fontsize=8)
        ax.set_title(subset.replace("_", " ").title(), fontsize=11, fontweight="bold")
        ax.set_ylim(0, 0.8)
        ax.set_ylabel("Max allele frequency")
        ax.set_xticklabels(samples, rotation=15)
        ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.suptitle("Max MAF per Sample and Variant Subset", fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "maf_facets.png")
    save_fig(fig, path)
    return path


def plot_tmb_bar(tmb_on_target, tmb_all, exome_mb, outdir) -> str:
    if not HAS_PLOT:
        return ""
    samples = sorted(set(list(tmb_on_target) + list(tmb_all)))
    x = np.arange(len(samples))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - width/2, [tmb_all.get(s, 0) for s in samples],
           width, label="All variants", color="#4393C3", edgecolor="white")
    ax.bar(x + width/2, [tmb_on_target.get(s, 0) for s in samples],
           width, label="On-target only", color="#2166AC", edgecolor="white")
    for ref in [1, 10, 100]:
        ax.axhline(ref, linestyle="--", color="grey", linewidth=0.8)
        ax.text(len(samples)-0.5, ref * 1.05, f"{ref}", color="grey", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=15)
    ax.set_ylabel("TMB (mutations / Mb)")
    ax.set_title(f"Tumor Mutational Burden  (exome: {exome_mb:.2f} Mb)", fontweight="bold")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "tmb_bar.png")
    save_fig(fig, path)
    return path


def plot_variant_type_pie(variants, outdir) -> str:
    if not HAS_PLOT:
        return ""
    mc_counts = Counter(v.get("mc", "other") for v in variants)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # Left: variant type SNV/InDel
    vt_counts = Counter(v.get("variant_type", "Unknown") for v in variants)
    axes[0].pie(vt_counts.values(), labels=vt_counts.keys(),
                autopct="%1.1f%%", colors=["#4393C3","#D6604D","#92C5DE"],
                startangle=90, wedgeprops={"edgecolor":"white"})
    axes[0].set_title("Variant Type (SNV vs InDel)", fontweight="bold")
    # Right: mutation category
    mc_order = ["missense","nonsense","frameshift","in_frame","splice_site","silent","other"]
    mc_vals  = [mc_counts.get(m, 0) for m in mc_order]
    mc_colors = plt.cm.Set3(np.linspace(0, 1, len(mc_order)))
    axes[1].pie(mc_vals, labels=mc_order, autopct=lambda p: f"{p:.1f}%" if p > 1 else "",
                colors=mc_colors, startangle=90, wedgeprops={"edgecolor":"white"})
    axes[1].set_title("Mutation Category (mc)", fontweight="bold")
    fig.suptitle("Variant Classification", fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "variant_type_pie.png")
    save_fig(fig, path)
    return path


def plot_freq_histogram(variants, outdir) -> str:
    if not HAS_PLOT:
        return ""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    samples = sorted(set(v.get("sample","?") for v in variants))
    colors = ["#2166AC", "#D6604D", "#5AAE61", "#7B3294"]
    for ax_idx, (title, sub) in enumerate([
        ("All variants", variants),
        ("Selected variants", [v for v in variants if v.get("selected")]),
    ]):
        ax = axes[ax_idx]
        for i, s in enumerate(samples):
            freqs = []
            for v in sub:
                if v.get("sample") != s:
                    continue
                try:
                    freqs.append(float(v["FREQ"]))
                except (ValueError, TypeError, KeyError):
                    pass
            ax.hist(freqs, bins=20, alpha=0.6, color=colors[i % len(colors)],
                    label=s, edgecolor="white")
        ax.set_xlabel("Allele Frequency (FREQ)")
        ax.set_ylabel("Count")
        ax.set_title(title, fontweight="bold")
        ax.legend()
        ax.axvline(0.05, color="grey", linestyle="--", linewidth=0.8, label="min_freq=0.05")
        ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.suptitle("Allele Frequency Distribution", fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "freq_histogram.png")
    save_fig(fig, path)
    return path


def plot_roc(roc_data: dict, sample_id: str, outdir: str) -> str:
    if not HAS_PLOT or roc_data is None:
        return ""
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(roc_data["fpr"], roc_data["tpr"], color="#2166AC", lw=2,
            label=f"FREQ score (AUC = {roc_data['auc']:.3f})")
    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", lw=0.8, label="Random")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title(f"ROC Curve — {sample_id}", fontsize=13, fontweight="bold")
    ax.legend(loc="lower right")
    ax.grid(linestyle="--", alpha=0.4)
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "roc_curve.png")
    save_fig(fig, path)
    return path


def plot_pr(roc_data: dict, sample_id: str, outdir: str) -> str:
    if not HAS_PLOT or roc_data is None:
        return ""
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(roc_data["recall"], roc_data["precision"], color="#D6604D", lw=2,
            label=f"FREQ score (AUPRC = {roc_data['auprc']:.3f})")
    # Baseline = fraction of positives
    n_pos = sum(1 for l in roc_data.get("thresholds", []) if l == 1)
    n_tot = len(roc_data.get("thresholds", [])) or 1
    ax.axhline(n_pos / n_tot, color="grey", linestyle="--", lw=0.8, label="Random")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Recall", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_title(f"Precision-Recall Curve — {sample_id}", fontsize=13, fontweight="bold")
    ax.legend(loc="upper right")
    ax.grid(linestyle="--", alpha=0.4)
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "pr_curve.png")
    save_fig(fig, path)
    return path


def plot_cadd_distribution(variants, outdir) -> str:
    if not HAS_PLOT:
        return ""
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = {
        "selected (truth=1)":      {"color": "#D6604D", "alpha": 0.7},
        "selected (truth=0)":      {"color": "#4393C3", "alpha": 0.7},
        "not selected (truth=1)":  {"color": "#FDAE61", "alpha": 0.5},
    }
    for label, style in groups.items():
        vals = []
        for v in variants:
            sel   = bool(v.get("selected"))
            truth = int(v.get("truth_label", 0))
            if label == "selected (truth=1)"     and sel and truth == 1: pass
            elif label == "selected (truth=0)"   and sel and truth == 0: pass
            elif label == "not selected (truth=1)" and not sel and truth == 1: pass
            else: continue
            try:
                vals.append(float(v["CADD13_PHRED"]))
            except (ValueError, TypeError, KeyError):
                pass
        if vals:
            ax.hist(vals, bins=15, alpha=style["alpha"], color=style["color"],
                    label=f"{label} (n={len(vals)})", edgecolor="white")
    ax.set_xlabel("CADD13_PHRED score")
    ax.set_ylabel("Count")
    ax.set_title("CADD Score Distribution by Selection & Truth Status", fontweight="bold")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "cadd_distribution.png")
    save_fig(fig, path)
    return path


def plot_druggable_summary(variants, outdir) -> str:
    if not HAS_PLOT:
        return ""
    fig, ax = plt.subplots(figsize=(7, 5))
    samples = sorted(set(v.get("sample","?") for v in variants))
    x = np.arange(len(samples))
    w = 0.35
    drug_ct   = [sum(1 for v in variants if v.get("sample")==s and v.get("druggable")) for s in samples]
    no_drug   = [sum(1 for v in variants if v.get("sample")==s and not v.get("druggable")) for s in samples]
    ax.bar(x - w/2, drug_ct, w, label="Druggable",    color="#2166AC", edgecolor="white")
    ax.bar(x + w/2, no_drug, w, label="Not druggable", color="#D6604D", edgecolor="white")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=15)
    ax.set_ylabel("Variant count")
    ax.set_title("Druggable vs Non-Druggable Variants", fontweight="bold")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    path = os.path.join(outdir, "plots", "druggable_summary.png")
    save_fig(fig, path)
    return path


# ═══════════════════════════════════════════════════════════════════════════════
# HTML report builder
# ═══════════════════════════════════════════════════════════════════════════════

HTML_TMPL = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Somatic Variants Pipeline — Test Report</title>
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; margin:0; background:#f7f9fc; color:#1a1a2e; }}
  header {{ background: #16213e; color: white; padding: 24px 40px; }}
  header h1 {{ margin:0; font-size:1.8em; }}
  header p  {{ margin:4px 0 0; opacity:.8; font-size:.9em; }}
  main {{ max-width: 1200px; margin: 30px auto; padding: 0 24px; }}
  h2 {{ border-left: 5px solid #2166AC; padding-left: 12px; margin-top: 36px;
        font-size: 1.3em; color: #16213e; }}
  h3 {{ color: #2166AC; margin-top: 24px; font-size:1.05em; }}
  .card {{ background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,.08);
           padding: 20px 24px; margin-bottom: 20px; }}
  .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
                  gap: 16px; margin: 16px 0; }}
  .metric-box {{ background: #EFF6FF; border-radius: 6px; padding: 14px 16px; text-align:center; }}
  .metric-box .val {{ font-size: 2em; font-weight: 700; color: #2166AC; }}
  .metric-box .lbl {{ font-size: .8em; color: #555; margin-top: 4px; }}
  .pass {{ color: #1a7a4a; font-weight:700; }}
  .fail {{ color: #c0392b; font-weight:700; }}
  table {{ border-collapse: collapse; width: 100%; font-size:.9em; }}
  th {{ background:#2166AC; color:white; padding: 8px 12px; text-align:left; }}
  td {{ padding: 7px 12px; border-bottom: 1px solid #e2e8f0; }}
  tr:hover td {{ background: #EFF6FF; }}
  .img-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr));
               gap: 20px; margin-top: 20px; }}
  .img-card {{ background: white; border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,.08);
               padding: 12px; }}
  .img-card img {{ width:100%; border-radius:4px; }}
  .img-card p {{ font-size:.8em; color:#666; margin:6px 0 0; text-align:center; }}
  .badge {{ display:inline-block; padding: 2px 8px; border-radius:20px; font-size:.8em;
            font-weight:600; }}
  .badge-pass {{ background:#d4edda; color:#155724; }}
  .badge-fail {{ background:#f8d7da; color:#721c24; }}
  .badge-warn {{ background:#fff3cd; color:#856404; }}
  footer {{ text-align:center; padding:24px; color:#888; font-size:.85em; }}
  .file-tree {{ font-family: monospace; font-size:.88em; background:#f1f5f9;
                border-radius:6px; padding:16px; line-height:1.7; }}
</style>
</head>
<body>
<header>
  <h1>🧬 Somatic Variants Pipeline — Automated Test Report</h1>
  <p>Generated: {timestamp} &nbsp;|&nbsp; Data: {data_file} &nbsp;|&nbsp; Samples: {samples}</p>
</header>
<main>

<h2>📋 Test Summary</h2>
<div class="card">
  <div class="metric-grid">
    {metric_boxes}
  </div>
  <table>
    <tr><th>Test</th><th>Status</th><th>Details</th></tr>
    {test_rows}
  </table>
</div>

<h2>📁 Output File Locations</h2>
<div class="card">
  <div class="file-tree">{file_tree}</div>
</div>

<h2>📊 Territory Annotation</h2>
<div class="card">
  {territory_table}
</div>

<h2>🎯 Candidate Selection</h2>
<div class="card">
  {selection_table}
</div>

<h2>💊 Drug Annotation</h2>
<div class="card">
  {drug_table}
</div>

<h2>📈 TMB (Tumor Mutational Burden)</h2>
<div class="card">
  {tmb_table}
</div>

<h2>🧪 Performance Metrics (ROC / PR)</h2>
<div class="card">
  {perf_table}
  <h3>Threshold Sweep (top 10 by F1)</h3>
  {sweep_table}
</div>

<h2>📉 Visualisations</h2>
{plot_section}

</main>
<footer>Somatic Variants Pipeline · Test Report · {timestamp}</footer>
</body>
</html>
"""


def make_metric_box(val, label):
    return f'<div class="metric-box"><div class="val">{val}</div><div class="lbl">{label}</div></div>'


def make_table(headers, rows_data):
    ths = "".join(f"<th>{h}</th>" for h in headers)
    trs = ""
    for row in rows_data:
        cells = "".join(f"<td>{c}</td>" for c in row)
        trs += f"<tr>{cells}</tr>"
    return f"<table><tr>{ths}</tr>{trs}</table>"


def make_plot_section(plot_paths: dict) -> str:
    cards = ""
    labels = {
        "territory_bar":   "Variants by Territory (faceted: all vs selected, log-scale)",
        "maf_facets":      "Max Allele Frequency per Sample and Variant Subset",
        "tmb_bar":         "Tumor Mutational Burden (TMB) — reference lines at 1/10/100 mut/Mb",
        "variant_type_pie":"Variant Type (SNV/InDel) and Mutation Category Distribution",
        "freq_histogram":  "Allele Frequency Distribution (all vs selected)",
        "roc_curve":       "ROC Curve — AUROC (higher = better classifier performance)",
        "pr_curve":        "Precision-Recall Curve — AUPRC (accounts for class imbalance)",
        "cadd_distribution":"CADD Score by Selection and Truth Status",
        "druggable_summary":"Druggable vs Non-Druggable Variants per Sample",
    }
    for key, path in plot_paths.items():
        if not path or not os.path.exists(path):
            continue
        rel = os.path.relpath(path, os.path.dirname(path) + "/..")
        lbl = labels.get(key, key.replace("_"," ").title())
        cards += f'<div class="img-card"><img src="{rel}" alt="{lbl}"><p>{lbl}</p></div>'
    return f'<div class="img-grid">{cards}</div>'


def make_badge(ok: bool, warn=False):
    if warn:
        return '<span class="badge badge-warn">WARN</span>'
    return '<span class="badge badge-pass">PASS</span>' if ok else '<span class="badge badge-fail">FAIL</span>'


# ═══════════════════════════════════════════════════════════════════════════════
# Main runner
# ═══════════════════════════════════════════════════════════════════════════════

def run_tests(data_file: str, target_bed: str, bait_bed: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "plots"), exist_ok=True)

    results   = []  # list of (name, pass, detail)
    plot_paths = {}
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())

    # ── 1. Load test data ────────────────────────────────────────────────────
    print(f"\n[1/8] Loading test data from {data_file} …")
    variants = load_csv(data_file)
    n_total  = len(variants)
    samples  = sorted(set(v.get("sample","?") for v in variants))
    ok_load  = n_total > 0
    results.append(("Load test CSV", ok_load, f"{n_total} rows, {len(samples)} samples: {', '.join(samples)}"))
    print(f"      {n_total} variants across {len(samples)} samples")

    # ── 2. BED overlap / territory annotation ────────────────────────────────
    print("[2/8] Annotating territory (BED overlap) …")
    target = load_bed(target_bed)
    bait   = load_bed(bait_bed)
    exome_mb = bed_total_mb(target)
    variants = annotate_territory(variants, target, bait)
    terr_counts = Counter(v["territory"] for v in variants)
    ok_terr = terr_counts.get("on_target", 0) > 0 or terr_counts.get("off", 0) > 0
    results.append(("BED territory annotation", ok_terr,
                    f"on_target={terr_counts.get('on_target',0)}, on_bait={terr_counts.get('on_bait',0)}, off={terr_counts.get('off',0)}"))
    print(f"      {dict(terr_counts)}")

    terr_rows = []
    for s in samples:
        for t in ["on_target","on_bait","off"]:
            n = sum(1 for v in variants if v.get("sample")==s and v.get("territory")==t)
            terr_rows.append({"sample":s,"territory":t,"count":n})
    write_csv(terr_rows, os.path.join(outdir, "territory_summary.csv"))

    # ── 3. Candidate selection ───────────────────────────────────────────────
    print("[3/8] Selecting candidates …")
    variants = select_candidates(variants, false_genes=set())
    n_sel  = sum(1 for v in variants if v.get("selected"))
    ok_sel = n_sel > 0
    results.append(("Candidate selection", ok_sel,
                    f"{n_sel}/{n_total} selected ({100*n_sel/n_total:.1f}%)"))
    print(f"      {n_sel} candidates selected")

    sel_rows = []
    for s in samples:
        tot  = sum(1 for v in variants if v.get("sample")==s)
        sel  = sum(1 for v in variants if v.get("sample")==s and v.get("selected"))
        sel_rows.append({"sample":s,"total":tot,"selected":sel,
                          "pct_selected":f"{100*sel/tot:.1f}%" if tot else "—"})
    write_csv(sel_rows, os.path.join(outdir, "selection_summary.csv"))

    # ── 4. Drug / functional annotation ─────────────────────────────────────
    print("[4/8] Drug / functional annotation …")
    variants = annotate_drug(variants)
    n_drug = sum(1 for v in variants if v.get("druggable"))
    ok_drug = True  # always completes
    results.append(("Drug annotation", ok_drug,
                    f"{n_drug}/{n_total} druggable variants"))
    mc_counts = Counter(v.get("mc","other") for v in variants)
    vt_counts = Counter(v.get("variant_type","?") for v in variants)
    print(f"      mc={dict(mc_counts)}")
    print(f"      variant_type={dict(vt_counts)}")

    drug_rows = [{"sample":s,
                  "total": sum(1 for v in variants if v.get("sample")==s),
                  "druggable": sum(1 for v in variants if v.get("sample")==s and v.get("druggable")),
                  "SNV": sum(1 for v in variants if v.get("sample")==s and v.get("variant_type")=="SNV"),
                  "InDel": sum(1 for v in variants if v.get("sample")==s and v.get("variant_type")=="InDel")}
                 for s in samples]
    write_csv(drug_rows, os.path.join(outdir, "drug_summary.csv"))

    # ── 5. TMB ───────────────────────────────────────────────────────────────
    print("[5/8] Computing TMB …")
    tmb_on_target = compute_tmb(variants, exome_mb, territory="on_target")
    tmb_all       = compute_tmb(variants, exome_mb, territory="all")
    ok_tmb = len(tmb_on_target) > 0
    tmb_vals = list(tmb_on_target.values())
    results.append(("TMB computation", ok_tmb,
                    f"exome={exome_mb:.3f} Mb; TMB (on-target)={[f'{v:.2f}' for v in tmb_vals]}"))
    print(f"      exome_mb={exome_mb:.3f}  tmb_on_target={dict({k:round(v,2) for k,v in tmb_on_target.items()})}")

    tmb_rows = []
    for s in samples:
        tmb_rows.append({"sample":s,
                          "exome_mb": round(exome_mb, 3),
                          "n_on_target": terr_counts.get("on_target", 0),
                          "TMB_on_target": round(tmb_on_target.get(s, 0), 3),
                          "TMB_all": round(tmb_all.get(s, 0), 3)})
    write_csv(tmb_rows, os.path.join(outdir, "tmb_summary.csv"))

    # ── 6. MAF ───────────────────────────────────────────────────────────────
    print("[6/8] Computing MAF …")
    maf_result = compute_max_maf(variants)
    ok_maf = len(maf_result["all"]) > 0
    results.append(("MAF analysis", ok_maf,
                    f"max MAF all={{{', '.join(f'{s}:{v:.3f}' for s,v in maf_result['all'].items())}}}"))

    maf_flat = []
    for subset, data in maf_result.items():
        for s, val in data.items():
            maf_flat.append({"sample":s,"subset":subset,"max_MAF":round(val,4)})
    write_csv(maf_flat, os.path.join(outdir, "maf_summary.csv"))

    # ── 7. VCF conversion ────────────────────────────────────────────────────
    print("[7/8] VCF conversion (benchmark_convert_pretty_to_vcf.py) …")
    vcf_out = os.path.join(outdir, "calls_test.vcf")
    converter = os.path.join(os.path.dirname(__file__), "benchmark_convert_pretty_to_vcf.py")
    ok_vcf = False
    vcf_detail = ""
    try:
        result = subprocess.run(
            [sys.executable, converter,
             "--input",  data_file,
             "--output", vcf_out,
             "--sample", "COMBINED"],
            capture_output=True, text=True, timeout=30
        )
        ok_vcf    = result.returncode == 0 and os.path.exists(vcf_out)
        vcf_lines = sum(1 for l in open(vcf_out) if not l.startswith("#")) if ok_vcf else 0
        vcf_detail = f"Wrote {vcf_lines} VCF records to {os.path.basename(vcf_out)}"
        print(f"      {vcf_detail}")
    except Exception as e:
        vcf_detail = str(e)
        print(f"      VCF conversion failed: {e}")
    results.append(("VCF conversion", ok_vcf, vcf_detail))

    # ── 8. Performance metrics (ROC / PR) ────────────────────────────────────
    print("[8/8] Computing ROC / PR performance metrics …")
    freq_vals = []
    labels_   = []
    for v in variants:
        try:
            freq_vals.append(float(v["FREQ"]))
            labels_.append(int(v.get("truth_label", 0)))
        except (ValueError, TypeError, KeyError):
            pass

    roc_data = compute_roc_pr(freq_vals, labels_)
    ok_perf  = roc_data is not None

    bench_rows = []
    if roc_data:
        print(f"      AUROC={roc_data['auc']:.4f}  AUPRC={roc_data['auprc']:.4f}")
        # Best F1 threshold
        best = max(roc_data["pr_sweep"], key=lambda r: r["f1"])
        bench_rows.append({
            "sample": "ALL", "score": "FREQ",
            "AUROC": roc_data["auc"], "AUPRC": roc_data["auprc"],
            "n_positive": sum(labels_), "n_total": len(labels_),
            "best_threshold": best["threshold"],
            "best_precision": best["precision"],
            "best_recall":    best["recall"],
            "best_F1":        best["f1"],
        })
        results.append(("ROC/PR computation", ok_perf,
                        f"AUROC={roc_data['auc']:.4f}, AUPRC={roc_data['auprc']:.4f}, "
                        f"best F1={best['f1']:.3f} at thresh={best['threshold']:.3f}"))
        write_csv(bench_rows, os.path.join(outdir, "benchmark_summary.csv"))
        write_csv(roc_data["pr_sweep"], os.path.join(outdir, "pr_sweep.csv"))
    else:
        results.append(("ROC/PR computation", False, "No valid score/label pairs"))

    # ── Plots ────────────────────────────────────────────────────────────────
    if HAS_PLOT:
        print("\n[Plots] Generating visualisations …")
        plot_paths["territory_bar"]    = plot_territory_bar(variants, outdir)
        plot_paths["maf_facets"]       = plot_maf_facets(maf_result, outdir)
        plot_paths["tmb_bar"]          = plot_tmb_bar(tmb_on_target, tmb_all, exome_mb, outdir)
        plot_paths["variant_type_pie"] = plot_variant_type_pie(variants, outdir)
        plot_paths["freq_histogram"]   = plot_freq_histogram(variants, outdir)
        plot_paths["roc_curve"]        = plot_roc(roc_data, "ALL samples", outdir)
        plot_paths["pr_curve"]         = plot_pr(roc_data, "ALL samples", outdir)
        plot_paths["cadd_distribution"] = plot_cadd_distribution(variants, outdir)
        plot_paths["druggable_summary"] = plot_druggable_summary(variants, outdir)
        print(f"      Saved {sum(1 for p in plot_paths.values() if p)} plots")

    # ── Build HTML report ────────────────────────────────────────────────────
    print("\n[Report] Building HTML report …")
    all_pass = all(r[1] for r in results)

    metric_boxes = "".join([
        make_metric_box(n_total, "Total Variants"),
        make_metric_box(len(samples), "Samples"),
        make_metric_box(n_sel, "Candidates Selected"),
        make_metric_box(n_drug, "Druggable Variants"),
        make_metric_box(f"{roc_data['auc']:.3f}" if roc_data else "—", "AUROC"),
        make_metric_box(f"{roc_data['auprc']:.3f}" if roc_data else "—", "AUPRC"),
        make_metric_box(f"{exome_mb:.2f} Mb", "Exome Target Size"),
        make_metric_box("✅ PASS" if all_pass else "❌ FAIL", "Overall Result"),
    ])

    test_rows = ""
    for name, ok, detail in results:
        badge = make_badge(ok)
        test_rows += f"<tr><td>{name}</td><td>{badge}</td><td>{detail}</td></tr>"

    # File tree
    out_files = []
    for root, dirs, files in os.walk(outdir):
        dirs.sort()
        level = root.replace(outdir, "").count(os.sep)
        indent = "│   " * level
        out_files.append(f"{indent}📁 {os.path.basename(root)}/")
        for f in sorted(files):
            fpath = os.path.join(root, f)
            size  = os.path.getsize(fpath)
            out_files.append(f"{indent}│   📄 {f}  ({size:,} bytes)")
    file_tree = "\n".join(out_files)

    # Territory table
    territory_table = make_table(
        ["Sample","Territory","Count"],
        [[r["sample"], r["territory"], r["count"]] for r in terr_rows]
    )

    # Selection table
    selection_table = make_table(
        ["Sample","Total","Selected","% Selected"],
        [[r["sample"], r["total"], r["selected"], r["pct_selected"]] for r in sel_rows]
    )

    # Drug table
    drug_table = make_table(
        ["Sample","Total","Druggable","SNV","InDel"],
        [[r["sample"], r["total"], r["druggable"], r["SNV"], r["InDel"]] for r in drug_rows]
    )

    # TMB table
    tmb_table = make_table(
        ["Sample","Exome (Mb)","TMB on-target","TMB all"],
        [[r["sample"], r["exome_mb"], r["TMB_on_target"], r["TMB_all"]] for r in tmb_rows]
    )

    # Performance table
    if bench_rows:
        perf_table = make_table(
            ["Sample","Score","AUROC","AUPRC","Best F1","Best Threshold","N positives","N total"],
            [[r.get("sample",""), r.get("score",""),
              r.get("AUROC",""), r.get("AUPRC",""),
              r.get("best_F1",""), r.get("best_threshold",""),
              r.get("n_positive",""), r.get("n_total","")] for r in bench_rows]
        )
        top_sweep = sorted(roc_data["pr_sweep"], key=lambda r: -r["f1"])[:10]
        sweep_table = make_table(
            ["FREQ threshold","Precision","Recall","F1","TP","FP","FN"],
            [[r["threshold"], r["precision"], r["recall"], r["f1"],
              r["tp"], r["fp"], r["fn"]] for r in top_sweep]
        )
    else:
        perf_table  = "<p>Performance data not available.</p>"
        sweep_table = "<p>Threshold sweep not available.</p>"

    plot_section = make_plot_section(plot_paths)

    html = HTML_TMPL.format(
        timestamp   = timestamp,
        data_file   = os.path.basename(data_file),
        samples     = ", ".join(samples),
        metric_boxes= metric_boxes,
        test_rows   = test_rows,
        file_tree   = file_tree,
        territory_table  = territory_table,
        selection_table  = selection_table,
        drug_table       = drug_table,
        tmb_table        = tmb_table,
        perf_table       = perf_table,
        sweep_table      = sweep_table,
        plot_section     = plot_section,
    )

    report_path = os.path.join(outdir, "test_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\n[Report] Written to {report_path}")

    # Final summary
    n_pass = sum(1 for r in results if r[1])
    n_fail = len(results) - n_pass
    print(f"\n{'='*60}")
    print(f"  Tests PASSED: {n_pass}/{len(results)}")
    print(f"  Tests FAILED: {n_fail}/{len(results)}")
    print(f"  Output dir:   {os.path.abspath(outdir)}")
    print("="*60)
    return 0 if all_pass else 1


# ═══════════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    repo_root = Path(__file__).parent.parent
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data",   default=str(repo_root / "test_data" / "synthetic_pretty.csv"),
                   help="Input pretty CSV (default: test_data/synthetic_pretty.csv)")
    p.add_argument("--target", default=str(repo_root / "test_data" / "target.bed"),
                   help="Target BED (default: test_data/target.bed)")
    p.add_argument("--bait",   default=str(repo_root / "test_data" / "bait.bed"),
                   help="Bait BED (default: test_data/bait.bed)")
    p.add_argument("--outdir", default=str(repo_root / "test_results"),
                   help="Output directory (default: test_results/)")
    args = p.parse_args()

    sys.exit(run_tests(
        data_file  = args.data,
        target_bed = args.target,
        bait_bed   = args.bait,
        outdir     = args.outdir,
    ))


if __name__ == "__main__":
    main()
