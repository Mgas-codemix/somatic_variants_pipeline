#!/usr/bin/env Rscript
# generate_example_plots.R
# Generates example QC plots (ggplot2) mimicking MultiQC / pipeline QC outputs.
# Usage: Rscript scripts/generate_example_plots.R [output_dir]
# Saves PNGs to output_dir (default: docs/test_report/plots).

suppressPackageStartupMessages({
  library(ggplot2)
})

cmd_args    <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(cmd_args) >= 1) cmd_args[1] else "docs/test_report/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, filename, width = 7, height = 4) {
  path <- file.path(output_dir, filename)
  ggsave(path, p, width = width, height = height, dpi = 150)
  message("Saved: ", path)
}

# ── 1. Per-base sequence quality ──────────────────────────────────────────────
set.seed(42)
pos  <- 1:150
qual <- 38 - (pos / 150) * 8 + rnorm(150, 0, 1)
qual <- pmax(pmin(qual, 40), 20)
df_qual <- data.frame(
  position = pos,
  median   = qual,
  q25      = qual - 2,
  q75      = qual + 2
)
p_qual <- ggplot(df_qual, aes(x = position)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#aed6f1", alpha = 0.6) +
  geom_line(aes(y = median), color = "#2e86c1", linewidth = 0.8) +
  geom_hline(yintercept = 30, linetype = "dashed", color = "#27ae60", linewidth = 0.5) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#e74c3c", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) +
  labs(
    title = "Per-Base Sequence Quality",
    x     = "Position in read (bp)",
    y     = "Phred quality score"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
save_plot(p_qual, "fastqc_per_base_quality.png")

# ── 2. Adapter content ────────────────────────────────────────────────────────
adapter_pct <- 0.5 * exp(seq(0, log(15), length.out = 150)) - 0.4
adapter_pct <- pmax(adapter_pct, 0)
df_adapter <- data.frame(
  position   = pos,
  adapter_r1 = adapter_pct,
  adapter_r2 = adapter_pct * 0.85
)
p_adapter <- ggplot(df_adapter, aes(x = position)) +
  geom_line(aes(y = adapter_r1, color = "R1"), linewidth = 0.8) +
  geom_line(aes(y = adapter_r2, color = "R2"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("R1" = "#2e86c1", "R2" = "#e67e22")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    title  = "Adapter Content",
    x      = "Position in read (bp)",
    y      = "% reads with adapter",
    color  = "Read"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
save_plot(p_adapter, "fastqc_adapter_content.png")

# ── 3. Alignment summary ──────────────────────────────────────────────────────
samples   <- paste0("Sample_", LETTERS[1:4])
aln_stats <- data.frame(
  sample    = rep(samples, 3),
  category  = rep(c("Mapped", "Duplicates", "Unmapped"), each = 4),
  pct       = c(97.2, 96.8, 98.1, 95.5,
                12.4, 14.1, 10.8, 16.3,
                2.8,  3.2,  1.9,  4.5)
)
p_aln <- ggplot(
  aln_stats[aln_stats$category != "Duplicates", ],
  aes(x = sample, y = pct, fill = category)
) +
  geom_bar(stat = "identity", position = "dodge", width = 0.65) +
  scale_fill_manual(values = c("Mapped" = "#27ae60", "Unmapped" = "#e74c3c")) +
  scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
  labs(
    title = "Alignment Summary",
    x     = NULL,
    y     = "% reads",
    fill  = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
save_plot(p_aln, "alignment_stats.png")

# ── 4. Insert size distribution ───────────────────────────────────────────────
ins_sizes <- data.frame(
  insert_size = rnorm(10000, mean = 185, sd = 30)
)
ins_sizes$insert_size <- round(pmax(ins_sizes$insert_size, 50))
p_ins <- ggplot(ins_sizes, aes(x = insert_size)) +
  geom_histogram(binwidth = 5, fill = "#8e44ad", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 185, linetype = "dashed", color = "#2c3e50", linewidth = 0.7) +
  labs(
    title = "Insert Size Distribution",
    x     = "Insert size (bp)",
    y     = "Count"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
save_plot(p_ins, "insert_size_distribution.png")

# ── 5. Duplication rate ───────────────────────────────────────────────────────
dup_data <- data.frame(
  sample  = samples,
  dup_pct = c(12.4, 14.1, 10.8, 16.3)
)
p_dup <- ggplot(dup_data, aes(x = sample, y = dup_pct, fill = sample)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(dup_pct, "%")), vjust = -0.4, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(limits = c(0, 25), labels = function(x) paste0(x, "%")) +
  labs(
    title = "Duplication Rate per Sample",
    x     = NULL,
    y     = "% duplicate reads"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
save_plot(p_dup, "duplication_rate.png")

# ── 6. Variant calling stats (Ti/Tv + counts) ─────────────────────────────────
vc_data <- data.frame(
  sample   = rep(samples, 2),
  category = rep(c("SNVs", "Indels"), each = 4),
  count    = c(1523, 1389, 1654, 1201,
               234,  198,  267,  187)
)
p_vc <- ggplot(vc_data, aes(x = sample, y = count, fill = category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65) +
  scale_fill_manual(values = c("SNVs" = "#2980b9", "Indels" = "#e67e22")) +
  labs(
    title = "Variant Counts by Type",
    x     = NULL,
    y     = "Number of variants",
    fill  = "Variant type"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
save_plot(p_vc, "variant_calling_stats.png")

message("\nAll example plots written to: ", output_dir)
