#!/usr/bin/env Rscript
# Annotate variants with on_target / on_bait / off (get_on_target logic).
# Performs real BED interval overlap using GenomicRanges when available,
# falling back to data.table interval arithmetic.

args <- commandArgs(trailingOnly = TRUE)
input_rds  <- gsub("^--input=",      "", args[grepl("^--input=",      args)])
target_bed <- gsub("^--target_bed=", "", args[grepl("^--target_bed=", args)])
bait_bed   <- gsub("^--bait_bed=",   "", args[grepl("^--bait_bed=",   args)])
output     <- gsub("^--output=",     "", args[grepl("^--output=",     args)])

library(data.table)

# Helper: load a BED file (tab-separated, no header); returns data.table with chr/start/end
read_bed <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  bed <- fread(path, header = FALSE, sep = "\t", select = 1:3,
               col.names = c("chr", "start", "end"), data.table = TRUE)
  bed[, chr   := as.character(chr)]
  bed[, start := as.integer(start)]
  bed[, end   := as.integer(end)]
  bed
}

# Helper: for each variant row, check overlap with a BED data.table using a fast interval join.
# Returns a logical vector (TRUE = overlaps at least one interval).
overlap_bed <- function(variants_dt, bed_dt) {
  if (is.null(bed_dt) || nrow(bed_dt) == 0L) return(rep(FALSE, nrow(variants_dt)))

  # Try GenomicRanges first (most accurate, handles strand/boundary correctly)
  if (requireNamespace("GenomicRanges", quietly = TRUE) &&
      requireNamespace("IRanges", quietly = TRUE)) {
    vgr  <- GenomicRanges::GRanges(
      seqnames = variants_dt$Chr,
      ranges   = IRanges::IRanges(start = as.integer(variants_dt$Start),
                                  end   = as.integer(variants_dt$End))
    )
  # BED format is 0-based half-open [start, end); convert to 1-based closed [start+1, end] for GRanges
    bgr  <- GenomicRanges::GRanges(
      seqnames = bed_dt$chr,
      ranges   = IRanges::IRanges(start = bed_dt$start + 1L, end = bed_dt$end)  # BED 0-based → 1-based
    )
    hits <- GenomicRanges::findOverlaps(vgr, bgr, ignore.strand = TRUE)
    result <- logical(nrow(variants_dt))
    result[S4Vectors::queryHits(hits)] <- TRUE
    return(result)
  }

  # Fallback: data.table non-equi join
  # BED is 0-based half-open [start, end) → convert to 1-based closed [start+1, end]
  vdt <- data.table(
    idx     = seq_len(nrow(variants_dt)),
    chr     = as.character(variants_dt$Chr),
    pos     = as.integer(variants_dt$Start),
    pos_end = as.integer(variants_dt$Start)
  )
  # BED [start, end) -> 1-based [start+1, end]
  bdt <- copy(bed_dt)
  bdt[, start1 := start + 1L]

  setkey(bdt, chr, start1, end)
  setkey(vdt, chr, pos, pos_end)
  hits <- foverlaps(vdt, bdt, by.x = c("chr", "pos", "pos_end"),
                    by.y = c("chr", "start1", "end"),
                    type = "within", nomatch = NULL)
  result <- logical(nrow(variants_dt))
  if (nrow(hits) > 0L) result[hits$idx] <- TRUE
  result
}

df         <- readRDS(input_rds)
target_dt  <- read_bed(target_bed)
bait_dt    <- read_bed(bait_bed)

if (all(c("Chr", "Start", "End") %in% names(df))) {
  vdt            <- as.data.table(df)
  vdt[, Chr   := as.character(Chr)]
  vdt[, Start := as.integer(Start)]
  vdt[, End   := as.integer(End)]

  df$on_target <- overlap_bed(vdt, target_dt)
  df$on_bait   <- overlap_bed(vdt, bait_dt)
  df$territory <- ifelse(df$on_target, "on_target",
                   ifelse(df$on_bait,  "on_bait", "off"))

  message(sprintf("Territory summary: on_target=%d  on_bait=%d  off=%d",
                  sum(df$on_target, na.rm = TRUE),
                  sum(df$on_bait,   na.rm = TRUE),
                  sum(df$territory == "off", na.rm = TRUE)))
} else {
  warning("Columns Chr/Start/End not found – setting territory='off' for all rows")
  df$on_target <- FALSE
  df$on_bait   <- FALSE
  df$territory <- "off"
}

saveRDS(df, output)
message("Wrote ", output)
