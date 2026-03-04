#!/usr/bin/env Rscript
# =============================================================================
#  Step 4a: Scale Factor Computation (called by 04_scale_deeptools.sh)
#
#  処理内容:
#    - グループごとの MergedBAM のピーク内リード数を csaw でカウント
#    - 先頭グループ (samples.tsv の出現順) を基準 (scale_factor=1.0) に正規化
#    - scale_factors.tsv を出力
#
#  環境変数:
#    PIPELINE_DIR, PIPELINE_SAMPLES, PIPELINE_BLACKLIST,
#    PIPELINE_THREADS, PIPELINE_PEAK_DIR, PIPELINE_MERGED_BAM_DIR
#
#  出力:
#    Peak_nomodel/scale_factors.tsv  (group <TAB> scale_factor)
# =============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(csaw)
  library(edgeR)
  library(tidyverse)
  library(data.table)
  library(BiocParallel)
})

# --- 環境変数読み込み ---
dir_proj       <- Sys.getenv("PIPELINE_DIR",          unset = getwd())
samples_tsv    <- Sys.getenv("PIPELINE_SAMPLES",      unset = file.path(dir_proj, "samples.tsv"))
blacklist_path <- Sys.getenv("PIPELINE_BLACKLIST",    unset = "")
n_threads      <- as.integer(Sys.getenv("PIPELINE_THREADS", unset = "4"))
peak_dir       <- Sys.getenv("PIPELINE_PEAK_DIR",    unset = file.path(dir_proj, "Peak_nomodel"))
merged_bam_dir <- Sys.getenv("PIPELINE_MERGED_BAM_DIR", unset = file.path(dir_proj, "MergedBAM"))

BPPARAM <- MulticoreParam(n_threads)

# --- ピークファイル読み込み ---
peak_bed <- file.path(peak_dir, "Peaks_250bp_th0.bed")
if (!file.exists(peak_bed)) {
  # フィルタ前のファイルを試みる
  peak_bed <- file.path(peak_dir, "Peak_250bp_noBL.bed")
}
stopifnot(file.exists(peak_bed))

peaks_dt <- fread(peak_bed, sep = "\t", header = FALSE, select = 1:3,
                  col.names = c("chrom", "start", "end"))
all.peaks <- GRanges(peaks_dt)
message("Peaks for scale factor computation: ", length(all.peaks))

# --- ブラックリスト読み込み ---
blacklist <- GRanges()
if (nzchar(blacklist_path) && file.exists(blacklist_path)) {
  bl_dt <- fread(blacklist_path, select = 1:3, col.names = c("chrom", "start", "end"))
  bl_dt[start == 0L, start := 1L]
  blacklist <- GRanges(bl_dt)
}

# --- MergedBAM グループ一覧 ---
samples_df <- read.table(samples_tsv, header = TRUE, sep = "\t",
                         comment.char = "#", stringsAsFactors = FALSE)
colnames(samples_df) <- c("sample_name", "group")
groups <- unique(samples_df$group)
message("Groups: ", paste(groups, collapse = ", "))

bam_files <- file.path(merged_bam_dir,
  paste0(groups, ".noBL.filtered.sorted.bam"))

missing <- bam_files[!file.exists(bam_files)]
if (length(missing) > 0) {
  stop("Merged BAM files not found:\n", paste(missing, collapse = "\n"))
}

# --- csaw カウント ---
message("Counting reads in peaks per group...")
is_human <- any(grepl("hg|human|GRCh", Sys.getenv("PIPELINE_GENOME", "hg38")))
chrom_numbers <- if (is_human) c(1:22, "X") else c(1:19, "X")
standard_chr <- paste0("chr", chrom_numbers)

param <- readParam(
  max.frag = 1000,
  pe = "both",
  discard = blacklist,
  restrict = standard_chr
)

peak.counts <- regionCounts(bam_files, all.peaks, param = param, BPPARAM = BPPARAM)
count_mat <- peak.counts@assays@data$counts

# --- 全リード数 (ピーク内) ---
peak_reads <- colSums(count_mat)

# --- 全リード数 (BAM 全体: samtools view -c) ---
message("Counting total mapped reads per group BAM...")
total_reads <- sapply(bam_files, function(bam) {
  as.numeric(system(paste("samtools view -c", shQuote(bam)), intern = TRUE))
})

# --- スケールファクター算出 ---
# CPM (= peak_reads / total_reads * 1e6)
# スケールファクター = CPM[ref] / CPM[i]  →  先頭グループが基準 (= 1.0)
cpm_per_group <- peak_reads / total_reads * 1e6
scale_factors  <- cpm_per_group[1] / cpm_per_group

message("Scale factors:")
for (i in seq_along(groups)) {
  message(sprintf("  %-20s  total_reads=%d  peak_reads=%d  scale=%.6f",
                  groups[i], total_reads[i], peak_reads[i], scale_factors[i]))
}

# --- TSV 出力 ---
out_tsv <- file.path(peak_dir, "scale_factors.tsv")
sf_df <- data.frame(group = groups, scale_factor = round(scale_factors, 6))
write.table(sf_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("Scale factors saved: ", out_tsv)
