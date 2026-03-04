#!/usr/bin/env Rscript
# =============================================================================
#  Step 3: Peak Counts & Threshold Filtering
#
#  処理内容:
#    1. Peak_250bp_noBL.bed の各ピークについて、全サンプル BAM の
#       リード数を csaw::regionCounts でカウント
#    2. aveLogCPM ヒストグラムを出力 → 低シグナルピーク除去の閾値を確認
#    3. PEAK_LOGCPM_THRESHOLD 未満のピークを除去
#    4. フィルタ後ピーク BED と生カウント行列を出力
#
#  環境変数 (run_pipeline.sh から自動 export):
#    PIPELINE_DIR, PIPELINE_SAMPLES, PIPELINE_BLACKLIST,
#    PIPELINE_THREADS, PIPELINE_THRESHOLD,
#    PIPELINE_PEAK_DIR, PIPELINE_BAM_DIR,
#    PIPELINE_FROM_SUBSTEP  ("" = full run, "b" = skip csaw; load RDS)
#
#  出力:
#    Peak_nomodel/Histogram_peaks.png
#    Peak_nomodel/Peaks_250bp_th0.bed   (フィルタ後 peak BED)
#    Peak_nomodel/PeakCounts_th0.txt    (フィルタ後 raw counts)
# =============================================================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(csaw)
  library(edgeR)
  library(ggplot2)
  library(tidyverse)
  library(data.table)
  library(BiocParallel)
})

# --- 環境変数から設定読み込み ---
dir_proj      <- Sys.getenv("PIPELINE_DIR",      unset = getwd())
samples_tsv   <- Sys.getenv("PIPELINE_SAMPLES",  unset = file.path(dir_proj, "samples.tsv"))
blacklist_path <- Sys.getenv("PIPELINE_BLACKLIST", unset = "")
n_threads     <- as.integer(Sys.getenv("PIPELINE_THREADS", unset = "4"))
peak_dir      <- Sys.getenv("PIPELINE_PEAK_DIR", unset = file.path(dir_proj, "Peak_nomodel"))
bam_dir       <- Sys.getenv("PIPELINE_BAM_DIR",  unset = file.path(dir_proj, "BAM"))
threshold_raw <- Sys.getenv("PIPELINE_THRESHOLD", unset = "")
from_sub      <- Sys.getenv("PIPELINE_FROM_SUBSTEP", unset = "")

# 閾値モードの判定
# "" → ヒストグラム出力のみ・停止  "auto" → 自動検出  数値 → そのまま使用
THRESHOLD_MODE <- if (nchar(trimws(threshold_raw)) == 0) {
  "interactive"
} else if (tolower(trimws(threshold_raw)) == "auto") {
  "auto"
} else {
  "fixed"
}
threshold <- if (THRESHOLD_MODE == "fixed") {
  val <- suppressWarnings(as.numeric(threshold_raw))
  if (is.na(val)) stop("PEAK_LOGCPM_THRESHOLD='", threshold_raw, "' は数値でありません。")
  val
} else {
  NA_real_
}

BPPARAM <- MulticoreParam(n_threads)

# --- ピークファイル読み込み ---
peak_bed <- file.path(peak_dir, "Peak_250bp_noBL.bed")
stopifnot(file.exists(peak_bed))
peaks_dt <- fread(peak_bed, sep = "\t", header = FALSE, select = 1:3,
                  col.names = c("chrom", "start", "end"))
all.peaks <- GRanges(peaks_dt)
message("Peaks loaded: ", length(all.peaks))

# --- ブラックリスト読み込み ---
blacklist <- GRanges()
if (nzchar(blacklist_path) && file.exists(blacklist_path)) {
  bl_dt <- fread(blacklist_path, select = 1:3, col.names = c("chrom", "start", "end"))
  bl_dt <- bl_dt[start > 0]       # start=0 対策
  bl_dt[start == 0L, start := 1L]
  blacklist <- GRanges(bl_dt)
}

# --- サンプルシートからBAMファイル取得 ---
samples_df <- read.table(samples_tsv, header = TRUE, sep = "\t",
                         comment.char = "#", stringsAsFactors = FALSE)
colnames(samples_df) <- c("sample_name", "group")

bam_files <- file.path(bam_dir,
  paste0(samples_df$sample_name, ".noDup.noMT.filt.sorted.bam"))

# BAMファイル存在確認
missing <- bam_files[!file.exists(bam_files)]
if (length(missing) > 0) {
  stop("BAM files not found:\n", paste(missing, collapse = "\n"))
}
message("BAM files: ", length(bam_files))

# --- csaw regionCounts ---
message("Counting reads in peaks (csaw)...")
# ゲノムの標準染色体を制限対象に設定
is_human <- any(grepl("hg|human|GRCh", Sys.getenv("PIPELINE_GENOME", "hg38")))
chrom_numbers <- if (is_human) c(1:22, "X") else c(1:19, "X")
standard_chr <- paste0("chr", chrom_numbers)

param <- readParam(
  max.frag = 1000,
  pe = "both",          # paired-end: use fragment (insert) for counting
  discard = blacklist,
  restrict = standard_chr
)

# --- csaw regionCounts (substep a) ---
rds_file <- file.path(peak_dir, "peak_counts.rds")

if (from_sub == "b") {
  if (!file.exists(rds_file)) {
    stop("RDS not found: ", rds_file,
         "\n  Run step 3 from the beginning first:  ./run_pipeline.sh --from 3")
  }
  message("[a] Skipped — loading cached peak counts from: ", rds_file)
  peak.counts <- readRDS(rds_file)
} else {
  message("[a] Counting reads in peaks (csaw)...")
  peak.counts <- regionCounts(bam_files, all.peaks, param = param, BPPARAM = BPPARAM)
  saveRDS(peak.counts, rds_file)
  message("    Cached to: ", rds_file)
}
peak.abundances <- aveLogCPM(asDGEList(peak.counts))

# --- ヒストグラム出力 (常に実行) ---
message("Plotting aveLogCPM histogram...")
hist_png <- file.path(peak_dir, "Histogram_peaks.png")

# aveLogCPM 分布を TSV でも保存 (閾値検討用)
logcpm_tsv <- file.path(peak_dir, "peak_logcpm.tsv")
fwrite(data.frame(peak_id = seq_along(peak.abundances),
                  aveLogCPM = round(peak.abundances, 4)),
       logcpm_tsv, sep = "\t")

df_hist <- data.frame(aveLogCPM = peak.abundances)

# auto モード: histogram の最初の valley を検出
if (THRESHOLD_MODE == "auto") {
  h <- hist(peak.abundances, breaks = 60, plot = FALSE)
  dens <- h$counts
  valley_idx <- which(diff(sign(diff(dens))) == 2)[1]
  if (!is.na(valley_idx) && length(valley_idx) > 0) {
    threshold <- h$breaks[valley_idx + 1]
    message("Auto threshold (histogram valley): ", round(threshold, 3))
  } else {
    threshold <- -1
    message("Auto threshold detection failed. Using default: ", threshold)
  }
}

# プロット描画 (interactive のときは threshold ラインなし)
vline_val <- if (THRESHOLD_MODE == "interactive") NA_real_ else threshold
p_hist <- ggplot(df_hist, aes(x = aveLogCPM)) +
  geom_histogram(binwidth = 0.2, fill = "#00DDAA", alpha = 0.6, color = "white") +
  { if (!is.na(vline_val))
      list(
        geom_vline(xintercept = vline_val, color = "red", linetype = "dashed", linewidth = 1),
        annotate("text", x = vline_val + 0.1, y = Inf, hjust = 0, vjust = 1.5,
                 label = paste0("threshold = ", round(vline_val, 2)),
                 color = "red", size = 4)
      )
  } +
  labs(title = "Peak Signal Distribution (aveLogCPM)",
       x = "Average log2 CPM", y = "Number of peaks") +
  theme_bw()
ggsave(hist_png, p_hist, width = 7, height = 5, dpi = 300)
message("Histogram saved: ", hist_png)

# =============================================================================
#  INTERACTIVE モード: 閾値を確認してから再実行するよう促して終了
# =============================================================================
if (THRESHOLD_MODE == "interactive") {
  msg <- paste0(
    "\n",
    "╔══════════════════════════════════════════════════════════════════╗\n",
    "║  ⏸  STEP 3 PAUSED — threshold not set                          ║\n",
    "╠══════════════════════════════════════════════════════════════════╣\n",
    "║  1. ヒストグラムを確認:                                          ║\n",
    "║       ", peak_dir, "/Histogram_peaks.png\n",
    "║  2. 低シグナルピークを除去する aveLogCPM の閾値を決定            ║\n",
    "║     (ヒストグラムの谷の値: 通常 -1〜1 程度)                      ║\n",
    "║  3. config.sh の PEAK_LOGCPM_THRESHOLD に値を設定:               ║\n",
    "║       PEAK_LOGCPM_THRESHOLD=\"-1\"   ← 例                         ║\n",
    "║  4. Step 3b に再実行 (カウントはキャッシュ使用):               ║\n",
    "║       ./run_pipeline.sh --from 3b                                ║\n",
    "╚══════════════════════════════════════════════════════════════════╝\n"
  )
  message(msg)
  # 終了コード 99 で停止 (run_pipeline.sh が検知して一時停止)
  quit(save = "no", status = 99)
}

# --- ピークフィルタリング ---
keep <- peak.abundances >= threshold
message("Peaks before filter: ", length(keep))
message("Peaks after  filter: ", sum(keep))
message("Peaks removed      : ", sum(!keep))

peak.counts.filt <- peak.counts[keep, ]
peaks.filt <- all.peaks[keep]

# --- カウント行列出力 ---
count_mat <- peak.counts.filt@assays@data$counts
sample_names <- samples_df$sample_name

out_counts <- data.frame(
  chr    = as.character(seqnames(peaks.filt)),
  start  = start(peaks.filt),
  end    = end(peaks.filt),
  width  = width(peaks.filt),
  strand = as.character(strand(peaks.filt)),
  as.data.frame(count_mat)
)
colnames(out_counts)[6:ncol(out_counts)] <- sample_names

counts_out <- file.path(peak_dir, "PeakCounts_th0.txt")
fwrite(out_counts, counts_out, sep = "\t")
message("Count matrix saved: ", counts_out)

# --- フィルタ後 peak BED 出力 ---
peaks_out <- file.path(peak_dir, "Peaks_250bp_th0.bed")
peaks_bed_df <- data.frame(
  chr   = as.character(seqnames(peaks.filt)),
  start = start(peaks.filt),
  end   = end(peaks.filt)
)
fwrite(peaks_bed_df, peaks_out, sep = "\t", col.names = FALSE)
message("Filtered peak BED saved: ", peaks_out)

message("")
message("Step 3 complete.")
message("  Histogram   : ", hist_png)
message("  Peak BED    : ", peaks_out)
message("  Count matrix: ", counts_out)
