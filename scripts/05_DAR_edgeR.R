#!/usr/bin/env Rscript
# =============================================================================
#  Step 5: DAR Detection (edgeR quasi-likelihood)
#
#  処理内容:
#    - PeakCounts_th0.txt を読み込み
#    - samples.tsv からグループ情報を取得
#    - edgeR QL (glmQLFit / glmQLFTest) で全ペアワイズ比較
#    - 結果 TSV と Volcano plot を出力
#
#  環境変数:
#    PIPELINE_DIR, PIPELINE_SAMPLES, PIPELINE_PEAK_DIR,
#    PIPELINE_DAR_FDR, PIPELINE_DAR_LFC
#
#  出力: Peak_nomodel/edgeR/
#    {G1}_vs_{G2}_DA_edgeR_Results.tsv
#    {G1}_vs_{G2}_edgeR_volcano.png
# =============================================================================
suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(tidyverse)
  library(data.table)
})

# --- 環境変数 ---
dir_proj    <- Sys.getenv("PIPELINE_DIR",      unset = getwd())
samples_tsv <- Sys.getenv("PIPELINE_SAMPLES",  unset = file.path(dir_proj, "samples.tsv"))
peak_dir    <- Sys.getenv("PIPELINE_PEAK_DIR", unset = file.path(dir_proj, "Peak_nomodel"))
dar_fdr     <- as.numeric(Sys.getenv("PIPELINE_DAR_FDR", unset = "0.05"))
dar_lfc     <- as.numeric(Sys.getenv("PIPELINE_DAR_LFC", unset = "1"))

out_dir <- file.path(peak_dir, "edgeR")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- サンプルシート ---
samples_df <- read.table(samples_tsv, header = TRUE, sep = "\t",
                         comment.char = "#", stringsAsFactors = FALSE)
colnames(samples_df) <- c("sample_name", "group")

# --- カウントデータ読み込み ---
counts_file <- file.path(peak_dir, "PeakCounts_th0.txt")
stopifnot(file.exists(counts_file))
d <- fread(counts_file, sep = "\t")

# 染色体・座標列を取得
coord_cols <- c("chr", "start", "end")
count_cols  <- samples_df$sample_name

# region_name (行名) = "chr_start_end"
region_name <- paste(d$chr, d$start, d$end, sep = "_")

# カウント行列 (サンプルシートと同順に並べ直し)
raw_mat <- as.matrix(d[, ..count_cols])
rownames(raw_mat) <- region_name

# グループファクター (samples.tsv の出現順を保持)
group_levels <- unique(samples_df$group)
class <- factor(samples_df$group, levels = group_levels)
message("Groups: ", paste(group_levels, collapse = " / "))
message("Samples: ", paste(samples_df$sample_name, collapse = ", "))

# =============================================================================
#  edgeR QL 統計モデル作成
# =============================================================================
design <- model.matrix(~ 0 + class)
colnames(design) <- levels(class)

y <- DGEList(counts = raw_mat, group = class)
keep <- filterByExpr(y, group = class)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
kept_rows <- rownames(y)

message("Peaks after filterByExpr: ", length(kept_rows))

# =============================================================================
#  全ペアワイズ比較
# =============================================================================
pairs <- t(combn(colnames(design), 2))

write_results <- function(tt, g1, g2) {
  # 座標を chr/start/end に分割
  pos_mat <- str_split_fixed(rownames(tt), "_", n = 3)
  df_out <- data.frame(
    chr   = pos_mat[, 1],
    start = as.integer(pos_mat[, 2]),
    end   = as.integer(pos_mat[, 3]),
    tt,
    stringsAsFactors = FALSE
  ) %>%
    arrange(FDR)

  # TSV 出力
  tsv_file <- file.path(out_dir, paste0(g1, "_vs_", g2, "_DA_edgeR_Results.tsv"))
  fwrite(df_out, tsv_file, sep = "\t")

  # --- Volcano plot ---
  df_vol <- df_out %>%
    mutate(sig = case_when(
      FDR < dar_fdr & logFC >  dar_lfc ~ "up",
      FDR < dar_fdr & logFC < -dar_lfc ~ "down",
      TRUE ~ "ns"
    ))
  n_up   <- sum(df_vol$sig == "up")
  n_down <- sum(df_vol$sig == "down")

  p <- ggplot(df_vol, aes(x = logFC, y = -log10(FDR + 1e-300), fill = sig)) +
    geom_point(pch = 21, size = 1.5, alpha = 0.7, stroke = 0.2) +
    scale_fill_manual(values = c("up" = "#F73719", "ns" = "#CACACA", "down" = "#33DFD4"),
                      labels = c(paste0("up (", n_up, ")"),
                                 "ns",
                                 paste0("down (", n_down, ")"))) +
    geom_vline(xintercept = c(-dar_lfc, dar_lfc), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(dar_fdr), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = dplyr::filter(df_vol, FDR < dar_fdr / 10, abs(logFC) > dar_lfc),
      aes(label = paste0(chr, ":", start, "-", end)),
      size = 2.5, max.overlaps = 8, segment.size = 0.2
    ) +
    labs(title = paste0(g1, " vs ", g2, "  (QLF)"),
         subtitle = paste0("FDR<", dar_fdr, " & |logFC|>", dar_lfc,
                           ": up=", n_up, ", down=", n_down),
         x = paste0("log2FC (", g1, " / ", g2, ")"),
         y = "-log10(FDR)",
         fill = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")

  png_file <- file.path(out_dir, paste0(g1, "_vs_", g2, "_edgeR_volcano.png"))
  ggsave(png_file, p, width = 6, height = 6, dpi = 300)
  message(sprintf("  %s vs %s: up=%d, down=%d  → %s", g1, g2, n_up, n_down,
                  basename(tsv_file)))
}

message("Running pairwise QLF tests (", nrow(pairs), " pairs)...")
for (i in seq_len(nrow(pairs))) {
  g1 <- pairs[i, 1]
  g2 <- pairs[i, 2]

  v <- setNames(rep(0, ncol(design)), colnames(design))
  v[g1] <-  1
  v[g2] <- -1

  qlf <- glmQLFTest(fit, contrast = v)
  tt  <- topTags(qlf, n = Inf)$table
  rownames(tt) <- kept_rows[match(rownames(tt), kept_rows)]

  write_results(tt, g1, g2)
}

message("")
message("Step 5 complete. Results: ", out_dir)
