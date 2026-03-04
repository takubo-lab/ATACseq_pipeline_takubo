#!/usr/bin/env Rscript
# =============================================================================
#  Step 7: PCA & Visualization
#
#  処理内容:
#    1. PeakCounts_th0.txt を読み込み、ライブラリサイズ補正 (CPM) + z-score
#    2. 高分散ピーク (上位 5000) を用いて PCA
#       - PC1 vs PC2、PC1 vs PC3 プロット (グループ着色)
#       - 各サンプル名ラベル付き
#    3. サンプル間相関 ヒートマップ (pearson)
#    4. (2グループ以上の場合) 各ペアワイズ比較の DAR ベン図
#
#  環境変数:
#    PIPELINE_DIR, PIPELINE_SAMPLES, PIPELINE_PEAK_DIR, PIPELINE_PLOTS_DIR,
#    PIPELINE_DAR_FDR, PIPELINE_DAR_LFC
#
#  出力: Plots/
#    PCA_PC1_PC2.png, PCA_PC1_PC3.png
#    Sample_correlation_heatmap.png
#    Scatter_{sampleA}_vs_{sampleB}.png  (各サンプルペア相関散布図)
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(tidyverse)
  library(data.table)
  library(RColorBrewer)
  library(pheatmap)
})

# --- 環境変数 ---
dir_proj   <- Sys.getenv("PIPELINE_DIR",       unset = getwd())
samples_tsv<- Sys.getenv("PIPELINE_SAMPLES",   unset = file.path(dir_proj, "samples.tsv"))
peak_dir   <- Sys.getenv("PIPELINE_PEAK_DIR",  unset = file.path(dir_proj, "Peak_nomodel"))
plots_dir  <- Sys.getenv("PIPELINE_PLOTS_DIR", unset = file.path(dir_proj, "Plots"))
dar_fdr    <- as.numeric(Sys.getenv("PIPELINE_DAR_FDR", unset = "0.05"))
dar_lfc    <- as.numeric(Sys.getenv("PIPELINE_DAR_LFC", unset = "1"))

dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# --- サンプルシート ---
samples_df <- read.table(samples_tsv, header = TRUE, sep = "\t",
                         comment.char = "#", stringsAsFactors = FALSE)
colnames(samples_df) <- c("sample_name", "group")
group_levels <- unique(samples_df$group)
n_groups     <- length(group_levels)

# グループカラーパレット
grp_colors <- setNames(
  colorRampPalette(brewer.pal(min(n_groups, 11), "Set1"))(n_groups),
  group_levels
)

# --- カウントデータ読み込み ---
counts_file <- file.path(peak_dir, "PeakCounts_th0.txt")
stopifnot(file.exists(counts_file))
d_raw <- fread(counts_file, sep = "\t")

# カウント行列: sample_name 列のみ抽出 (samples.tsv と同順)
count_cols   <- samples_df$sample_name
mat_raw      <- as.matrix(d_raw[, ..count_cols])
rownames(mat_raw) <- paste(d_raw$chr, d_raw$start, d_raw$end, sep = "_")

message("Count matrix: ", nrow(mat_raw), " peaks × ", ncol(mat_raw), " samples")

# --- CPM 正規化 ---
lib_sizes <- colSums(mat_raw)
cpm_mat   <- sweep(mat_raw, 2, lib_sizes / 1e6, FUN = "/")

# --- 高分散ピーク (上位 5000) 選択 ---
n_top <- min(5000, nrow(cpm_mat))
row_var <- apply(cpm_mat, 1, var)
top_idx <- order(row_var, decreasing = TRUE)[1:n_top]
cpm_top <- cpm_mat[top_idx, , drop = FALSE]

# --- z-score 正規化 ---
zscore_mat <- t(scale(t(cpm_top)))  # 行方向 (ピークごと) にスケーリング
zscore_mat <- zscore_mat[complete.cases(zscore_mat), ]

message("Features for PCA: ", nrow(zscore_mat))

# =============================================================================
#  PCA
# =============================================================================
pca_res <- prcomp(t(zscore_mat), center = FALSE, scale. = FALSE)
pct_var <- round(100 * summary(pca_res)$importance[2, ], 1)

pca_df <- data.frame(
  pca_res$x[, 1:min(3, ncol(pca_res$x))],
  sample = samples_df$sample_name,
  group  = factor(samples_df$group, levels = group_levels)
)

make_pca_plot <- function(df, xpc, ypc, colors) {
  xvar <- pct_var[xpc]
  yvar <- pct_var[ypc]
  x_col <- paste0("PC", xpc)
  y_col <- paste0("PC", ypc)
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                 fill = group, label = sample)) +
    geom_point(size = 6, pch = 21, alpha = 0.85, stroke = 0.4, color = "black") +
    geom_text_repel(size = 3, max.overlaps = 20, segment.size = 0.3) +
    scale_fill_manual(values = colors) +
    labs(x = paste0("PC", xpc, ": ", xvar, "%"),
         y = paste0("PC", ypc, ": ", yvar, "%"),
         fill = "Group") +
    coord_fixed() +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          panel.grid = element_blank())
}

p12 <- make_pca_plot(pca_df, 1, 2, grp_colors)
ggsave(file.path(plots_dir, "PCA_PC1_PC2.png"), p12, width = 6, height = 6, dpi = 300)
message("Saved: PCA_PC1_PC2.png")

if (ncol(pca_res$x) >= 3) {
  p13 <- make_pca_plot(pca_df, 1, 3, grp_colors)
  ggsave(file.path(plots_dir, "PCA_PC1_PC3.png"), p13, width = 6, height = 6, dpi = 300)
  message("Saved: PCA_PC1_PC3.png")
}

# =============================================================================
#  サンプル間相関ヒートマップ
# =============================================================================
message("Generating sample correlation heatmap...")
cor_mat <- cor(cpm_mat, method = "pearson")

# グループアノテーションカラー
annot_df  <- data.frame(Group = samples_df$group, row.names = samples_df$sample_name)
annot_col <- list(Group = grp_colors)

png(file.path(plots_dir, "Sample_correlation_heatmap.png"),
    width = 8, height = 7, units = "in", res = 300)
pheatmap::pheatmap(
  cor_mat,
  annotation_row = annot_df,
  annotation_col = annot_df,
  annotation_colors = annot_col,
  color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  breaks = seq(0.8, 1, length.out = 101),
  display_numbers = TRUE, number_format = "%.3f", fontsize_number = 7,
  cluster_rows = TRUE, cluster_cols = TRUE,
  main = "Sample Pearson Correlation (CPM)"
)
dev.off()
message("Saved: Sample_correlation_heatmap.png")

# =============================================================================
#  サンプル間散布図 (グループ内レプリカ確認用)
# =============================================================================
plot_scatter_log <- function(cpm_mat, s1, s2, out_path) {
  df <- data.frame(x = log2(cpm_mat[, s1] + 1),
                   y = log2(cpm_mat[, s2] + 1))
  r  <- round(cor(df$x, df$y, method = "pearson"), 3)
  p  <- ggplot(df, aes(x, y)) +
    geom_point(pch = 21, fill = "steelblue", alpha = 0.25, size = 0.8) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", r), size = 5) +
    labs(x = paste0("log2(CPM+1)  ", s1),
         y = paste0("log2(CPM+1)  ", s2)) +
    theme_bw(base_size = 11)
  ggsave(out_path, p, width = 5, height = 5, dpi = 200)
}

# グループ内ペア (レプリカ確認)
for (grp in group_levels) {
  samps_in_grp <- samples_df$sample_name[samples_df$group == grp]
  if (length(samps_in_grp) >= 2) {
    pairs_in_grp <- combn(samps_in_grp, 2, simplify = FALSE)
    for (pr in pairs_in_grp) {
      s1 <- pr[1]; s2 <- pr[2]
      out_png <- file.path(plots_dir,
        paste0("Scatter_", s1, "_vs_", s2, ".png"))
      plot_scatter_log(cpm_mat, s1, s2, out_png)
      message("Saved: Scatter_", s1, "_vs_", s2, ".png")
    }
  }
}

# =============================================================================
#  (オプション) DAR ベン図 — ggvenn (インストール済みの場合のみ)
# =============================================================================
if (n_groups >= 2 && requireNamespace("ggvenn", quietly = TRUE)) {
  library(ggvenn)
  edgeR_dir <- file.path(peak_dir, "edgeR")
  tsv_files <- list.files(edgeR_dir, pattern = "_DA_edgeR_Results.tsv", full.names = TRUE)

  if (length(tsv_files) >= 2) {
    message("Generating Venn diagrams for DAR sets...")

    # up DAR sets per comparison
    dar_sets_up   <- list()
    dar_sets_down <- list()

    for (f in tsv_files) {
      comp_name <- sub("_DA_edgeR_Results.tsv", "", basename(f))
      dt <- fread(f, sep = "\t")
      key <- paste(dt$chr, dt$start, dt$end, sep = "_")
      up_keys   <- key[dt$FDR < dar_fdr & dt$logFC >  dar_lfc]
      down_keys <- key[dt$FDR < dar_fdr & dt$logFC < -dar_lfc]
      dar_sets_up[[comp_name]]   <- up_keys
      dar_sets_down[[comp_name]] <- down_keys
    }

    # 2〜4 要素までベン図を描画
    if (length(dar_sets_up) <= 4) {
      p_venn_up <- ggvenn(dar_sets_up, show_percentage = FALSE) +
        ggtitle(paste0("Up-DAR overlap (FDR<", dar_fdr, ", logFC>", dar_lfc, ")"))
      ggsave(file.path(plots_dir, "Venn_up_DAR.png"), p_venn_up,
             width = 7, height = 5, dpi = 300)

      p_venn_dn <- ggvenn(dar_sets_down, show_percentage = FALSE) +
        ggtitle(paste0("Down-DAR overlap (FDR<", dar_fdr, ", logFC<-", dar_lfc, ")"))
      ggsave(file.path(plots_dir, "Venn_down_DAR.png"), p_venn_dn,
             width = 7, height = 5, dpi = 300)
      message("Saved: Venn_up_DAR.png, Venn_down_DAR.png")
    } else {
      message("More than 4 comparisons — skipping Venn diagram (too many sets).")
    }
  }
}

message("")
message("Step 7 complete. Plots: ", plots_dir)
