# ============================================================
# plot_utils.R — 可視化設定の読み込みユーティリティ
# ============================================================
# 使用方法:
#   source("scripts/plot_utils.R")
#   cfg <- load_plot_config()
#   cfg$colors$group_palette   # カラーパレット取得
#   theme_pipeline()           # 統一テーマを適用
# ============================================================

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

#' プロット設定ファイルを読み込む
#'
#' plot_config.yml (個人設定) → plot_config.default.yml (デフォルト) の順に探索。
#' 個人設定が存在すればそちらを優先し、存在しなければデフォルトを使用する。
#'
#' @param config_dir 設定ファイルがあるディレクトリ (デフォルト: プロジェクトルート)
#' @return list  設定内容
load_plot_config <- function(config_dir = NULL) {

  # プロジェクトルートを推定
  if (is.null(config_dir)) {
    # スクリプトの場所 or 作業ディレクトリから探索
    candidates <- c(
      Sys.getenv("PIPELINE_ROOT", unset = ""),
      getwd(),
      dirname(sys.frame(1)$ofile %||% ".")
    )
    for (d in candidates) {
      if (d != "" && file.exists(file.path(d, "plot_config.default.yml"))) {
        config_dir <- d
        break
      }
    }
    if (is.null(config_dir)) config_dir <- getwd()
  }

  user_file    <- file.path(config_dir, "plot_config.yml")
  default_file <- file.path(config_dir, "plot_config.default.yml")

  if (file.exists(user_file)) {
    message("[plot_utils] 個人設定を読み込み: ", user_file)
    cfg <- yaml::read_yaml(user_file)
  } else if (file.exists(default_file)) {
    message("[plot_utils] デフォルト設定を読み込み: ", default_file)
    message("[plot_utils] 個人カスタマイズ: cp plot_config.default.yml plot_config.yml")
    cfg <- yaml::read_yaml(default_file)
  } else {
    warning("[plot_utils] 設定ファイルが見つかりません。ハードコード値を使用します。")
    cfg <- .default_config()
  }

  # 構造を検証して欠損キーをデフォルトで補完
  cfg <- .fill_defaults(cfg)

  return(cfg)
}


#' パイプライン統一 ggplot2 テーマ
#'
#' @param cfg 設定リスト (省略時は自動読み込み)
#' @return ggplot2 theme オブジェクト
theme_pipeline <- function(cfg = NULL) {
  if (is.null(cfg)) cfg <- load_plot_config()

  base_size   <- cfg$theme$base_size   %||% 12
  base_family <- cfg$theme$base_family %||% "sans"
  theme_name  <- cfg$theme$ggplot_theme %||% "theme_bw"

  # テーマ関数を取得
  theme_fn <- tryCatch(
    get(theme_name, envir = asNamespace("ggplot2")),
    error = function(e) ggplot2::theme_bw
  )

  theme_fn(base_size = base_size, base_family = base_family) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      strip.text   = element_text(face = "bold"),
      legend.position = "right"
    )
}


#' グループカラーを取得
#'
#' group_colors (名前付き) が定義されていればそちらを、
#' なければ group_palette から順に割り当てる。
#'
#' @param groups  グループ名のベクトル
#' @param cfg     設定リスト
#' @return named character vector  グループ名 → 色
get_group_colors <- function(groups, cfg = NULL) {
  if (is.null(cfg)) cfg <- load_plot_config()

  groups <- unique(groups)

  # 名前付き色が定義されている場合
  if (!is.null(cfg$colors$group_colors)) {
    named <- cfg$colors$group_colors
    # 足りない分はパレットから補完
    missing <- setdiff(groups, names(named))
    if (length(missing) > 0) {
      palette <- cfg$colors$group_palette
      extra <- palette[seq_along(missing)]
      names(extra) <- missing
      named <- c(named, extra)
    }
    return(unlist(named[groups]))
  }

  # パレットから割り当て
  palette <- cfg$colors$group_palette
  if (length(groups) > length(palette)) {
    palette <- rep_len(palette, length(groups))
  }
  cols <- palette[seq_along(groups)]
  names(cols) <- groups
  return(unlist(cols))
}


#' 図を保存するラッパー
#'
#' @param plot     ggplot オブジェクト
#' @param filename ファイル名 (拡張子なし)
#' @param type     図の種類 ("pca", "volcano", "heatmap", "barplot")
#' @param outdir   出力ディレクトリ
#' @param cfg      設定リスト
save_plot <- function(plot, filename, type = "pca", outdir = ".", cfg = NULL) {
  if (is.null(cfg)) cfg <- load_plot_config()

  fmt    <- cfg$output$format %||% "png"
  dpi    <- cfg$output$dpi    %||% 300
  size   <- cfg$figure_size[[type]] %||% list(width = 8, height = 6)
  width  <- size$width
  height <- size$height

  filepath <- file.path(outdir, paste0(filename, ".", fmt))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  ggsave(
    filename = filepath,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    units    = "in"
  )
  message("[plot_utils] 保存: ", filepath, " (", width, "x", height, " in, ", dpi, " dpi)")
  return(invisible(filepath))
}


# ----------------------------------------------------------
# 内部ヘルパー
# ----------------------------------------------------------

#' %||% 演算子 (NULL合体)
`%||%` <- function(x, y) if (is.null(x)) y else x

#' ハードコードのフォールバック設定
.default_config <- function() {
  list(
    theme = list(base_family = "sans", base_size = 12, ggplot_theme = "theme_bw"),
    figure_size = list(
      pca     = list(width = 8, height = 6),
      volcano = list(width = 8, height = 6),
      heatmap = list(width = 10, height = 8),
      barplot = list(width = 10, height = 6)
    ),
    output = list(format = "png", dpi = 300),
    colors = list(
      group_palette = c("#E64B35","#4DBBD5","#00A087","#3C5488",
                        "#F39B7F","#8491B4","#91D1C2","#DC9100"),
      heatmap = list(low = "#313695", mid = "#FFFFBF", high = "#A50026"),
      volcano = list(up = "#E64B35", down = "#4DBBD5", ns = "grey70")
    ),
    pca     = list(point_size = 4, label_size = 3, show_labels = TRUE, ellipse = TRUE, ellipse_alpha = 0.1),
    volcano = list(point_size = 1.5, label_top_n = 15, alpha_ns = 0.3),
    heatmap = list(show_rownames = FALSE, show_colnames = TRUE, cluster_rows = TRUE,
                   cluster_cols = TRUE, scale = "row", top_n = 50)
  )
}

#' 欠損キーをデフォルトで補完
.fill_defaults <- function(cfg) {
  defaults <- .default_config()

  fill_recursive <- function(user, def) {
    for (key in names(def)) {
      if (is.null(user[[key]])) {
        user[[key]] <- def[[key]]
      } else if (is.list(def[[key]]) && is.list(user[[key]])) {
        user[[key]] <- fill_recursive(user[[key]], def[[key]])
      }
    }
    return(user)
  }

  fill_recursive(cfg, defaults)
}
