#!/usr/bin/env Rscript
# =============================================================================
#  Step 8b: fgsea-like peak-set enrichment for ATAC-seq DAR results
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fgsea)
  library(ggplot2)
})

dir_proj <- Sys.getenv("PIPELINE_DIR", unset = getwd())
peak_dir <- Sys.getenv("PIPELINE_PEAK_DIR", unset = file.path(dir_proj, "Peak_nomodel"))
manifest_file <- Sys.getenv("PIPELINE_PEAK_SET_MANIFEST", unset = file.path(peak_dir, "PeakSetLibrary", "region_sets_manifest.tsv"))
results_dir <- Sys.getenv("PIPELINE_FGSEA_RESULTS_DIR", unset = file.path(peak_dir, "PeakSetEnrichment"))
fgsea_min <- as.integer(Sys.getenv("PIPELINE_FGSEA_MIN_SIZE", unset = "15"))
fgsea_max <- as.integer(Sys.getenv("PIPELINE_FGSEA_MAX_SIZE", unset = "500"))
fgsea_padj <- as.numeric(Sys.getenv("PIPELINE_FGSEA_PADJ", unset = "0.05"))

make_peak_id <- function(chr, start, end) {
  paste0(chr, ":", start, "-", end)
}

read_region_set <- function(bed_path) {
  dt <- fread(bed_path, sep = "\t", header = FALSE, select = 1:3,
              col.names = c("chr", "start", "end"))
  unique(make_peak_id(dt$chr, dt$start, dt$end))
}

script_dir <- tryCatch({
  normalizePath(dirname(sys.frame(1)$ofile), mustWork = TRUE)
}, error = function(e) {
  normalizePath(file.path(Sys.getenv("PIPELINE_ROOT", unset = file.path(dir_proj, "ATACseq_pipeline_takubo")), "scripts"),
                mustWork = TRUE)
})

source(file.path(script_dir, "plot_utils.R"))
pcfg <- load_plot_config(Sys.getenv("PIPELINE_ROOT", unset = file.path(dir_proj, "ATACseq_pipeline_takubo")))

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(manifest_file, sep = "\t")
manifest <- manifest %>%
  mutate(
    pathway = paste(source, set_name, sep = "::"),
    peak_count = as.integer(peak_count)
  ) %>%
  filter(peak_count > 0)

if (nrow(manifest) == 0) {
  stop("No normalized peak sets found in ", manifest_file)
}

pathways <- lapply(manifest$normalized_bed, read_region_set)
names(pathways) <- manifest$pathway
manifest$normalized_peak_count <- lengths(pathways)

edgeR_dir <- file.path(peak_dir, "edgeR")
dar_files <- list.files(edgeR_dir, pattern = "_DA_edgeR_Results\\.tsv$", full.names = TRUE)
if (length(dar_files) == 0) {
  stop("No edgeR DAR result files found in ", edgeR_dir)
}

build_summary_plot <- function(res_tbl, comparison_name, out_dir, cfg, padj_cutoff) {
  sig_tbl <- res_tbl %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    mutate(direction_class = if_else(NES >= 0, "positive", "negative"))

  if (nrow(sig_tbl) == 0) {
    return(invisible(NULL))
  }

  top_n <- cfg$fgsea$top_n_plot %||% 20
  top_pos <- sig_tbl %>% filter(direction_class == "positive") %>% arrange(padj) %>% head(top_n %/% 2)
  top_neg <- sig_tbl %>% filter(direction_class == "negative") %>% arrange(padj) %>% head(top_n %/% 2)
  plot_tbl <- bind_rows(top_pos, top_neg) %>%
    arrange(NES) %>%
    mutate(pathway_label = paste(source, set_name, sep = " | "))

  if (nrow(plot_tbl) == 0) {
    return(invisible(NULL))
  }

  volcano_cols <- cfg$colors$volcano
  p <- ggplot(plot_tbl, aes(x = reorder(pathway_label, NES), y = NES,
                            fill = direction_class)) +
    geom_col(width = cfg$fgsea$bar_width %||% 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(
      positive = unname(volcano_cols$up %||% "#E64B35"),
      negative = unname(volcano_cols$down %||% "#4DBBD5")
    ), guide = "none") +
    labs(
      title = paste0("Peak-set enrichment: ", comparison_name),
      subtitle = paste0("padj < ", padj_cutoff),
      x = NULL,
      y = "Normalized Enrichment Score (NES)"
    ) +
    theme_pipeline(cfg) +
    theme(axis.text.y = element_text(size = 8))

  save_plot(
    p,
    paste0("peakset_fgsea_barplot_", comparison_name),
    type = "fgsea_bar",
    outdir = out_dir,
    cfg = cfg
  )
}

all_results <- vector("list", length(dar_files))

for (i in seq_along(dar_files)) {
  dar_file <- dar_files[[i]]
  comparison_name <- sub("_DA_edgeR_Results\\.tsv$", "", basename(dar_file))
  comparison_dir <- file.path(results_dir, comparison_name)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

  dar_tbl <- fread(dar_file, sep = "\t") %>%
    mutate(
      peak_id = make_peak_id(chr, start, end)
    )

  stats <- dar_tbl$logFC
  names(stats) <- dar_tbl$peak_id
  stats <- stats[!is.na(stats)]
  stats <- sort(stats, decreasing = TRUE)

  fgsea_res <- fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    eps = 0,
    minSize = fgsea_min,
    maxSize = fgsea_max
  )

  manifest_join <- manifest[, setdiff(names(manifest), "comparison"), with = FALSE]
  res_tbl <- as_tibble(fgsea_res) %>%
    mutate(
      comparison = comparison_name,
      leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1)),
      leadingEdgeSize = lengths(fgsea_res$leadingEdge)
    ) %>%
    left_join(as_tibble(manifest_join), by = "pathway") %>%
    arrange(padj, desc(abs(NES)))

  fwrite(res_tbl, file.path(comparison_dir, paste0(comparison_name, "_peakset_fgsea.tsv")), sep = "\t")

  for (source_name in unique(res_tbl$source)) {
    source_tbl <- res_tbl %>% filter(source == source_name)
    fwrite(
      source_tbl,
      file.path(comparison_dir, paste0(comparison_name, "_peakset_fgsea_", source_name, ".tsv")),
      sep = "\t"
    )
  }

  build_summary_plot(res_tbl, comparison_name, comparison_dir, pcfg, fgsea_padj)
  all_results[[i]] <- as.data.frame(res_tbl)
}

result_tables <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, all_results)

if (length(result_tables) == 0) {
  stop("No peak-set enrichment result tables were generated.")
}

combined_df <- bind_rows(result_tables) %>%
  arrange(comparison, padj, desc(abs(NES)))
fwrite(combined_df, file.path(results_dir, "peakset_fgsea_all_comparisons.tsv"), sep = "\t")

summary_tbl <- combined_df %>%
  group_by(comparison, source) %>%
  summarise(
    tested_sets = n(),
    significant_sets = sum(!is.na(padj) & padj < fgsea_padj),
    .groups = "drop"
  )
fwrite(summary_tbl, file.path(results_dir, "peakset_fgsea_summary.tsv"), sep = "\t")

message("Step 8 complete. Results: ", results_dir)