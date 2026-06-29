options(stringsAsFactors = FALSE)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) == 0) return(invisible(NULL))

  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr", "multcompView", "circlize"))
install_if_missing(c("ComplexHeatmap"), bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(multcompView)
  library(circlize)
  library(ComplexHeatmap)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

group_levels <- c("5503_sandy", "5503_clay", "6007_sandy", "6007_clay")

category_colors <- c(
  "Environment-sensitive in 5503" = "#d95f02",
  "Environment-sensitive in 6007" = "#e41a1c",
  "Stable in 5503 across soils" = "#7570b3",
  "Stable in 6007 across soils" = "#e7298a"
)

category_levels <- names(category_colors)

fungi_project_dir <- "C:/Users/joyde/OneDrive/Área de Trabalho/Desktop NIOO/Amplicons/ITS_analysis/funfun"

out_dir <- file.path(project_dir, "thesis_heatmap_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

read_table_safely <- function(path) {
  stop_if_missing(path)
  first_line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  delim <- if (grepl(";", first_line, fixed = TRUE) && !grepl(",", first_line, fixed = TRUE)) ";" else ","
  readr::read_delim(path, delim = delim, show_col_types = FALSE, trim_ws = TRUE)
}

normalize_pair_key <- function(g1, g2) {
  paste(sort(c(g1, g2)), collapse = "__")
}

make_letters_generic <- function(df_mod, feature_col) {
  groups_present <- unique(c(df_mod$group1, df_mod$group2))
  groups_present <- intersect(group_levels, groups_present)

  if (length(groups_present) == 0) {
    return(tibble(GroupLabel = group_levels, Letters = ""))
  }

  if (length(groups_present) == 1) {
    out <- tibble(GroupLabel = group_levels, Letters = "")
    out$Letters[out$GroupLabel == groups_present[1]] <- "a"
    return(out)
  }

  pmat <- matrix(1, nrow = length(groups_present), ncol = length(groups_present))
  rownames(pmat) <- groups_present
  colnames(pmat) <- groups_present

  for (i in seq_len(nrow(df_mod))) {
    g1 <- df_mod$group1[i]
    g2 <- df_mod$group2[i]
    pval <- if ("p.adj" %in% names(df_mod) && !is.na(df_mod$`p.adj`[i])) df_mod$`p.adj`[i] else df_mod$p[i]
    if (is.na(g1) || is.na(g2) || is.na(pval)) next
    if (g1 %in% groups_present && g2 %in% groups_present) {
      pmat[g1, g2] <- pval
      pmat[g2, g1] <- pval
    }
  }

  letters_obj <- multcompView::multcompLetters(pmat, threshold = 0.05)
  out <- tibble(GroupLabel = group_levels, Letters = "")
  idx <- match(names(letters_obj$Letters), out$GroupLabel)
  out$Letters[idx[!is.na(idx)]] <- unname(letters_obj$Letters)
  out
}

build_letters_table <- function(selected_df, dunn_df, feature_col, label_col) {
  feats <- unique(selected_df[[label_col]])
  bind_rows(lapply(feats, function(feat_label) {
    df_mod <- dunn_df %>% filter(.data[[label_col]] == feat_label)
    letters_df <- make_letters_generic(df_mod, feature_col)
    tibble(FeatureLabel = feat_label, GroupLabel = letters_df$GroupLabel, Letters = letters_df$Letters)
  }))
}

calc_clr_limit <- function(...) {
  mats <- list(...)
  vals <- unlist(lapply(mats, function(x) as.numeric(x)), use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(2)
  lim <- stats::quantile(abs(vals), probs = 0.98, na.rm = TRUE)
  lim <- max(1, min(4, as.numeric(lim)))
  ceiling(lim * 2) / 2
}

clean_bacteria_label <- function(module, module_name) {
  x <- ifelse(is.na(module_name) | module_name == "", module, paste0(module, " | ", module_name))
  x <- sub(" =>.*$", "", x)
  x <- sub(",.*$", "", x)
  x
}

clean_fungi_label <- function(pathway_label) {
  x <- pathway_label
  x <- sub("^(ko\\d{5}) \\| \\d{5}\\s+", "\\1 | ", x)
  x <- sub(",.*$", "", x)
  x
}

build_bacteria_prep <- function(compartment = c("roots", "stalk")) {
  compartment <- match.arg(compartment)

  means_file <- file.path(project_dir, "kegg_modules_analysis", "tables", "module_group_means_clr.csv")
  dunn_file <- file.path(project_dir, "kegg_modules_analysis", "posthoc",
                         if (compartment == "roots") "dunn_kegg_modules_roots.csv" else "dunn_kegg_modules_stalk.csv")
  selected_file <- file.path(
    project_dir, "kegg_modules_analysis", "final_figure", "final_plot_priority_env_first", "tables",
    if (compartment == "roots") "roots_top10_modules_used_in_heatmap.csv" else "stalk_top10_modules_used_in_heatmap.csv"
  )

  stop_if_missing(means_file)
  stop_if_missing(dunn_file)
  stop_if_missing(selected_file)

  means_df <- read_table_safely(means_file)
  dunn_df <- read_table_safely(dunn_file)
  selected_df <- read_table_safely(selected_file) %>%
    mutate(
      FeatureLabel = clean_bacteria_label(Module, Module_Name),
      Category = factor(Category, levels = category_levels)
    ) %>%
    arrange(Category, Kruskal_rank_value, desc(CLR_Variance))

  pattern <- if (compartment == "roots") "root" else "stalk"

  means_sub <- means_df %>%
    filter(str_detect(SampleTypeClean, pattern)) %>%
    filter(Module %in% selected_df$Module) %>%
    mutate(
      FeatureLabel = clean_bacteria_label(Module, Module_Name),
      GroupLabel = factor(GroupLabel, levels = group_levels)
    ) %>%
    select(Module, FeatureLabel, GroupLabel, Mean_CLR) %>%
    distinct()

  mat <- means_sub %>%
    select(FeatureLabel, GroupLabel, Mean_CLR) %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- mat$FeatureLabel
  mat$FeatureLabel <- NULL
  mat <- as.matrix(mat)

  ann_df <- selected_df %>%
    select(Module, FeatureLabel, Category, CLR_Variance) %>%
    distinct() %>%
    filter(FeatureLabel %in% rownames(mat))
  ann_df <- ann_df[match(rownames(mat), ann_df$FeatureLabel), , drop = FALSE]

  dunn_df <- dunn_df %>%
    mutate(FeatureLabel = clean_bacteria_label(Module, Module_Name))

  letters_mat <- build_letters_table(ann_df, dunn_df, "Module", "FeatureLabel") %>%
    pivot_wider(names_from = GroupLabel, values_from = Letters, values_fill = "") %>%
    as.data.frame()
  rownames(letters_mat) <- letters_mat$FeatureLabel
  letters_mat$FeatureLabel <- NULL
  letters_mat <- as.matrix(letters_mat)
  letters_mat <- letters_mat[rownames(mat), colnames(mat), drop = FALSE]

  list(mat = mat, ann_df = ann_df, letters_mat = letters_mat)
}

build_fungi_prep <- function(compartment = c("roots", "stalk")) {
  compartment <- match.arg(compartment)

  means_file <- file.path(fungi_project_dir, "funfun_results", "funfun_kegg_pathway_analysis", "tables", "funfun_pathway_group_means_clr.csv")
  dunn_file <- file.path(fungi_project_dir, "funfun_results", "funfun_kegg_pathway_analysis", "posthoc",
                         if (compartment == "roots") "dunn_funfun_pathways_roots.csv" else "dunn_funfun_pathways_stalk.csv")
  selected_file <- file.path(
    fungi_project_dir, "funfun_results", "funfun_kegg_pathway_analysis", "exploratory_classification", "tables",
    if (compartment == "roots") "roots_top10_pathways_heatmap.csv" else "stalk_top10_pathways_heatmap.csv"
  )

  stop_if_missing(means_file)
  stop_if_missing(dunn_file)
  stop_if_missing(selected_file)

  means_df <- read_table_safely(means_file)
  dunn_df <- read_table_safely(dunn_file)
  selected_df <- read_table_safely(selected_file) %>%
    mutate(
      FeatureLabel = clean_fungi_label(PathwayLabel),
      Category = factor(Category, levels = category_levels)
    ) %>%
    arrange(Category, desc(CLR_Variance))

  pattern <- if (compartment == "roots") "root" else "stalk"

  means_sub <- means_df %>%
    filter(str_detect(SampleTypeClean, pattern)) %>%
    filter(PathwayLabel %in% selected_df$PathwayLabel) %>%
    mutate(
      FeatureLabel = clean_fungi_label(PathwayLabel),
      GroupLabel = factor(GroupLabel, levels = group_levels)
    ) %>%
    select(PathwayLabel, FeatureLabel, GroupLabel, Mean_CLR) %>%
    distinct()

  mat <- means_sub %>%
    select(FeatureLabel, GroupLabel, Mean_CLR) %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame()
  rownames(mat) <- mat$FeatureLabel
  mat$FeatureLabel <- NULL
  mat <- as.matrix(mat)

  ann_df <- selected_df %>%
    select(PathwayLabel, FeatureLabel, Category, CLR_Variance) %>%
    distinct() %>%
    filter(FeatureLabel %in% rownames(mat))
  ann_df <- ann_df[match(rownames(mat), ann_df$FeatureLabel), , drop = FALSE]

  dunn_df <- dunn_df %>%
    mutate(FeatureLabel = clean_fungi_label(PathwayLabel))

  letters_mat <- build_letters_table(ann_df, dunn_df, "PathwayLabel", "FeatureLabel") %>%
    pivot_wider(names_from = GroupLabel, values_from = Letters, values_fill = "") %>%
    as.data.frame()
  rownames(letters_mat) <- letters_mat$FeatureLabel
  letters_mat$FeatureLabel <- NULL
  letters_mat <- as.matrix(letters_mat)
  letters_mat <- letters_mat[rownames(mat), colnames(mat), drop = FALSE]

  list(mat = mat, ann_df = ann_df, letters_mat = letters_mat)
}

make_single_heatmap <- function(prepped, compartment_title = NULL, clr_limit, row_fontsize = 12, letter_fontsize = 10) {
  mat <- prepped$mat
  ann_df <- prepped$ann_df
  letters_mat <- prepped$letters_mat
  present_categories <- unique(as.character(ann_df$Category))

  left_ha <- rowAnnotation(
    Category = ann_df$Category,
    col = list(Category = category_colors),
    show_annotation_name = FALSE,
    width = grid::unit(6, "mm")
  )

  hm <- Heatmap(
    mat,
    name = "Mean CLR",
    col = circlize::colorRamp2(c(-clr_limit, 0, clr_limit), c("#2166ac", "white", "#b2182b")),
    width = grid::unit(5.8, "cm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = row_fontsize),
    row_names_max_width = grid::unit(7.5, "cm"),
    column_names_gp = grid::gpar(fontsize = 11),
    row_split = factor(ann_df$Category, levels = present_categories),
    row_gap = grid::unit(2, "mm"),
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title = compartment_title,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    left_annotation = left_ha,
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- letters_mat[i, j]
      if (!is.na(lab) && nzchar(lab)) {
        grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = letter_fontsize, col = "black"))
      }
    }
  )

  hm
}

draw_panel <- function(root_prep, stalk_prep, figure_title = NULL, file_stub) {
  clr_limit <- calc_clr_limit(root_prep$mat, stalk_prep$mat)
  is_fungal_panel <- identical(file_stub, "fungal_thesis_panel")

  ht_root <- make_single_heatmap(root_prep, NULL, clr_limit)
  ht_stalk <- make_single_heatmap(stalk_prep, NULL, clr_limit)
  ht <- ht_root %v% ht_stalk

  lgd_category <- Legend(
    title = "Category",
    at = names(category_colors),
    legend_gp = grid::gpar(fill = unname(category_colors)),
    labels_gp = grid::gpar(fontsize = 10),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold")
  )

  lgd_clr <- Legend(
    title = "Mean CLR",
    col_fun = circlize::colorRamp2(c(-clr_limit, 0, clr_limit), c("#2166ac", "white", "#b2182b")),
    at = c(-clr_limit, 0, clr_limit),
    labels_gp = grid::gpar(fontsize = 10),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold")
  )

  pd <- NULL
  if (!is_fungal_panel) {
    pd <- packLegend(
      lgd_category,
      lgd_clr,
      direction = "vertical",
      gap = grid::unit(5, "mm")
    )
  }

  png(file.path(out_dir, paste0(file_stub, ".png")), width = 3200, height = 5200, res = 300)
  draw(
    ht,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    padding = grid::unit(c(8, 12, 28, 6), "mm")
  )
  if (!is.null(figure_title) && nzchar(figure_title)) {
    grid::grid.text(
      figure_title,
      x = grid::unit(0.03, "npc"),
      y = grid::unit(0.985, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 15, fontface = "bold")
    )
  }
  if (is_fungal_panel) {
    draw(
      lgd_clr,
      x = grid::unit(0.04, "npc"),
      y = grid::unit(0.92, "npc"),
      just = c("left", "top")
    )
    draw(
      lgd_category,
      x = grid::unit(0.04, "npc"),
      y = grid::unit(0.02, "npc"),
      just = c("left", "bottom")
    )
  } else {
    draw(
      pd,
      x = grid::unit(0.985, "npc"),
      y = grid::unit(0.03, "npc"),
      just = c("right", "bottom")
    )
  }
  dev.off()

  pdf(file.path(out_dir, paste0(file_stub, ".pdf")), width = 8.2, height = 11.4)
  draw(
    ht,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    padding = grid::unit(c(8, 12, 28, 6), "mm")
  )
  if (!is.null(figure_title) && nzchar(figure_title)) {
    grid::grid.text(
      figure_title,
      x = grid::unit(0.03, "npc"),
      y = grid::unit(0.985, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 15, fontface = "bold")
    )
  }
  if (is_fungal_panel) {
    draw(
      lgd_clr,
      x = grid::unit(0.04, "npc"),
      y = grid::unit(0.92, "npc"),
      just = c("left", "top")
    )
    draw(
      lgd_category,
      x = grid::unit(0.04, "npc"),
      y = grid::unit(0.02, "npc"),
      just = c("left", "bottom")
    )
  } else {
    draw(
      pd,
      x = grid::unit(0.985, "npc"),
      y = grid::unit(0.03, "npc"),
      just = c("right", "bottom")
    )
  }
  dev.off()
}

bacteria_roots <- build_bacteria_prep("roots")
bacteria_stalk <- build_bacteria_prep("stalk")
fungi_roots <- build_fungi_prep("roots")
fungi_stalk <- build_fungi_prep("stalk")

draw_panel(
  bacteria_roots,
  bacteria_stalk,
  NULL,
  "bacterial_thesis_panel"
)

draw_panel(
  fungi_roots,
  fungi_stalk,
  NULL,
  "fungal_thesis_panel"
)

message("")
message("Painéis gerados com sucesso.")
message("Saídas em: ", out_dir)
message("- ", file.path(out_dir, "bacterial_thesis_panel.png"))
message("- ", file.path(out_dir, "bacterial_thesis_panel.pdf"))
message("- ", file.path(out_dir, "fungal_thesis_panel.png"))
message("- ", file.path(out_dir, "fungal_thesis_panel.pdf"))
