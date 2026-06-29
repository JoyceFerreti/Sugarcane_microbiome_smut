# ============================================================
# Classificacao exploratoria de pathways do FunFun
# Baseado em:
# - pathways com Kruskal_p <= 0.05
# - Dunn post-hoc para comparacoes par-a-par
# - classificacao em categorias biologicas
# - UpSet plots
# - heatmaps focados com letras
# ============================================================

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

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr", "ComplexUpset", "ggplot2", "multcompView", "circlize"))
install_if_missing(c("ComplexHeatmap"), bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ComplexUpset)
  library(ggplot2)
  library(multcompView)
  library(circlize)
  library(ComplexHeatmap)
})

# -----------------------------
# 1. Configuracao
# -----------------------------
project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

base_dir <- file.path(project_dir, "funfun_kegg_pathway_analysis")
tables_dir <- file.path(base_dir, "tables")
posthoc_dir <- file.path(base_dir, "posthoc")

out_dir <- file.path(base_dir, "exploratory_classification")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

pathway_group_means_file <- file.path(tables_dir, "funfun_pathway_group_means_clr.csv")
kruskal_roots_file <- file.path(tables_dir, "kruskal_funfun_pathways_roots.csv")
kruskal_stalk_file <- file.path(tables_dir, "kruskal_funfun_pathways_stalk.csv")
dunn_roots_file <- file.path(posthoc_dir, "dunn_funfun_pathways_roots.csv")
dunn_stalk_file <- file.path(posthoc_dir, "dunn_funfun_pathways_stalk.csv")

top_n_heatmap <- 10
top_n_table <- 20
group_levels <- c("5503_sandy", "5503_clay", "6007_sandy", "6007_clay")

priority_order <- c(
  "Environment_sensitive_5503",
  "Environment_sensitive_6007",
  "Stable_5503_across_soils",
  "Stable_6007_across_soils"
)

category_labels <- c(
  Environment_sensitive_5503 = "Environment-sensitive in 5503",
  Environment_sensitive_6007 = "Environment-sensitive in 6007",
  Stable_5503_across_soils = "Stable in 5503 across soils",
  Stable_6007_across_soils = "Stable in 6007 across soils"
)

category_colors <- c(
  "Environment-sensitive in 5503" = "#d95f02",
  "Environment-sensitive in 6007" = "#e41a1c",
  "Stable in 5503 across soils" = "#7570b3",
  "Stable in 6007 across soils" = "#e7298a"
)

# -----------------------------
# 2. Funcoes auxiliares
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

for (f in c(pathway_group_means_file, kruskal_roots_file, kruskal_stalk_file, dunn_roots_file, dunn_stalk_file)) {
  stop_if_missing(f)
}

normalize_pair_key <- function(g1, g2) {
  pair <- sort(c(g1, g2))
  paste(pair, collapse = "__")
}

get_pair_status <- function(dunn_df, alpha = 0.05) {
  if (nrow(dunn_df) == 0) return(tibble())

  # Exploratiorio: prioriza p.adj quando disponivel; se ausente, usa p
  dunn_df %>%
    mutate(
      PairKey = mapply(normalize_pair_key, group1, group2),
      P_used = ifelse(!is.na(p.adj), p.adj, p),
      Significant = !is.na(P_used) & P_used <= alpha
    ) %>%
    select(Compartimento, PathwayID, PathwayName, PathwayLabel, PairKey, Significant, P_used) %>%
    distinct()
}

classify_pathways <- function(dunn_df, kruskal_df, compartment_name) {
  pair_df <- get_pair_status(dunn_df, alpha = 0.05)

  sig_paths <- kruskal_df %>%
    filter(!is.na(Kruskal_p) & Kruskal_p <= 0.05) %>%
    select(PathwayID, PathwayName, PathwayLabel) %>%
    distinct()

  if (nrow(sig_paths) == 0) return(tibble())

  expected_pairs <- tibble(
    PairKey = c(
      normalize_pair_key("5503_sandy", "5503_clay"),
      normalize_pair_key("6007_sandy", "6007_clay"),
      normalize_pair_key("5503_sandy", "6007_sandy"),
      normalize_pair_key("5503_clay", "6007_clay")
    ),
    ContrastType = c(
      "same_genotype_5503_across_soils",
      "same_genotype_6007_across_soils",
      "same_soil_sandy_across_genotypes",
      "same_soil_clay_across_genotypes"
    )
  )

  expanded <- sig_paths %>%
    tidyr::crossing(expected_pairs) %>%
    left_join(pair_df, by = c("PathwayID", "PathwayName", "PathwayLabel", "PairKey")) %>%
    mutate(
      Compartimento = compartment_name,
      Significant = dplyr::coalesce(Significant, FALSE),
      NoDifference = !Significant
    )

  expanded %>%
    select(Compartimento, PathwayID, PathwayName, PathwayLabel, ContrastType, Significant, NoDifference) %>%
    pivot_wider(
      names_from = ContrastType,
      values_from = c(Significant, NoDifference),
      values_fill = FALSE
    ) %>%
    mutate(
      Stable_5503_across_soils = NoDifference_same_genotype_5503_across_soils,
      Stable_6007_across_soils = NoDifference_same_genotype_6007_across_soils,
      Environment_sensitive_5503 = Significant_same_genotype_5503_across_soils,
      Environment_sensitive_6007 = Significant_same_genotype_6007_across_soils,
      Common_in_sandy_between_genotypes = NoDifference_same_soil_sandy_across_genotypes,
      Common_in_clay_between_genotypes = NoDifference_same_soil_clay_across_genotypes,
      Genotype_different_in_sandy = Significant_same_soil_sandy_across_genotypes,
      Genotype_different_in_clay = Significant_same_soil_clay_across_genotypes
    )
}

assign_primary_category <- function(class_df) {
  class_df %>%
    rowwise() %>%
    mutate(
      PrimaryCategoryKey = {
        hits <- priority_order[unlist(c_across(all_of(priority_order)))]
        if (length(hits) == 0) NA_character_ else hits[1]
      },
      Category = dplyr::recode(PrimaryCategoryKey, !!!category_labels)
    ) %>%
    ungroup() %>%
    filter(!is.na(Category))
}

prepare_rank_table <- function(class_df, means_df, compartment_pattern) {
  means_sub <- means_df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern))

  variability_df <- means_sub %>%
    group_by(PathwayID, PathwayName, PathwayLabel) %>%
    summarise(CLR_Variance = var(Mean_CLR, na.rm = TRUE), .groups = "drop")

  class_primary <- assign_primary_category(class_df)

  class_primary %>%
    left_join(variability_df, by = c("PathwayID", "PathwayName", "PathwayLabel")) %>%
    mutate(CLR_Variance = ifelse(is.na(CLR_Variance), 0, CLR_Variance)) %>%
    arrange(factor(Category, levels = unname(category_labels)), desc(CLR_Variance))
}

select_top_by_category <- function(rank_df, n) {
  rank_df %>%
    group_by(Category) %>%
    arrange(desc(CLR_Variance), .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup()
}

make_letters_for_pathway <- function(df_mod, group_levels) {
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
    pval <- ifelse(!is.na(df_mod$p.adj[i]), df_mod$p.adj[i], df_mod$p[i])
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

build_letters_table <- function(selected_df, dunn_df) {
  paths <- unique(selected_df$PathwayLabel)

  bind_rows(lapply(paths, function(path_lab) {
    df_mod <- dunn_df %>% filter(PathwayLabel == path_lab)
    letters_df <- make_letters_for_pathway(df_mod, group_levels)
    tibble(PathwayLabel = path_lab, GroupLabel = letters_df$GroupLabel, Letters = letters_df$Letters)
  }))
}

build_heatmap_objects <- function(selected_df, dunn_df, means_df, compartment_pattern) {
  means_sub <- means_df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    filter(PathwayLabel %in% selected_df$PathwayLabel) %>%
    mutate(GroupLabel = factor(GroupLabel, levels = group_levels)) %>%
    select(PathwayLabel, GroupLabel, Mean_CLR) %>%
    distinct()

  mat <- means_sub %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame()

  rownames(mat) <- mat$PathwayLabel
  mat$PathwayLabel <- NULL
  mat <- as.matrix(mat)

  ann_df <- selected_df %>%
    distinct(PathwayID, PathwayName, PathwayLabel, Category, CLR_Variance) %>%
    filter(PathwayLabel %in% rownames(mat))

  ann_df <- ann_df[match(rownames(mat), ann_df$PathwayLabel), , drop = FALSE]

  letters_long <- build_letters_table(selected_df, dunn_df)
  letters_mat <- letters_long %>%
    pivot_wider(names_from = GroupLabel, values_from = Letters, values_fill = "") %>%
    as.data.frame()

  rownames(letters_mat) <- letters_mat$PathwayLabel
  letters_mat$PathwayLabel <- NULL
  letters_mat <- as.matrix(letters_mat)
  letters_mat <- letters_mat[rownames(mat), colnames(mat), drop = FALSE]

  list(mat = mat, ann_df = ann_df, letters_mat = letters_mat)
}

plot_focus_heatmap <- function(prepped, title_txt, file_name) {
  mat <- prepped$mat
  ann_df <- prepped$ann_df
  letters_mat <- prepped$letters_mat

  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))

  row_ha <- rowAnnotation(
    Category = ann_df$Category,
    col = list(Category = category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )

  hm <- Heatmap(
    mat,
    name = "Mean CLR",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
    width = grid::unit(8, "cm"),
    row_dend_width = grid::unit(10, "mm"),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 10),
    row_split = factor(ann_df$Category, levels = unname(category_labels)),
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title = title_txt,
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- letters_mat[i, j]
      if (!is.na(lab) && nzchar(lab)) {
        grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 8, col = "black"))
      }
    }
  )

  ht <- row_ha + hm

  lgd_category <- Legend(
    title = "Category",
    at = names(category_colors),
    legend_gp = grid::gpar(fill = unname(category_colors)),
    labels_gp = grid::gpar(fontsize = 10),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold")
  )

  lgd_clr <- Legend(
    title = "Mean CLR",
    col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
    at = c(-2, -1, 0, 1, 2),
    labels_gp = grid::gpar(fontsize = 10),
    title_gp = grid::gpar(fontsize = 11, fontface = "bold")
  )

  pd <- packLegend(
    lgd_category,
    lgd_clr,
    direction = "vertical",
    gap = grid::unit(6, "mm")
  )

  png(file.path(plot_dir, file_name), width = 5200, height = 3600, res = 300)
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE
  )
  draw(
    pd,
    x = grid::unit(0.985, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("right", "top")
  )
  dev.off()
}

plot_upset <- function(class_df, prefix, title_txt) {
  set_cols <- c(
    "Environment_sensitive_5503",
    "Environment_sensitive_6007",
    "Stable_5503_across_soils",
    "Stable_6007_across_soils",
    "Common_in_sandy_between_genotypes",
    "Common_in_clay_between_genotypes",
    "Genotype_different_in_sandy",
    "Genotype_different_in_clay"
  )

  upset_df <- class_df %>%
    select(PathwayLabel, all_of(set_cols)) %>%
    distinct()

  if (nrow(upset_df) == 0) return(invisible(NULL))

  p <- ComplexUpset::upset(
    upset_df,
    intersect = set_cols,
    name = "Pathways",
    width_ratio = 0.18
  ) +
    ggtitle(title_txt)

  ggsave(
    filename = file.path(plot_dir, paste0(prefix, "_upset.png")),
    plot = p,
    width = 14,
    height = 8,
    dpi = 300
  )
}

# -----------------------------
# 3. Leitura
# -----------------------------
pathway_group_means <- read_csv(pathway_group_means_file, show_col_types = FALSE)
kruskal_roots <- read_csv(kruskal_roots_file, show_col_types = FALSE)
kruskal_stalk <- read_csv(kruskal_stalk_file, show_col_types = FALSE)
dunn_roots <- read_csv(dunn_roots_file, show_col_types = FALSE)
dunn_stalk <- read_csv(dunn_stalk_file, show_col_types = FALSE)

# -----------------------------
# 4. Classificacao
# -----------------------------
roots_class <- classify_pathways(dunn_roots, kruskal_roots, "Rhizosphere+Root")
stalk_class <- classify_pathways(dunn_stalk, kruskal_stalk, "Stalk")

write_csv(roots_class, file.path(table_dir, "roots_pathway_classification.csv"))
write_csv(stalk_class, file.path(table_dir, "stalk_pathway_classification.csv"))

# -----------------------------
# 5. UpSet
# -----------------------------
plot_upset(
  roots_class,
  "roots",
  "FunFun KEGG pathways classification intersections - roots"
)

plot_upset(
  stalk_class,
  "stalk",
  "FunFun KEGG pathways classification intersections - stalk"
)

# -----------------------------
# 6. Heatmaps focados
# -----------------------------
roots_rank <- prepare_rank_table(roots_class, pathway_group_means, "root")
stalk_rank <- prepare_rank_table(stalk_class, pathway_group_means, "stalk")

roots_top20 <- select_top_by_category(roots_rank, top_n_table)
stalk_top20 <- select_top_by_category(stalk_rank, top_n_table)
roots_top10 <- select_top_by_category(roots_rank, top_n_heatmap)
stalk_top10 <- select_top_by_category(stalk_rank, top_n_heatmap)

write_csv(roots_top20, file.path(table_dir, "roots_top20_pathways_by_category.csv"))
write_csv(stalk_top20, file.path(table_dir, "stalk_top20_pathways_by_category.csv"))
write_csv(roots_top10, file.path(table_dir, "roots_top10_pathways_heatmap.csv"))
write_csv(stalk_top10, file.path(table_dir, "stalk_top10_pathways_heatmap.csv"))

roots_prep <- build_heatmap_objects(roots_top10, dunn_roots, pathway_group_means, "root")
stalk_prep <- build_heatmap_objects(stalk_top10, dunn_stalk, pathway_group_means, "stalk")

plot_focus_heatmap(
  roots_prep,
  "FunFun KEGG pathways highlighting environmental sensitivity and genotype stability - roots",
  "heatmap_funfun_pathways_roots_letters.png"
)

plot_focus_heatmap(
  stalk_prep,
  "FunFun KEGG pathways highlighting environmental sensitivity and genotype stability - stalk",
  "heatmap_funfun_pathways_stalk_letters.png"
)

# -----------------------------
# 7. Final
# -----------------------------
message("")
message("Classificacao exploratoria dos pathways do FunFun concluida.")
message("Tabelas: ", table_dir)
message("Graficos: ", plot_dir)
