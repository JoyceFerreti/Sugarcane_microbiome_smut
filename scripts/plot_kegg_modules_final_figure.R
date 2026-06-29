# ============================================================
# FIGURA FINAL - KEGG modules
# Usa resultados ja calculados do pipeline anterior:
# - Kruskal-Wallis
# - Dunn post-hoc
# - classificacao dos modulos
#
# Gera heatmaps focados para:
# - Rhizosphere+Root
# - Stalk
#
# Categorias usadas:
# - Stable_5503_across_soils
# - Stable_6007_across_soils
# - Environment_sensitive_5503
# - Environment_sensitive_6007
#
# Criterio de selecao:
# - top 10 por categoria no heatmap
# - top 20 por categoria em tabela
# - ordenar por menor Kruskal_FDR
# - desempatar por maior variancia de CLR
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# -----------------------------
# 0. Diretorio de trabalho
# -----------------------------
project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

# -----------------------------
# 1. Pacotes
# -----------------------------
needed_cran <- c("readr", "dplyr", "tidyr", "tibble", "stringr", "circlize", "multcompView")
needed_bioc <- c("ComplexHeatmap")

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

install_if_missing(needed_cran)
install_if_missing(needed_bioc, bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(circlize)
  library(multcompView)
  library(ComplexHeatmap)
})

# -----------------------------
# 2. Parametros da figura
# -----------------------------
top_n_heatmap <- 10
top_n_table <- 20
group_levels <- c("5503_sandy", "5503_clay", "6007_sandy", "6007_clay")

focus_categories <- c(
  "Stable_5503_across_soils",
  "Stable_6007_across_soils",
  "Environment_sensitive_5503",
  "Environment_sensitive_6007"
)

category_labels <- c(
  Stable_5503_across_soils = "Stable in 5503 across soils",
  Stable_6007_across_soils = "Stable in 6007 across soils",
  Environment_sensitive_5503 = "Environment-sensitive in 5503",
  Environment_sensitive_6007 = "Environment-sensitive in 6007"
)

category_colors <- c(
  "Stable in 5503 across soils" = "#7570b3",
  "Stable in 6007 across soils" = "#e7298a",
  "Environment-sensitive in 5503" = "#d95f02",
  "Environment-sensitive in 6007" = "#e41a1c"
)

# -----------------------------
# 3. Caminhos
# -----------------------------
base_dir <- file.path(project_dir, "picrust2_results", "kegg_modules_analysis")
tables_dir <- file.path(base_dir, "tables")
class_dir <- file.path(base_dir, "classification_upset", "tables")
posthoc_dir <- file.path(base_dir, "posthoc")

out_dir <- file.path(base_dir, "final_figure")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

module_group_means_file <- file.path(tables_dir, "module_group_means_clr.csv")
kruskal_roots_file <- file.path(tables_dir, "kruskal_kegg_modules_roots.csv")
kruskal_stalk_file <- file.path(tables_dir, "kruskal_kegg_modules_stalk.csv")
roots_class_file <- file.path(class_dir, "roots_module_classification.csv")
stalk_class_file <- file.path(class_dir, "stalk_module_classification.csv")
dunn_roots_file <- file.path(posthoc_dir, "dunn_kegg_modules_roots.csv")
dunn_stalk_file <- file.path(posthoc_dir, "dunn_kegg_modules_stalk.csv")

required_files <- c(
  module_group_means_file,
  kruskal_roots_file,
  kruskal_stalk_file,
  roots_class_file,
  stalk_class_file,
  dunn_roots_file,
  dunn_stalk_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    paste0(
      "Os seguintes arquivos nao foram encontrados:\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}

# -----------------------------
# 4. Leitura
# -----------------------------
module_group_means <- read_csv(module_group_means_file, show_col_types = FALSE)
kruskal_roots <- read_csv(kruskal_roots_file, show_col_types = FALSE)
kruskal_stalk <- read_csv(kruskal_stalk_file, show_col_types = FALSE)
roots_class <- read_csv(roots_class_file, show_col_types = FALSE)
stalk_class <- read_csv(stalk_class_file, show_col_types = FALSE)
dunn_roots <- read_csv(dunn_roots_file, show_col_types = FALSE)
dunn_stalk <- read_csv(dunn_stalk_file, show_col_types = FALSE)

# -----------------------------
# 5. Funcoes auxiliares
# -----------------------------
make_module_label <- function(module, module_name) {
  ifelse(
    is.na(module_name) | module_name == "",
    module,
    paste0(module, " | ", module_name)
  )
}

short_module_name <- function(x) {
  sub(" =>.*$", "", x)
}


prepare_rank_table <- function(class_df, kruskal_df, compartment_pattern) {
  means_df <- module_group_means %>%
    mutate(
      SampleTypeClean = str_to_lower(Sample.Type),
      GroupLabel = paste(str_remove(Treatment, "^IAC-"), str_to_lower(Soil), sep = "_")
    ) %>%
    filter(str_detect(SampleTypeClean, compartment_pattern))
  
  variability_df <- means_df %>%
    group_by(Module, Module_Name) %>%
    summarise(CLR_Variance = var(Mean_CLR, na.rm = TRUE), .groups = "drop")
  
  class_long <- class_df %>%
    select(Module, Module_Name, all_of(focus_categories)) %>%
    pivot_longer(
      cols = all_of(focus_categories),
      names_to = "CategoryKey",
      values_to = "InCategory"
    ) %>%
    filter(InCategory) %>%
    mutate(Category = recode(CategoryKey, !!!category_labels))
  
  kruskal_rank <- kruskal_df %>%
    transmute(
      Module,
      Module_Name,
      Kruskal_p,
      Kruskal_FDR = ifelse(is.na(Kruskal_FDR), 1, Kruskal_FDR),
      Kruskal_rank_value = ifelse(is.na(Kruskal_FDR), Kruskal_p, Kruskal_FDR)
    )
  
  class_long %>%
    left_join(kruskal_rank, by = c("Module", "Module_Name")) %>%
    left_join(variability_df, by = c("Module", "Module_Name")) %>%
    mutate(
      CLR_Variance = ifelse(is.na(CLR_Variance), 0, CLR_Variance),
      Kruskal_rank_value = ifelse(is.na(Kruskal_rank_value), 1, Kruskal_rank_value)
    ) %>%
    arrange(Category, Kruskal_rank_value, desc(CLR_Variance))
}

select_top_by_category <- function(rank_df, n) {
  rank_df %>%
    group_by(Category) %>%
    arrange(Kruskal_rank_value, desc(CLR_Variance), .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup()
}

make_letters_for_module <- function(df_mod, group_levels) {
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
    pval <- df_mod$p.adj[i]
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
  modules <- unique(selected_df$Module)
  bind_rows(lapply(modules, function(mod) {
    df_mod <- dunn_df %>% filter(Module == mod)
    letters_df <- make_letters_for_module(df_mod, group_levels)
    tibble(Module = mod, GroupLabel = letters_df$GroupLabel, Letters = letters_df$Letters)
  }))
}

build_heatmap_objects <- function(selected_df, dunn_df, compartment_pattern) {
  means_df <- module_group_means %>%
    mutate(
      SampleTypeClean = str_to_lower(Sample.Type),
      GroupLabel = paste(str_remove(Treatment, "^IAC-"), str_to_lower(Soil), sep = "_")
    ) %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    filter(Module %in% selected_df$Module) %>%
    mutate(
      GroupLabel = factor(GroupLabel, levels = group_levels),
      ModuleLabel = make_module_label(Module, Module_Name),
      ModuleLabel = short_module_name(ModuleLabel)
    ) %>%
    select(Module, Module_Name, ModuleLabel, GroupLabel, Mean_CLR) %>%
    distinct()
  
  mat <- means_df %>%
    select(ModuleLabel, GroupLabel, Mean_CLR) %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame()
  
  rownames(mat) <- mat$ModuleLabel
  mat$ModuleLabel <- NULL
  mat <- as.matrix(mat)
  
  ann_df <- selected_df %>%
    mutate(
      ModuleLabel = make_module_label(Module, Module_Name),
      ModuleLabel = short_module_name(ModuleLabel)
    ) %>%
    distinct(Module, Module_Name, Category, ModuleLabel, Kruskal_p, Kruskal_FDR, CLR_Variance) %>%
    filter(ModuleLabel %in% rownames(mat))
  
  ann_df <- ann_df[match(rownames(mat), ann_df$ModuleLabel), , drop = FALSE]
  
  letters_long <- build_letters_table(selected_df, dunn_df) %>%
    left_join(
      ann_df %>% select(Module, ModuleLabel),
      by = "Module"
    )
  
  letters_mat <- letters_long %>%
    select(ModuleLabel, GroupLabel, Letters) %>%
    pivot_wider(names_from = GroupLabel, values_from = Letters, values_fill = "") %>%
    as.data.frame()
  
  rownames(letters_mat) <- letters_mat$ModuleLabel
  letters_mat$ModuleLabel <- NULL
  letters_mat <- as.matrix(letters_mat)
  letters_mat <- letters_mat[rownames(mat), colnames(mat), drop = FALSE]
  
  list(mat = mat, ann_df = ann_df, letters_mat = letters_mat)
}

plot_focus_heatmap <- function(prepped, title_txt, png_file, pdf_file) {
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
    width = grid::unit(9, "cm"),
    row_dend_width = grid::unit(12, "mm"),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    row_split = factor(ann_df$Category, levels = unique(ann_df$Category)),
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title = title_txt,
    show_heatmap_legend = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- letters_mat[i, j]
      if (!is.na(lab) && nzchar(lab)) {
        grid::grid.text(lab, x, y, gp = grid::gpar(fontsize = 9, col = "black"))
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
  
  png(png_file, width = 5200, height = 3600, res = 300)
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
  
  pdf(pdf_file, width = 17, height = 12)
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

# -----------------------------
# 6. Selecionar modulos
# -----------------------------
roots_rank <- prepare_rank_table(roots_class, kruskal_roots, "root")
stalk_rank <- prepare_rank_table(stalk_class, kruskal_stalk, "stalk")

roots_top20 <- select_top_by_category(roots_rank, top_n_table)
stalk_top20 <- select_top_by_category(stalk_rank, top_n_table)
roots_top10 <- select_top_by_category(roots_rank, top_n_heatmap)
stalk_top10 <- select_top_by_category(stalk_rank, top_n_heatmap)

write_csv(roots_top20, file.path(table_dir, "roots_top20_modules_by_focus_category.csv"))
write_csv(stalk_top20, file.path(table_dir, "stalk_top20_modules_by_focus_category.csv"))
write_csv(roots_top10, file.path(table_dir, "roots_top10_modules_used_in_heatmap.csv"))
write_csv(stalk_top10, file.path(table_dir, "stalk_top10_modules_used_in_heatmap.csv"))

roots_letters <- build_letters_table(roots_top10, dunn_roots)
stalk_letters <- build_letters_table(stalk_top10, dunn_stalk)

write_csv(roots_letters, file.path(table_dir, "roots_dunn_letters_for_heatmap.csv"))
write_csv(stalk_letters, file.path(table_dir, "stalk_dunn_letters_for_heatmap.csv"))

# -----------------------------
# 7. Plotar heatmaps
# -----------------------------
roots_prep <- build_heatmap_objects(roots_top10, dunn_roots, "root")
stalk_prep <- build_heatmap_objects(stalk_top10, dunn_stalk, "stalk")

plot_focus_heatmap(
  roots_prep,
  "KEGG modules highlighting genotype stability and environmental sensitivity - roots",
  file.path(plot_dir, "heatmap_kegg_modules_focus_roots_letters.png"),
  file.path(plot_dir, "heatmap_kegg_modules_focus_roots_letters.pdf")
)

plot_focus_heatmap(
  stalk_prep,
  "KEGG modules highlighting genotype stability and environmental sensitivity - stalk",
  file.path(plot_dir, "heatmap_kegg_modules_focus_stalk_letters.png"),
  file.path(plot_dir, "heatmap_kegg_modules_focus_stalk_letters.pdf")
)

message("")
message("Figura final gerada com sucesso.")
message("Diretorio de trabalho: ", project_dir)
message("Saidas em: ", out_dir)
message("")
message("Arquivos principais:")
message("- ", file.path(plot_dir, "heatmap_kegg_modules_focus_roots_letters.png"))
message("- ", file.path(plot_dir, "heatmap_kegg_modules_focus_stalk_letters.png"))
message("- ", file.path(plot_dir, "heatmap_kegg_modules_focus_roots_letters.pdf"))
message("- ", file.path(plot_dir, "heatmap_kegg_modules_focus_stalk_letters.pdf"))
