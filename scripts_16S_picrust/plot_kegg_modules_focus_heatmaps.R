# ============================================================
# Heatmaps focados de modulos KEGG
# Categorias:
# - Stable_5503_across_soils
# - Stable_6007_across_soils
# - Environment_sensitive_5503
# - Environment_sensitive_6007
# Criterio:
# - ordenar por menor Kruskal_FDR
# - desempatar por maior variancia de CLR entre grupos
# Saidas:
# - top 10 por categoria para heatmap
# - top 20 por categoria para tabelas
# ============================================================

options(stringsAsFactors = FALSE)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) == 0) {
    return(invisible(NULL))
  }

  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE)
  } else {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr", "circlize"))
install_if_missing(c("ComplexHeatmap"), bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(circlize)
  library(ComplexHeatmap)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

top_n_heatmap <- 10
top_n_table <- 20

base_dir <- file.path(project_dir, "picrust2_results", "kegg_modules_analysis")
tables_dir <- file.path(base_dir, "tables")
class_dir <- file.path(base_dir, "classification_upset", "tables")
out_dir <- file.path(base_dir, "focused_category_heatmaps")
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

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

for (f in c(
  module_group_means_file,
  kruskal_roots_file,
  kruskal_stalk_file,
  roots_class_file,
  stalk_class_file
)) {
  stop_if_missing(f)
}

module_group_means <- read_csv(module_group_means_file, show_col_types = FALSE)
kruskal_roots <- read_csv(kruskal_roots_file, show_col_types = FALSE)
kruskal_stalk <- read_csv(kruskal_stalk_file, show_col_types = FALSE)
roots_class <- read_csv(roots_class_file, show_col_types = FALSE)
stalk_class <- read_csv(stalk_class_file, show_col_types = FALSE)

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

make_module_label <- function(module, module_name) {
  ifelse(
    is.na(module_name) | module_name == "",
    module,
    paste0(module, " | ", module_name)
  )
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

build_heatmap_matrix <- function(selected_df, compartment_pattern) {
  means_df <- module_group_means %>%
    mutate(
      SampleTypeClean = str_to_lower(Sample.Type),
      GroupLabel = paste(str_remove(Treatment, "^IAC-"), str_to_lower(Soil), sep = "_")
    ) %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    filter(Module %in% selected_df$Module) %>%
    mutate(
      GroupLabel = factor(GroupLabel, levels = c("5503_sandy", "5503_clay", "6007_sandy", "6007_clay")),
      ModuleLabel = make_module_label(Module, Module_Name)
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
    mutate(ModuleLabel = make_module_label(Module, Module_Name)) %>%
    distinct(Module, Module_Name, Category, ModuleLabel, Kruskal_p, Kruskal_FDR, CLR_Variance) %>%
    filter(ModuleLabel %in% rownames(mat))

  ann_df <- ann_df[match(rownames(mat), ann_df$ModuleLabel), , drop = FALSE]

  list(mat = mat, ann_df = ann_df)
}

plot_focus_heatmap <- function(prepped, title_txt, file_name) {
  mat <- prepped$mat
  ann_df <- prepped$ann_df

  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    return(invisible(NULL))
  }

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
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 10),
    row_split = factor(ann_df$Category, levels = unique(ann_df$Category)),
    row_title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title = title_txt,
    heatmap_legend_param = list(title = "Mean CLR")
  )

  png(file.path(plot_dir, file_name), width = 2200, height = 2600, res = 300)
  draw(
    row_ha + hm,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  dev.off()
}

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

roots_prep <- build_heatmap_matrix(roots_top10, "root")
stalk_prep <- build_heatmap_matrix(stalk_top10, "stalk")

plot_focus_heatmap(
  roots_prep,
  "KEGG modules highlighting genotype stability and environmental sensitivity - roots",
  "heatmap_kegg_modules_focus_roots.png"
)

plot_focus_heatmap(
  stalk_prep,
  "KEGG modules highlighting genotype stability and environmental sensitivity - stalk",
  "heatmap_kegg_modules_focus_stalk.png"
)

message("")
message("Heatmaps focados concluidos.")
message("Tabelas: ", table_dir)
message("Graficos: ", plot_dir)
message("")
message("Arquivos principais:")
message("- roots heatmap: ", file.path(plot_dir, "heatmap_kegg_modules_focus_roots.png"))
message("- stalk heatmap: ", file.path(plot_dir, "heatmap_kegg_modules_focus_stalk.png"))
message("- roots top20: ", file.path(table_dir, "roots_top20_modules_by_focus_category.csv"))
message("- stalk top20: ", file.path(table_dir, "stalk_top20_modules_by_focus_category.csv"))
