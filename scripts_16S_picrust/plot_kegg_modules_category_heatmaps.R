# ============================================================
# Heatmaps de modulos KEGG com categorias biologicas
# Reutiliza os resultados ja gerados em:
# picrust2_results/kegg_modules_analysis/
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

base_dir <- file.path(project_dir, "picrust2_results", "kegg_modules_analysis")
tables_dir <- file.path(base_dir, "tables")
class_dir <- file.path(base_dir, "classification_upset", "tables")
plot_dir <- file.path(base_dir, "category_heatmaps")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

module_group_means_file <- file.path(tables_dir, "module_group_means_clr.csv")
roots_class_file <- file.path(class_dir, "roots_module_classification.csv")
stalk_class_file <- file.path(class_dir, "stalk_module_classification.csv")

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

stop_if_missing(module_group_means_file)
stop_if_missing(roots_class_file)
stop_if_missing(stalk_class_file)

module_group_means <- read_csv(module_group_means_file, show_col_types = FALSE)
roots_class <- read_csv(roots_class_file, show_col_types = FALSE)
stalk_class <- read_csv(stalk_class_file, show_col_types = FALSE)

category_priority <- c(
  "Stable_both_genotypes_across_soils",
  "Common_both_soils_between_genotypes",
  "Stable_5503_across_soils",
  "Stable_6007_across_soils",
  "Common_in_sandy_between_genotypes",
  "Common_in_clay_between_genotypes",
  "Environment_sensitive_5503",
  "Environment_sensitive_6007",
  "Genotype_different_in_sandy",
  "Genotype_different_in_clay"
)

category_labels <- c(
  Stable_both_genotypes_across_soils = "Stable in both genotypes across soils",
  Common_both_soils_between_genotypes = "Common between genotypes in both soils",
  Stable_5503_across_soils = "Stable in 5503 across soils",
  Stable_6007_across_soils = "Stable in 6007 across soils",
  Common_in_sandy_between_genotypes = "Common between genotypes in sandy soil",
  Common_in_clay_between_genotypes = "Common between genotypes in clay soil",
  Environment_sensitive_5503 = "Environment-sensitive in 5503",
  Environment_sensitive_6007 = "Environment-sensitive in 6007",
  Genotype_different_in_sandy = "Genotype-different in sandy soil",
  Genotype_different_in_clay = "Genotype-different in clay soil"
)

category_colors <- c(
  "Stable in both genotypes across soils" = "#1b9e77",
  "Common between genotypes in both soils" = "#66a61e",
  "Stable in 5503 across soils" = "#7570b3",
  "Stable in 6007 across soils" = "#e7298a",
  "Common between genotypes in sandy soil" = "#e6ab02",
  "Common between genotypes in clay soil" = "#a6761d",
  "Environment-sensitive in 5503" = "#d95f02",
  "Environment-sensitive in 6007" = "#e41a1c",
  "Genotype-different in sandy soil" = "#377eb8",
  "Genotype-different in clay soil" = "#4daf4a"
)

assign_principal_category <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      PrincipalCategoryKey = {
        hits <- category_priority[unlist(c_across(all_of(category_priority)))]
        if (length(hits) == 0) NA_character_ else hits[1]
      },
      PrincipalCategory = dplyr::recode(PrincipalCategoryKey, !!!category_labels)
    ) %>%
    ungroup()
}

make_module_label <- function(module, module_name) {
  ifelse(
    is.na(module_name) | module_name == "",
    module,
    paste0(module, " | ", module_name)
  )
}

prepare_heatmap_data <- function(group_means_df, class_df, compartment_pattern) {
  class_df <- assign_principal_category(class_df) %>%
    filter(!is.na(PrincipalCategory)) %>%
    mutate(ModuleLabel = make_module_label(Module, Module_Name))

  means_df <- group_means_df %>%
    mutate(
      SampleTypeClean = str_to_lower(Sample.Type),
      GroupLabel = paste(str_remove(Treatment, "^IAC-"), str_to_lower(Soil), sep = "_"),
      ModuleLabel = make_module_label(Module, Module_Name)
    ) %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    mutate(
      GroupLabel = factor(GroupLabel, levels = c("5503_sandy", "5503_clay", "6007_sandy", "6007_clay"))
    ) %>%
    filter(Module %in% class_df$Module) %>%
    select(Module, ModuleLabel, GroupLabel, Mean_CLR) %>%
    distinct()

  mat <- means_df %>%
    select(ModuleLabel, GroupLabel, Mean_CLR) %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame()

  rownames(mat) <- mat$ModuleLabel
  mat$ModuleLabel <- NULL
  mat <- as.matrix(mat)

  ann_df <- class_df %>%
    select(Module, Module_Name, PrincipalCategory, all_of(category_priority)) %>%
    mutate(ModuleLabel = make_module_label(Module, Module_Name)) %>%
    distinct() %>%
    filter(ModuleLabel %in% rownames(mat))

  ann_df <- ann_df[match(rownames(mat), ann_df$ModuleLabel), , drop = FALSE]

  list(mat = mat, ann_df = ann_df)
}

plot_category_heatmap <- function(prepped, file_name, title_txt) {
  mat <- prepped$mat
  ann_df <- prepped$ann_df

  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    return(invisible(NULL))
  }

  row_ha <- rowAnnotation(
    Category = ann_df$PrincipalCategory,
    col = list(Category = category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )

  hm <- Heatmap(
    mat,
    name = "CLR",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 10),
    row_split = factor(ann_df$PrincipalCategory, levels = unique(ann_df$PrincipalCategory)),
    column_title = title_txt,
    heatmap_legend_param = list(title = "Mean CLR")
  )

  png(file.path(plot_dir, file_name), width = 2400, height = 3200, res = 300)
  draw(
    row_ha + hm,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = FALSE
  )
  dev.off()
}

roots_prep <- prepare_heatmap_data(module_group_means, roots_class, "root")
stalk_prep <- prepare_heatmap_data(module_group_means, stalk_class, "stalk")

plot_category_heatmap(
  roots_prep,
  "heatmap_kegg_modules_roots_by_category.png",
  "KEGG modules by biological category - roots"
)

plot_category_heatmap(
  stalk_prep,
  "heatmap_kegg_modules_stalk_by_category.png",
  "KEGG modules by biological category - stalk"
)

roots_export <- assign_principal_category(roots_class)
stalk_export <- assign_principal_category(stalk_class)

write_csv(roots_export, file.path(plot_dir, "roots_module_classification_with_principal_category.csv"))
write_csv(stalk_export, file.path(plot_dir, "stalk_module_classification_with_principal_category.csv"))

message("")
message("Heatmaps por categoria biologica concluidos.")
message("Saida: ", plot_dir)
message("- roots: ", file.path(plot_dir, "heatmap_kegg_modules_roots_by_category.png"))
message("- stalk: ", file.path(plot_dir, "heatmap_kegg_modules_stalk_by_category.png"))
