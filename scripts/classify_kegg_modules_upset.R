# ============================================================
# Classificacao biologica de modulos KEGG a partir do Dunn post-hoc
# Gera:
# - listas de modulos por padrao de contraste
# - tabelas com IDs e nomes dos modulos
# - graficos UpSet para roots e stalk
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

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr", "ComplexUpset", "ggplot2"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ComplexUpset)
  library(ggplot2)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

base_dir <- file.path(project_dir, "picrust2_results", "kegg_modules_analysis")
tables_dir <- file.path(base_dir, "tables")
posthoc_dir <- file.path(base_dir, "posthoc")
out_dir <- file.path(base_dir, "classification_upset")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

dunn_roots_file <- file.path(posthoc_dir, "dunn_kegg_modules_roots.csv")
dunn_stalk_file <- file.path(posthoc_dir, "dunn_kegg_modules_stalk.csv")
kruskal_roots_file <- file.path(tables_dir, "kruskal_kegg_modules_roots.csv")
kruskal_stalk_file <- file.path(tables_dir, "kruskal_kegg_modules_stalk.csv")

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

stop_if_missing(dunn_roots_file)
stop_if_missing(dunn_stalk_file)
stop_if_missing(kruskal_roots_file)
stop_if_missing(kruskal_stalk_file)

dunn_roots <- read_csv(dunn_roots_file, show_col_types = FALSE)
dunn_stalk <- read_csv(dunn_stalk_file, show_col_types = FALSE)
kruskal_roots <- read_csv(kruskal_roots_file, show_col_types = FALSE)
kruskal_stalk <- read_csv(kruskal_stalk_file, show_col_types = FALSE)

groups_expected <- c("5503_sandy", "5503_clay", "6007_sandy", "6007_clay")

normalize_pair_key <- function(g1, g2) {
  pair <- sort(c(g1, g2))
  paste(pair, collapse = "__")
}

get_pair_status <- function(dunn_df, alpha = 0.05) {
  if (nrow(dunn_df) == 0) {
    return(tibble())
  }

  dunn_df %>%
    mutate(
      PairKey = mapply(normalize_pair_key, group1, group2),
      Significant = !is.na(p.adj) & p.adj <= alpha
    ) %>%
    select(Compartimento, Module, Module_Name, PairKey, Significant) %>%
    distinct()
}

classify_modules <- function(dunn_df, kruskal_df, compartment_name) {
  pair_df <- get_pair_status(dunn_df)

  sig_mods <- kruskal_df %>%
    filter((!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) | (!is.na(Kruskal_p) & Kruskal_p <= 0.05)) %>%
    select(Module, Module_Name) %>%
    distinct()

  if (nrow(sig_mods) == 0) {
    return(tibble())
  }

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

  expanded <- sig_mods %>%
    tidyr::crossing(expected_pairs) %>%
    left_join(pair_df, by = c("Module", "Module_Name", "PairKey")) %>%
    mutate(
      Compartimento = compartment_name,
      Significant = dplyr::coalesce(Significant, FALSE),
      NoDifference = !Significant
    )

  summary_tbl <- expanded %>%
    select(Compartimento, Module, Module_Name, ContrastType, Significant, NoDifference) %>%
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
      Genotype_different_in_clay = Significant_same_soil_clay_across_genotypes,
      Stable_both_genotypes_across_soils =
        NoDifference_same_genotype_5503_across_soils &
        NoDifference_same_genotype_6007_across_soils,
      Common_both_soils_between_genotypes =
        NoDifference_same_soil_sandy_across_genotypes &
        NoDifference_same_soil_clay_across_genotypes
    )

  summary_tbl
}

write_category_tables <- function(class_df, prefix) {
  category_cols <- c(
    "Stable_5503_across_soils",
    "Stable_6007_across_soils",
    "Environment_sensitive_5503",
    "Environment_sensitive_6007",
    "Common_in_sandy_between_genotypes",
    "Common_in_clay_between_genotypes",
    "Genotype_different_in_sandy",
    "Genotype_different_in_clay",
    "Stable_both_genotypes_across_soils",
    "Common_both_soils_between_genotypes"
  )

  for (col in category_cols) {
    out <- class_df %>%
      filter(.data[[col]]) %>%
      select(Compartimento, Module, Module_Name, all_of(category_cols))

    write_csv(out, file.path(table_dir, paste0(prefix, "_", col, ".csv")))
  }
}

plot_upset <- function(class_df, prefix, title_txt) {
  set_cols <- c(
    "Stable_5503_across_soils",
    "Stable_6007_across_soils",
    "Environment_sensitive_5503",
    "Environment_sensitive_6007",
    "Common_in_sandy_between_genotypes",
    "Common_in_clay_between_genotypes",
    "Genotype_different_in_sandy",
    "Genotype_different_in_clay"
  )

  upset_df <- class_df %>%
    select(Module, Module_Name, all_of(set_cols))

  if (nrow(upset_df) == 0) {
    return(invisible(NULL))
  }

  p <- ComplexUpset::upset(
    upset_df,
    intersect = set_cols,
    name = "Modules",
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

roots_class <- classify_modules(dunn_roots, kruskal_roots, "Rhizosphere+Root")
stalk_class <- classify_modules(dunn_stalk, kruskal_stalk, "Stalk")

write_csv(roots_class, file.path(table_dir, "roots_module_classification.csv"))
write_csv(stalk_class, file.path(table_dir, "stalk_module_classification.csv"))

write_category_tables(roots_class, "roots")
write_category_tables(stalk_class, "stalk")

plot_upset(
  roots_class,
  "roots",
  "KEGG modules classification intersections - roots"
)

plot_upset(
  stalk_class,
  "stalk",
  "KEGG modules classification intersections - stalk"
)

message("")
message("Classificacao de modulos concluida.")
message("Tabelas: ", table_dir)
message("Graficos UpSet: ", plot_dir)
message("")
message("Arquivos principais:")
message("- roots: ", file.path(table_dir, "roots_module_classification.csv"))
message("- stalk: ", file.path(table_dir, "stalk_module_classification.csv"))
message("- upset roots: ", file.path(plot_dir, "roots_upset.png"))
message("- upset stalk: ", file.path(plot_dir, "stalk_upset.png"))
