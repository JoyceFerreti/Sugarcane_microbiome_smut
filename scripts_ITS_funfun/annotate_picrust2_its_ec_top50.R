# ============================================================
# Anotacao exploratoria de EC para ITS
# Seleciona os top 50 EC por compartimento (roots e stalk)
# com base na variancia de CLR entre grupos e busca
# anotacoes funcionais para apoiar os heatmaps exploratorios.
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

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr"))
install_if_missing(c("KEGGREST"), bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(KEGGREST)
})

# -----------------------------
# 1. Configuracao
# -----------------------------
project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

ec_group_means_file <- file.path(project_dir, "picrust2_results", "EC_analysis", "tables", "ec_group_means_clr.csv")

analysis_dir <- file.path(project_dir, "picrust2_results", "EC_analysis", "annotation")
tables_dir <- file.path(analysis_dir, "tables")
cache_dir <- file.path(analysis_dir, "cache")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

top_n <- 50

# -----------------------------
# 2. Funcoes auxiliares
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

annotate_ec_vector <- function(ec_ids, cache_file) {
  if (file.exists(cache_file)) {
    cached <- read_csv(cache_file, show_col_types = FALSE)
  } else {
    cached <- tibble(EC = character(), KEGG_Name = character(), KEGG_Definition = character())
  }

  missing_ec <- setdiff(unique(ec_ids), cached$EC)

  if (length(missing_ec) > 0) {
    message("Buscando anotacoes de EC no KEGG...")

    new_ann <- lapply(missing_ec, function(ec) {
      res <- tryCatch(KEGGREST::keggGet(paste0("ec:", ec)), error = function(e) NULL)

      if (is.null(res) || length(res) == 0) {
        return(tibble(EC = ec, KEGG_Name = NA_character_, KEGG_Definition = NA_character_))
      }

      entry <- res[[1]]
      name_val <- if (!is.null(entry$NAME)) paste(entry$NAME, collapse = "; ") else NA_character_
      def_val <- if (!is.null(entry$DEFINITION)) entry$DEFINITION else NA_character_

      tibble(EC = ec, KEGG_Name = name_val, KEGG_Definition = def_val)
    }) %>%
      bind_rows()

    cached <- bind_rows(cached, new_ann) %>%
      distinct(EC, .keep_all = TRUE)

    write_csv(cached, cache_file)
  }

  cached
}

select_top_ec <- function(df, compartment_pattern, n = 50) {
  df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    group_by(EC) %>%
    summarise(CLR_Variance = var(Mean_CLR, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(CLR_Variance)) %>%
    slice_head(n = n)
}

# -----------------------------
# 3. Leitura
# -----------------------------
stop_if_missing(ec_group_means_file)

ec_group_means <- read_csv(ec_group_means_file, show_col_types = FALSE)

# -----------------------------
# 4. Selecionar top 50 por compartimento
# -----------------------------
top50_roots <- select_top_ec(ec_group_means, "root", n = top_n) %>%
  mutate(Compartimento = "Rhizosphere+Root")

top50_stalk <- select_top_ec(ec_group_means, "stalk", n = top_n) %>%
  mutate(Compartimento = "Stalk")

top100_all <- bind_rows(top50_roots, top50_stalk)

write_csv(top50_roots, file.path(tables_dir, "top50_EC_roots_by_variance.csv"))
write_csv(top50_stalk, file.path(tables_dir, "top50_EC_stalk_by_variance.csv"))

# -----------------------------
# 5. Anotacao
# -----------------------------
cache_file <- file.path(cache_dir, "EC_annotation_cache.csv")
ec_annotation <- annotate_ec_vector(top100_all$EC, cache_file)

top50_roots_annot <- top50_roots %>%
  left_join(ec_annotation, by = "EC")

top50_stalk_annot <- top50_stalk %>%
  left_join(ec_annotation, by = "EC")

write_csv(top50_roots_annot, file.path(tables_dir, "top50_EC_roots_annotated.csv"))
write_csv(top50_stalk_annot, file.path(tables_dir, "top50_EC_stalk_annotated.csv"))

# -----------------------------
# 6. Tabela expandida para heatmap exploratorio
# -----------------------------
roots_heatmap_table <- ec_group_means %>%
  filter(EC %in% top50_roots$EC, str_detect(SampleTypeClean, "root")) %>%
  left_join(ec_annotation, by = "EC")

stalk_heatmap_table <- ec_group_means %>%
  filter(EC %in% top50_stalk$EC, str_detect(SampleTypeClean, "stalk")) %>%
  left_join(ec_annotation, by = "EC")

write_csv(roots_heatmap_table, file.path(tables_dir, "roots_heatmap_EC_top50_annotated_long.csv"))
write_csv(stalk_heatmap_table, file.path(tables_dir, "stalk_heatmap_EC_top50_annotated_long.csv"))

# -----------------------------
# 7. Final
# -----------------------------
message("")
message("Anotacao exploratoria de EC concluida.")
message("Tabelas: ", tables_dir)
message("Cache: ", cache_file)
message("")
message("Arquivos principais:")
message("- ", file.path(tables_dir, "top50_EC_roots_annotated.csv"))
message("- ", file.path(tables_dir, "top50_EC_stalk_annotated.csv"))
