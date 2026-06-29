# ============================================================
# Analise de KO (KEGG Orthologs) a partir do PICRUSt2
# Projeto: 16S / PICRUSt2
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

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr", "pheatmap"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
})

# Ative isto apenas se quiser tentar buscar nomes dos KO pela API do KEGG.
fetch_kegg_names <- TRUE

if (fetch_kegg_names) {
  install_if_missing(c("KEGGREST"), bioc = TRUE)
  suppressPackageStartupMessages(library(KEGGREST))
}

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

metadata_file <- file.path(project_dir, "mapp_plant.csv")
ko_file <- file.path(project_dir, "picrust2_results", "summary_tables", "KO_predictions.csv")

analysis_dir <- file.path(project_dir, "picrust2_results", "kegg_analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

plots_dir <- file.path(analysis_dir, "plots")
tables_dir <- file.path(analysis_dir, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

stop_if_missing(metadata_file)
stop_if_missing(ko_file)

metadata <- read_delim(metadata_file, delim = ";", show_col_types = FALSE, trim_ws = TRUE)
ko_tbl <- read_csv(ko_file, show_col_types = FALSE)

ko_id_col <- names(ko_tbl)[1]
names(ko_tbl)[1] <- "KO"

ko_tbl <- ko_tbl %>%
  mutate(
    KO = str_trim(KO),
    KO = str_remove(KO, "^ko:"),
    KO = str_extract(KO, "K\\d{5}")
  ) %>%
  filter(!is.na(KO))

sample_cols <- intersect(names(ko_tbl)[-1], metadata$Sample)

if (length(sample_cols) == 0) {
  stop("Nenhuma coluna de amostra da tabela de KO foi encontrada nos metadados.", call. = FALSE)
}

ko_long <- ko_tbl %>%
  select(KO, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(metadata, by = "Sample")

ko_by_treatment <- ko_long %>%
  group_by(KO, Treatment) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

ko_by_treatment_sampletype <- ko_long %>%
  mutate(Group = paste(Treatment, Sample.Type, sep = "__")) %>%
  group_by(KO, Group) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

ko_by_all_groups <- ko_long %>%
  mutate(Group = paste(Soil, Treatment, Sample.Type, sep = "__")) %>%
  group_by(KO, Group) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

write_csv(ko_by_treatment, file.path(tables_dir, "KO_mean_by_treatment.csv"))
write_csv(ko_by_treatment_sampletype, file.path(tables_dir, "KO_mean_by_treatment_sampletype.csv"))
write_csv(ko_by_all_groups, file.path(tables_dir, "KO_mean_by_soil_treatment_sampletype.csv"))

select_top_kos <- function(df, n = 50) {
  df %>%
    group_by(KO) %>%
    summarise(GlobalMean = mean(MeanAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(GlobalMean)) %>%
    slice_head(n = n) %>%
    pull(KO)
}

make_matrix <- function(df, group_col, top_kos) {
  df %>%
    filter(KO %in% top_kos) %>%
    select(KO, !!rlang::sym(group_col), MeanAbundance) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(group_col), values_from = MeanAbundance, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("KO") %>%
    as.matrix()
}

top50_treatment <- select_top_kos(ko_by_treatment, n = 50)
top50_treatment_sampletype <- select_top_kos(ko_by_treatment_sampletype, n = 50)
top50_all_groups <- select_top_kos(ko_by_all_groups, n = 50)

mat_treatment <- make_matrix(ko_by_treatment, "Treatment", top50_treatment)
mat_treatment_sampletype <- make_matrix(ko_by_treatment_sampletype, "Group", top50_treatment_sampletype)
mat_all_groups <- make_matrix(ko_by_all_groups, "Group", top50_all_groups)

write_csv(
  tibble(KO = rownames(mat_treatment)),
  file.path(tables_dir, "top50_KO_by_treatment_ids.csv")
)

write_csv(
  tibble(KO = rownames(mat_treatment_sampletype)),
  file.path(tables_dir, "top50_KO_by_treatment_sampletype_ids.csv")
)

write_csv(
  tibble(KO = rownames(mat_all_groups)),
  file.path(tables_dir, "top50_KO_by_all_groups_ids.csv")
)

pheatmap(
  mat_treatment,
  scale = "row",
  clustering_method = "complete",
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Top 50 KO por Treatment",
  filename = file.path(plots_dir, "heatmap_KO_top50_by_treatment.png"),
  width = 9,
  height = 12
)

pheatmap(
  mat_treatment_sampletype,
  scale = "row",
  clustering_method = "complete",
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Top 50 KO por Treatment e Sample.Type",
  filename = file.path(plots_dir, "heatmap_KO_top50_by_treatment_sampletype.png"),
  width = 12,
  height = 14
)

pheatmap(
  mat_all_groups,
  scale = "row",
  clustering_method = "complete",
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 50 KO por Soil, Treatment e Sample.Type",
  filename = file.path(plots_dir, "heatmap_KO_top50_by_all_groups.png"),
  width = 14,
  height = 16
)

annotate_kos <- function(kos) {
  kos <- unique(na.omit(str_extract(kos, "K\\d{5}")))

  if (!fetch_kegg_names) {
    return(tibble(KO = kos, KEGG_Name = NA_character_, KEGG_Definition = NA_character_))
  }

  out <- lapply(kos, function(k) {
    res <- tryCatch(KEGGREST::keggGet(paste0("ko:", k)), error = function(e) NULL)
    if (is.null(res) || length(res) == 0) {
      return(tibble(KO = k, KEGG_Name = NA_character_, KEGG_Definition = NA_character_))
    }

    entry <- res[[1]]
    name_val <- if (!is.null(entry$NAME)) paste(entry$NAME, collapse = "; ") else NA_character_
    def_val <- if (!is.null(entry$DEFINITION)) entry$DEFINITION else NA_character_

    tibble(KO = k, KEGG_Name = name_val, KEGG_Definition = def_val)
  })

  bind_rows(out)
}

make_ko_label_map <- function(annotation_df) {
  annotation_df %>%
    mutate(
      KEGG_Name = dplyr::coalesce(KEGG_Name, "Unknown function"),
      KO_Label = ifelse(
        is.na(KEGG_Name) | KEGG_Name == "Unknown function",
        KO,
        paste0(KO, " | ", KEGG_Name)
      )
    ) %>%
    select(KO, KO_Label, KEGG_Name, KEGG_Definition) %>%
    distinct()
}

apply_heatmap_labels <- function(mat, label_map) {
  label_vec <- label_map$KO_Label
  names(label_vec) <- label_map$KO

  rownames(mat) <- ifelse(
    rownames(mat) %in% names(label_vec),
    unname(label_vec[rownames(mat)]),
    rownames(mat)
  )

  mat
}

all_kos <- unique(ko_long$KO)
ko_annotation_all <- annotate_kos(all_kos)
ko_label_map <- make_ko_label_map(ko_annotation_all)

ko_long_annotated <- ko_long %>%
  left_join(ko_label_map, by = "KO")

write_csv(ko_long, file.path(tables_dir, "KO_long_with_metadata.csv"))
write_csv(ko_long_annotated, file.path(tables_dir, "KO_long_with_metadata_annotated.csv"))
write_csv(ko_annotation_all, file.path(tables_dir, "KO_annotation_all_ids.csv"))

all_top_kos <- unique(c(rownames(mat_treatment), rownames(mat_treatment_sampletype), rownames(mat_all_groups)))
ko_annotation_top <- ko_label_map %>%
  filter(KO %in% all_top_kos)

write_csv(ko_annotation_top, file.path(tables_dir, "KO_annotation_top_ids.csv"))

mat_treatment_labeled <- apply_heatmap_labels(mat_treatment, ko_label_map)
mat_treatment_sampletype_labeled <- apply_heatmap_labels(mat_treatment_sampletype, ko_label_map)
mat_all_groups_labeled <- apply_heatmap_labels(mat_all_groups, ko_label_map)

pheatmap(
  mat_treatment_labeled,
  scale = "row",
  clustering_method = "complete",
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Top 50 KO por Treatment",
  filename = file.path(plots_dir, "heatmap_KO_top50_by_treatment_annotated.png"),
  width = 12,
  height = 14
)

pheatmap(
  mat_treatment_sampletype_labeled,
  scale = "row",
  clustering_method = "complete",
  fontsize_row = 6,
  fontsize_col = 9,
  main = "Top 50 KO por Treatment e Sample.Type",
  filename = file.path(plots_dir, "heatmap_KO_top50_by_treatment_sampletype_annotated.png"),
  width = 14,
  height = 16
)

pheatmap(
  mat_all_groups_labeled,
  scale = "row",
  clustering_method = "complete",
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 50 KO por Soil, Treatment e Sample.Type",
  filename = file.path(plots_dir, "heatmap_KO_top50_by_all_groups_annotated.png"),
  width = 16,
  height = 18
)

message("")
message("Analise de KO/KEGG concluida.")
message("Tabelas: ", tables_dir)
message("Graficos: ", plots_dir)
message("")
message("Se quiser nomes dos KO, ative fetch_kegg_names <- TRUE e rode novamente.")
