# ============================================================
# Analise estatistica de Pathways (PICRUSt2 - ITS)
# Fluxo:
# - le pathway_predictions.csv
# - junta com metadados
# - aplica CLR
# - separa roots e stalk
# - Kruskal-Wallis por pathway
# - Dunn post-hoc
# - heatmaps gerais e significativos
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

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr", "pheatmap", "rstatix"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(rstatix)
})

# -----------------------------
# 1. Configuracao
# -----------------------------
project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

metadata_file <- file.path(project_dir, "mapp.csv")
path_file <- file.path(project_dir, "picrust2_results", "summary_tables", "pathway_predictions.csv")

analysis_dir <- file.path(project_dir, "picrust2_results", "pathway_analysis")
tables_dir <- file.path(analysis_dir, "tables")
plots_dir <- file.path(analysis_dir, "plots")
posthoc_dir <- file.path(analysis_dir, "posthoc")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(posthoc_dir, recursive = TRUE, showWarnings = FALSE)

top_n_heatmap <- 40
pseudocount <- 1
group_levels <- c("5503_sandy", "5503_clay", "6007_sandy", "6007_clay")

# -----------------------------
# 2. Funcoes auxiliares
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

clr_transform_matrix <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat)
  mat <- mat + pseudocount
  apply(mat, 2, function(x) {
    log_x <- log(x)
    log_x - mean(log_x, na.rm = TRUE)
  })
}

make_group_label <- function(treatment, soil) {
  genotype <- str_remove(treatment, "^IAC-")
  soil_clean <- str_to_lower(soil)
  paste(genotype, soil_clean, sep = "_")
}

plot_pathway_heatmap <- function(mat, title_txt, file_name) {
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))

  pheatmap(
    mat,
    scale = "row",
    clustering_method = "complete",
    fontsize_row = 6,
    fontsize_col = 10,
    main = title_txt,
    filename = file.path(plots_dir, file_name),
    width = 10,
    height = 14
  )
}

make_heatmap_matrix <- function(df, compartment_pattern, top_n = 40) {
  df_sub <- df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern))

  if (nrow(df_sub) == 0) return(NULL)

  top_features <- df_sub %>%
    group_by(Pathway) %>%
    summarise(VarAcrossGroups = var(Mean_CLR, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(VarAcrossGroups)) %>%
    slice_head(n = top_n) %>%
    pull(Pathway)

  df_sub %>%
    filter(Pathway %in% top_features) %>%
    select(Pathway, GroupLabel, Mean_CLR) %>%
    distinct() %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Pathway") %>%
    as.matrix()
}

make_significant_heatmap_matrix <- function(df_group_means, kruskal_df, compartment_pattern) {
  if (is.null(kruskal_df) || nrow(kruskal_df) == 0) return(NULL)

  sig_path <- kruskal_df %>%
    filter(!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) %>%
    pull(Pathway)

  if (length(sig_path) == 0) {
    sig_path <- kruskal_df %>%
      filter(!is.na(Kruskal_p) & Kruskal_p <= 0.05) %>%
      pull(Pathway)
  }

  sig_path <- unique(sig_path)
  if (length(sig_path) == 0) return(NULL)

  df_sub <- df_group_means %>%
    filter(
      str_detect(SampleTypeClean, compartment_pattern),
      Pathway %in% sig_path
    ) %>%
    select(Pathway, GroupLabel, Mean_CLR) %>%
    distinct()

  if (nrow(df_sub) == 0) return(NULL)

  df_sub %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Pathway") %>%
    as.matrix()
}

kruskal_by_compartment <- function(df, compartment_pattern) {
  df_sub <- df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    mutate(Group = make_group_label(Treatment, Soil))

  if (nrow(df_sub) == 0) return(tibble())

  out <- lapply(split(df_sub, df_sub$Pathway), function(x) {
    n_groups <- dplyr::n_distinct(x$Group)
    if (n_groups < 2) return(NULL)

    test_res <- tryCatch(
      kruskal.test(CLR_Abundance ~ Group, data = x),
      error = function(e) NULL
    )

    if (is.null(test_res)) return(NULL)

    tibble(
      Pathway = x$Pathway[1],
      Compartimento = unique(x$Sample.Type)[1],
      Groups_Compared = paste(sort(unique(x$Group)), collapse = "; "),
      Kruskal_p = unname(test_res$p.value)
    )
  })

  bind_rows(out) %>%
    mutate(Kruskal_FDR = p.adjust(Kruskal_p, method = "fdr")) %>%
    arrange(Kruskal_p)
}

run_dunn_for_compartment <- function(long_df, sig_features, compartment_pattern, compartment_name) {
  df_sub <- long_df %>%
    filter(
      str_detect(SampleTypeClean, compartment_pattern),
      Pathway %in% sig_features
    ) %>%
    mutate(Group = make_group_label(Treatment, Soil))

  if (nrow(df_sub) == 0 || length(sig_features) == 0) return(tibble())

  out <- lapply(split(df_sub, df_sub$Pathway), function(x) {
    if (dplyr::n_distinct(x$Group) < 2) return(NULL)

    res <- tryCatch(
      rstatix::dunn_test(
        data = x,
        CLR_Abundance ~ Group,
        p.adjust.method = "BH"
      ),
      error = function(e) NULL
    )

    if (is.null(res) || nrow(res) == 0) return(NULL)

    res %>%
      mutate(
        Pathway = x$Pathway[1],
        Compartimento = compartment_name
      ) %>%
      select(Compartimento, Pathway, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)
  })

  bind_rows(out)
}

# -----------------------------
# 3. Leitura
# -----------------------------
stop_if_missing(metadata_file)
stop_if_missing(path_file)

metadata <- read_delim(metadata_file, delim = ";", show_col_types = FALSE, trim_ws = TRUE)
path_tbl <- read_csv(path_file, show_col_types = FALSE)

names(path_tbl)[1] <- "Pathway"

sample_cols <- intersect(names(path_tbl)[-1], metadata$Sample)
if (length(sample_cols) == 0) {
  stop("Nenhuma coluna de amostra da tabela de pathways foi encontrada nos metadados.", call. = FALSE)
}

metadata_use <- metadata %>%
  filter(Sample %in% sample_cols) %>%
  mutate(
    GroupLabel = make_group_label(Treatment, Soil),
    SampleTypeClean = str_to_lower(Sample.Type)
  )

# -----------------------------
# 4. CLR por amostra
# -----------------------------
path_sample_mat <- path_tbl %>%
  select(Pathway, all_of(sample_cols)) %>%
  as.data.frame()

rownames(path_sample_mat) <- path_sample_mat$Pathway
path_sample_mat$Pathway <- NULL
path_sample_mat <- as.matrix(path_sample_mat)

path_sample_clr <- clr_transform_matrix(path_sample_mat, pseudocount = pseudocount)

path_sample_clr_df <- as.data.frame(path_sample_clr) %>%
  rownames_to_column("Pathway")

write_csv(path_sample_clr_df, file.path(tables_dir, "pathway_sample_clr_matrix.csv"))

path_clr_long <- as.data.frame(path_sample_clr) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(
    cols = -Pathway,
    names_to = "Sample",
    values_to = "CLR_Abundance"
  ) %>%
  left_join(metadata_use, by = "Sample")

write_csv(path_clr_long, file.path(tables_dir, "pathway_clr_long_with_metadata.csv"))

path_group_means <- path_clr_long %>%
  group_by(Pathway, Sample.Type, SampleTypeClean, Treatment, Soil, GroupLabel) %>%
  summarise(Mean_CLR = mean(CLR_Abundance, na.rm = TRUE), .groups = "drop")

write_csv(path_group_means, file.path(tables_dir, "pathway_group_means_clr.csv"))

# -----------------------------
# 5. Heatmaps gerais
# -----------------------------
roots_mat <- make_heatmap_matrix(path_group_means, "root", top_n = top_n_heatmap)
stalk_mat <- make_heatmap_matrix(path_group_means, "stalk", top_n = top_n_heatmap)

plot_pathway_heatmap(
  roots_mat,
  "Pathways (CLR) - roots",
  "heatmap_pathways_roots.png"
)

plot_pathway_heatmap(
  stalk_mat,
  "Pathways (CLR) - stalk",
  "heatmap_pathways_stalk.png"
)

# -----------------------------
# 6. Kruskal-Wallis
# -----------------------------
kruskal_roots <- kruskal_by_compartment(path_clr_long, "root")
kruskal_stalk <- kruskal_by_compartment(path_clr_long, "stalk")

write_csv(kruskal_roots, file.path(tables_dir, "kruskal_pathways_roots.csv"))
write_csv(kruskal_stalk, file.path(tables_dir, "kruskal_pathways_stalk.csv"))

# -----------------------------
# 7. Heatmaps significativos
# -----------------------------
roots_sig_mat <- make_significant_heatmap_matrix(path_group_means, kruskal_roots, "root")
stalk_sig_mat <- make_significant_heatmap_matrix(path_group_means, kruskal_stalk, "stalk")

plot_pathway_heatmap(
  roots_sig_mat,
  "Significant pathways (Kruskal-Wallis) - roots",
  "heatmap_pathways_roots_significant.png"
)

plot_pathway_heatmap(
  stalk_sig_mat,
  "Significant pathways (Kruskal-Wallis) - stalk",
  "heatmap_pathways_stalk_significant.png"
)

# -----------------------------
# 8. Dunn post-hoc
# -----------------------------
sig_roots <- kruskal_roots %>%
  filter((!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) | (!is.na(Kruskal_p) & Kruskal_p <= 0.05)) %>%
  pull(Pathway) %>%
  unique()

sig_stalk <- kruskal_stalk %>%
  filter((!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) | (!is.na(Kruskal_p) & Kruskal_p <= 0.05)) %>%
  pull(Pathway) %>%
  unique()

dunn_roots <- run_dunn_for_compartment(path_clr_long, sig_roots, "root", "Rhizosphere+Root")
dunn_stalk <- run_dunn_for_compartment(path_clr_long, sig_stalk, "stalk", "Stalk")

write_csv(dunn_roots, file.path(posthoc_dir, "dunn_pathways_roots.csv"))
write_csv(dunn_stalk, file.path(posthoc_dir, "dunn_pathways_stalk.csv"))

# -----------------------------
# 9. Final
# -----------------------------
message("")
message("Analise estatistica de pathways concluida.")
message("Tabelas: ", tables_dir)
message("Graficos: ", plots_dir)
message("Pos-teste: ", posthoc_dir)
