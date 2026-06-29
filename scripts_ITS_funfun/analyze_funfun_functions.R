# ============================================================
# Analise funcional de fungos a partir do FunFun
# Fluxo:
# - le Results.tsv do FunFun (Function x ASV)
# - le ASV_counts.csv (ASV x Sample)
# - combina as duas tabelas para gerar Function x Sample
# - junta com metadados
# - aplica CLR
# - separa roots e stalk
# - Kruskal-Wallis por funcao
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

funfun_file <- file.path(project_dir, "funfun_results", "Results.tsv")
counts_file <- file.path(project_dir, "ASV_counts.csv")
metadata_file <- file.path(project_dir, "mapp.csv")

analysis_dir <- file.path(project_dir, "funfun_analysis")
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

clean_function_name <- function(x) {
  x %>%
    str_replace_all("\\[.*?\\]", "") %>%
    str_squish()
}

plot_function_heatmap <- function(mat, title_txt, file_name) {
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))

  pheatmap(
    mat,
    scale = "row",
    clustering_method = "complete",
    fontsize_row = 6,
    fontsize_col = 10,
    main = title_txt,
    filename = file.path(plots_dir, file_name),
    width = 11,
    height = 14
  )
}

make_heatmap_matrix <- function(df, compartment_pattern, top_n = 40) {
  df_sub <- df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern))

  if (nrow(df_sub) == 0) return(NULL)

  top_features <- df_sub %>%
    group_by(Function) %>%
    summarise(VarAcrossGroups = var(Mean_CLR, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(VarAcrossGroups)) %>%
    slice_head(n = top_n) %>%
    pull(Function)

  df_sub %>%
    filter(Function %in% top_features) %>%
    select(Function, GroupLabel, Mean_CLR) %>%
    distinct() %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Function") %>%
    as.matrix()
}

make_significant_heatmap_matrix <- function(df_group_means, kruskal_df, compartment_pattern) {
  if (is.null(kruskal_df) || nrow(kruskal_df) == 0) return(NULL)

  sig_fun <- kruskal_df %>%
    filter(!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) %>%
    pull(Function)

  if (length(sig_fun) == 0) {
    sig_fun <- kruskal_df %>%
      filter(!is.na(Kruskal_p) & Kruskal_p <= 0.05) %>%
      pull(Function)
  }

  sig_fun <- unique(sig_fun)
  if (length(sig_fun) == 0) return(NULL)

  df_sub <- df_group_means %>%
    filter(
      str_detect(SampleTypeClean, compartment_pattern),
      Function %in% sig_fun
    ) %>%
    select(Function, GroupLabel, Mean_CLR) %>%
    distinct()

  if (nrow(df_sub) == 0) return(NULL)

  df_sub %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Function") %>%
    as.matrix()
}

kruskal_by_compartment <- function(df, compartment_pattern) {
  df_sub <- df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    mutate(Group = make_group_label(Treatment, Soil))

  if (nrow(df_sub) == 0) return(tibble())

  out <- lapply(split(df_sub, df_sub$Function), function(x) {
    n_groups <- dplyr::n_distinct(x$Group)
    if (n_groups < 2) return(NULL)

    test_res <- tryCatch(
      kruskal.test(CLR_Abundance ~ Group, data = x),
      error = function(e) NULL
    )

    if (is.null(test_res)) return(NULL)

    tibble(
      Function = x$Function[1],
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
      Function %in% sig_features
    ) %>%
    mutate(Group = make_group_label(Treatment, Soil))

  if (nrow(df_sub) == 0 || length(sig_features) == 0) return(tibble())

  out <- lapply(split(df_sub, df_sub$Function), function(x) {
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
        Function = x$Function[1],
        Compartimento = compartment_name
      ) %>%
      select(Compartimento, Function, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)
  })

  bind_rows(out)
}

# -----------------------------
# 3. Leitura dos dados
# -----------------------------
stop_if_missing(funfun_file)
stop_if_missing(counts_file)
stop_if_missing(metadata_file)

funfun_raw <- read_tsv(funfun_file, show_col_types = FALSE)
counts_raw <- read_delim(counts_file, delim = ";", show_col_types = FALSE, trim_ws = TRUE)
metadata <- read_delim(metadata_file, delim = ";", show_col_types = FALSE, trim_ws = TRUE)

names(funfun_raw)[1] <- "Function"
funfun_raw <- funfun_raw %>%
  mutate(Function = clean_function_name(Function))

names(counts_raw)[1] <- "ASV"

# -----------------------------
# 4. Preparar matriz funcional por amostra
# -----------------------------
funfun_df <- funfun_raw %>%
  mutate(across(-Function, ~replace_na(as.numeric(.x), 0)))

funfun_mat <- funfun_df %>%
  as.data.frame()

rownames(funfun_mat) <- funfun_mat$Function
funfun_mat$Function <- NULL
funfun_mat <- as.matrix(funfun_mat)

counts_df <- counts_raw %>%
  mutate(across(-ASV, as.numeric)) %>%
  as.data.frame()

rownames(counts_df) <- counts_df$ASV
counts_df$ASV <- NULL
counts_mat <- as.matrix(counts_df)

common_asvs <- intersect(colnames(funfun_mat), rownames(counts_mat))
if (length(common_asvs) == 0) {
  stop("Nenhum ASV em comum entre a saida do FunFun e a tabela de contagens.", call. = FALSE)
}

funfun_mat <- funfun_mat[, common_asvs, drop = FALSE]
counts_mat <- counts_mat[common_asvs, , drop = FALSE]

# Function x Sample = (Function x ASV) %*% (ASV x Sample)
function_sample_mat <- funfun_mat %*% counts_mat

function_sample_df <- as.data.frame(function_sample_mat) %>%
  rownames_to_column("Function")

write_csv(function_sample_df, file.path(tables_dir, "funfun_function_sample_matrix.csv"))

# -----------------------------
# 5. CLR por amostra
# -----------------------------
function_sample_clr <- clr_transform_matrix(function_sample_mat, pseudocount = pseudocount)

function_sample_clr_df <- as.data.frame(function_sample_clr) %>%
  rownames_to_column("Function")

write_csv(function_sample_clr_df, file.path(tables_dir, "funfun_function_sample_clr_matrix.csv"))

sample_cols <- colnames(function_sample_mat)
metadata_use <- metadata %>%
  filter(Sample %in% sample_cols) %>%
  mutate(
    GroupLabel = make_group_label(Treatment, Soil),
    SampleTypeClean = str_to_lower(Sample.Type)
  )

function_clr_long <- as.data.frame(function_sample_clr) %>%
  rownames_to_column("Function") %>%
  pivot_longer(
    cols = -Function,
    names_to = "Sample",
    values_to = "CLR_Abundance"
  ) %>%
  left_join(metadata_use, by = "Sample")

write_csv(function_clr_long, file.path(tables_dir, "funfun_function_clr_long_with_metadata.csv"))

function_group_means <- function_clr_long %>%
  group_by(Function, Sample.Type, SampleTypeClean, Treatment, Soil, GroupLabel) %>%
  summarise(Mean_CLR = mean(CLR_Abundance, na.rm = TRUE), .groups = "drop")

write_csv(function_group_means, file.path(tables_dir, "funfun_function_group_means_clr.csv"))

# -----------------------------
# 6. Heatmaps gerais
# -----------------------------
roots_mat <- make_heatmap_matrix(function_group_means, "root", top_n = top_n_heatmap)
stalk_mat <- make_heatmap_matrix(function_group_means, "stalk", top_n = top_n_heatmap)

plot_function_heatmap(
  roots_mat,
  "FunFun functions (CLR) - roots",
  "heatmap_funfun_functions_roots.png"
)

plot_function_heatmap(
  stalk_mat,
  "FunFun functions (CLR) - stalk",
  "heatmap_funfun_functions_stalk.png"
)

# -----------------------------
# 7. Kruskal-Wallis
# -----------------------------
kruskal_roots <- kruskal_by_compartment(function_clr_long, "root")
kruskal_stalk <- kruskal_by_compartment(function_clr_long, "stalk")

write_csv(kruskal_roots, file.path(tables_dir, "kruskal_funfun_functions_roots.csv"))
write_csv(kruskal_stalk, file.path(tables_dir, "kruskal_funfun_functions_stalk.csv"))

# -----------------------------
# 8. Heatmaps significativos
# -----------------------------
roots_sig_mat <- make_significant_heatmap_matrix(function_group_means, kruskal_roots, "root")
stalk_sig_mat <- make_significant_heatmap_matrix(function_group_means, kruskal_stalk, "stalk")

plot_function_heatmap(
  roots_sig_mat,
  "Significant FunFun functions (Kruskal-Wallis) - roots",
  "heatmap_funfun_functions_roots_significant.png"
)

plot_function_heatmap(
  stalk_sig_mat,
  "Significant FunFun functions (Kruskal-Wallis) - stalk",
  "heatmap_funfun_functions_stalk_significant.png"
)

# -----------------------------
# 9. Dunn post-hoc
# -----------------------------
sig_roots <- kruskal_roots %>%
  filter((!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) | (!is.na(Kruskal_p) & Kruskal_p <= 0.05)) %>%
  pull(Function) %>%
  unique()

sig_stalk <- kruskal_stalk %>%
  filter((!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) | (!is.na(Kruskal_p) & Kruskal_p <= 0.05)) %>%
  pull(Function) %>%
  unique()

dunn_roots <- run_dunn_for_compartment(function_clr_long, sig_roots, "root", "Rhizosphere+Root")
dunn_stalk <- run_dunn_for_compartment(function_clr_long, sig_stalk, "stalk", "Stalk")

write_csv(dunn_roots, file.path(posthoc_dir, "dunn_funfun_functions_roots.csv"))
write_csv(dunn_stalk, file.path(posthoc_dir, "dunn_funfun_functions_stalk.csv"))

# -----------------------------
# 10. Final
# -----------------------------
message("")
message("Analise funcional do FunFun concluida.")
message("Tabelas: ", tables_dir)
message("Graficos: ", plots_dir)
message("Pos-teste: ", posthoc_dir)
