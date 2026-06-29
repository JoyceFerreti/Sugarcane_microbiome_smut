# ============================================================
# PICRUSt2 completo para ITS/fungos
# Projeto: ITS / ASVs
# Arquivos esperados nesta pasta:
#   - ASVs.fa
#   - ASV_counts.csv
#   - mapp.csv
#   - ASV_tax.csv (opcional)
#
# Fluxo:
# - prepara tabela BIOM
# - roda PICRUSt2
# - importa EC e pathways
# - salva tabelas-resumo e graficos simples
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

install_if_missing(c("readr", "dplyr", "tibble", "stringr", "ggplot2", "tidyr"))
install_if_missing(c("biomformat"), bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(biomformat)
})

# -----------------------------
# 1. Configuracao
# -----------------------------
project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

counts_file   <- file.path(project_dir, "ASV_counts.csv")
metadata_file <- file.path(project_dir, "mapp.csv")
fasta_file    <- file.path(project_dir, "ASVs.fa")
taxonomy_file <- file.path(project_dir, "ASV_tax.csv")

output_dir   <- file.path(project_dir, "picrust2_results")
input_dir    <- file.path(output_dir, "input_prepared")
pipeline_dir <- file.path(output_dir, "pipeline_out")
summary_dir  <- file.path(output_dir, "summary_tables")
plots_dir    <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

runner_mode <- "local"
threads <- 8
run_picrust2 <- TRUE

# -----------------------------
# 2. Funcoes auxiliares
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

read_fasta_ids <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  headers <- x[startsWith(x, ">")]
  sub("^>", "", headers)
}

write_tsv_plain <- function(df, path) {
  write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

run_cmd <- function(command, args = character()) {
  output <- suppressWarnings(system2(command, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  if (!identical(status, 0L)) {
    stop(
      paste0(
        "Falha ao executar comando: ", command, "\n",
        paste(output, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  invisible(output)
}

check_pipeline_outputs <- function(out_dir) {
  expected_files <- c(
    file.path(out_dir, "EC_metagenome_out", "pred_metagenome_unstrat.tsv.gz"),
    file.path(out_dir, "pathways_out", "path_abun_unstrat.tsv.gz")
  )

  missing_files <- expected_files[!file.exists(expected_files)]

  if (length(missing_files) > 0) {
    stop(
      paste0(
        "O comando do PICRUSt2 terminou, mas os arquivos finais nao foram encontrados.\n",
        "Arquivos ausentes:\n",
        paste(missing_files, collapse = "\n")
      ),
      call. = FALSE
    )
  }
}

run_picrust2_pipeline <- function(table_biom, seqs_fasta, out_dir, threads) {
  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }

  args <- c(
    "-s", shQuote(seqs_fasta),
    "-i", shQuote(table_biom),
    "-o", shQuote(out_dir),
    "-p", as.character(threads)
  )

  run_cmd("picrust2_pipeline.py", args)
  invisible(NULL)
}

read_picrust_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  readr::read_tsv(path, show_col_types = FALSE)
}

top_features_long <- function(df, metadata, label, n = 20) {
  if (is.null(df)) return(NULL)

  feature_col <- names(df)[1]
  sample_cols <- intersect(names(df)[-1], metadata$Sample)
  if (length(sample_cols) == 0) return(NULL)

  long_df <- df %>%
    select(all_of(c(feature_col, sample_cols))) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "Abundance"
    ) %>%
    left_join(metadata, by = "Sample") %>%
    rename(Feature = all_of(feature_col)) %>%
    group_by(Feature, Treatment, Sample.Type, Soil) %>%
    summarise(MeanAbundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

  top_features <- long_df %>%
    group_by(Feature) %>%
    summarise(GlobalMean = mean(MeanAbundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(GlobalMean)) %>%
    slice_head(n = n) %>%
    pull(Feature)

  long_df %>%
    filter(Feature %in% top_features) %>%
    mutate(ResultType = label)
}

save_feature_plot <- function(df, file_name, title_txt) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  p <- ggplot(df, aes(x = Treatment, y = MeanAbundance, fill = Soil)) +
    geom_col(position = "dodge") +
    facet_wrap(~ Feature, scales = "free_y") +
    theme_bw(base_size = 11) +
    labs(
      title = title_txt,
      x = "Treatment",
      y = "Predicted mean abundance"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 7)
    )

  ggsave(
    filename = file.path(plots_dir, file_name),
    plot = p,
    width = 15,
    height = 10,
    dpi = 300
  )
}

# -----------------------------
# 3. Leitura e validacao
# -----------------------------
stop_if_missing(counts_file)
stop_if_missing(metadata_file)
stop_if_missing(fasta_file)

counts_raw <- read_delim(
  counts_file,
  delim = ";",
  show_col_types = FALSE,
  trim_ws = TRUE
)

metadata <- read_delim(
  metadata_file,
  delim = ";",
  show_col_types = FALSE,
  trim_ws = TRUE
)

taxonomy <- NULL
if (file.exists(taxonomy_file)) {
  taxonomy <- read_delim(
    taxonomy_file,
    delim = ";",
    show_col_types = FALSE,
    trim_ws = TRUE
  )
}

if (!"ASV" %in% names(counts_raw)) {
  stop("A tabela de contagens precisa ter uma coluna chamada 'ASV'.", call. = FALSE)
}

if (!"Sample" %in% names(metadata)) {
  stop("A tabela de metadados precisa ter uma coluna chamada 'Sample'.", call. = FALSE)
}

sample_cols <- setdiff(names(counts_raw), "ASV")

missing_in_metadata <- setdiff(sample_cols, metadata$Sample)
if (length(missing_in_metadata) > 0) {
  stop(
    "As seguintes amostras da tabela de contagens nao estao nos metadados: ",
    paste(missing_in_metadata, collapse = ", "),
    call. = FALSE
  )
}

fasta_ids <- read_fasta_ids(fasta_file)
count_ids <- counts_raw$ASV

missing_in_fasta <- setdiff(count_ids, fasta_ids)
if (length(missing_in_fasta) > 0) {
  stop(
    "Existem ASVs na tabela de contagens que nao foram encontradas no FASTA: ",
    paste(head(missing_in_fasta, 20), collapse = ", "),
    call. = FALSE
  )
}

# -----------------------------
# 4. Preparar tabela BIOM
# -----------------------------
counts_matrix <- counts_raw %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  as.data.frame()

rownames(counts_matrix) <- counts_matrix$ASV
counts_matrix$ASV <- NULL

feature_table <- as.matrix(counts_matrix)

if (any(is.na(feature_table))) {
  stop("Foram encontrados valores NA na tabela de contagens.", call. = FALSE)
}

biom_obj <- biomformat::make_biom(data = feature_table)
biom_file <- file.path(input_dir, "asv_table.biom")
biomformat::write_biom(biom_obj, biom_file)

prepared_tsv <- file.path(input_dir, "asv_table.tsv")
prepared_df <- tibble::rownames_to_column(as.data.frame(feature_table), var = "ASV")
write_tsv_plain(prepared_df, prepared_tsv)

metadata_out <- file.path(input_dir, "metadata_clean.tsv")
write_tsv_plain(metadata, metadata_out)

if (!is.null(taxonomy)) {
  taxonomy_out <- file.path(input_dir, "taxonomy.tsv")
  write_tsv_plain(taxonomy, taxonomy_out)
}

message("Arquivos preparados em: ", input_dir)

# -----------------------------
# 5. Rodar PICRUSt2
# -----------------------------
if (run_picrust2) {
  message("Iniciando PICRUSt2 para ITS...")
  run_picrust2_pipeline(
    table_biom = biom_file,
    seqs_fasta = fasta_file,
    out_dir = pipeline_dir,
    threads = threads
  )
  check_pipeline_outputs(pipeline_dir)
  message("PICRUSt2 finalizado.")
} else {
  message("Execucao do PICRUSt2 pulada. Reutilizando resultados existentes.")
}

# -----------------------------
# 6. Importar resultados
# -----------------------------
ec_file <- file.path(pipeline_dir, "EC_metagenome_out", "pred_metagenome_unstrat.tsv.gz")
path_file <- file.path(pipeline_dir, "pathways_out", "path_abun_unstrat.tsv.gz")

ec_tbl <- read_picrust_table(ec_file)
path_tbl <- read_picrust_table(path_file)

if (!is.null(ec_tbl)) {
  write_csv(ec_tbl, file.path(summary_dir, "EC_predictions.csv"))
}

if (!is.null(path_tbl)) {
  write_csv(path_tbl, file.path(summary_dir, "pathway_predictions.csv"))
}

# -----------------------------
# 7. Resumos simples
# -----------------------------
ec_long <- top_features_long(ec_tbl, metadata, "EC", n = 20)
path_long <- top_features_long(path_tbl, metadata, "Pathway", n = 20)

if (!is.null(ec_long)) {
  write_csv(ec_long, file.path(summary_dir, "EC_top20_group_means.csv"))
}

if (!is.null(path_long)) {
  write_csv(path_long, file.path(summary_dir, "Pathway_top20_group_means.csv"))
}

save_feature_plot(ec_long, "EC_top20_by_group.png", "Top 20 predicted ECs by group")
save_feature_plot(path_long, "Pathway_top20_by_group.png", "Top 20 predicted pathways by group")

# -----------------------------
# 8. Final
# -----------------------------
message("")
message("Processo concluido.")
message("Entrada preparada: ", input_dir)
message("Saida do pipeline: ", pipeline_dir)
message("Tabelas-resumo: ", summary_dir)
message("Graficos: ", plots_dir)
