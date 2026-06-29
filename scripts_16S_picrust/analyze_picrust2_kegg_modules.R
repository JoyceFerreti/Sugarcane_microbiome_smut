# ============================================================
# Analise de modulos KEGG a partir de KO (PICRUSt2)
# Gera:
# - mapeamento KO -> modulo KEGG
# - agregacao por modulo
# - CLR por amostra
# - heatmaps separados para roots e stalk
# - Kruskal-Wallis por modulo dentro de cada compartimento
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
install_if_missing(c("KEGGREST"), bioc = TRUE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(KEGGREST)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

metadata_file <- file.path(project_dir, "mapp_plant.csv")
ko_file <- file.path(project_dir, "picrust2_results", "summary_tables", "KO_predictions.csv")

analysis_dir <- file.path(project_dir, "picrust2_results", "kegg_modules_analysis")
tables_dir <- file.path(analysis_dir, "tables")
plots_dir <- file.path(analysis_dir, "plots")
cache_dir <- file.path(analysis_dir, "cache")

dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

top_n_modules <- 40
pseudocount <- 1

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

stop_if_missing(metadata_file)
stop_if_missing(ko_file)

clr_transform_matrix <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat)
  mat <- mat + pseudocount
  apply(mat, 2, function(x) {
    log_x <- log(x)
    log_x - mean(log_x, na.rm = TRUE)
  })
}

clean_ko_ids <- function(x) {
  x %>%
    str_trim() %>%
    str_remove("^ko:") %>%
    str_extract("K\\d{5}")
}

clean_module_ids <- function(x) {
  x %>%
    str_trim() %>%
    str_remove("^md:") %>%
    str_extract("M\\d{5}")
}

make_group_label <- function(treatment, soil) {
  genotype <- str_remove(treatment, "^IAC-")
  soil_clean <- str_to_lower(soil)
  paste(genotype, soil_clean, sep = "_")
}

get_kegg_module_map <- function(cache_dir) {
  link_cache <- file.path(cache_dir, "ko_to_module_map.csv")
  module_cache <- file.path(cache_dir, "module_names.csv")

  if (file.exists(link_cache) && file.exists(module_cache)) {
    ko_module_map <- read_csv(link_cache, show_col_types = FALSE)
    module_names <- read_csv(module_cache, show_col_types = FALSE)
    return(list(ko_module_map = ko_module_map, module_names = module_names))
  }

  message("Baixando mapeamento KO -> modulo do KEGG...")

  raw_links <- KEGGREST::keggLink("module", "ko")
  ko_module_map <- tibble(
    KO_raw = names(raw_links),
    Module_raw = unname(raw_links)
  ) %>%
    transmute(
      KO = clean_ko_ids(KO_raw),
      Module = clean_module_ids(Module_raw)
    ) %>%
    filter(!is.na(KO), !is.na(Module)) %>%
    distinct()

  raw_module_names <- KEGGREST::keggList("module")
  module_names <- tibble(
    Module_raw = names(raw_module_names),
    Module_Name = unname(raw_module_names)
  ) %>%
    transmute(
      Module = clean_module_ids(Module_raw),
      Module_Name = Module_Name
    ) %>%
    filter(!is.na(Module)) %>%
    distinct()

  write_csv(ko_module_map, link_cache)
  write_csv(module_names, module_cache)

  list(ko_module_map = ko_module_map, module_names = module_names)
}

metadata <- read_delim(metadata_file, delim = ";", show_col_types = FALSE, trim_ws = TRUE)
ko_tbl <- read_csv(ko_file, show_col_types = FALSE)

names(ko_tbl)[1] <- "KO"
ko_tbl <- ko_tbl %>%
  mutate(KO = clean_ko_ids(KO)) %>%
  filter(!is.na(KO))

sample_cols <- intersect(names(ko_tbl)[-1], metadata$Sample)

if (length(sample_cols) == 0) {
  stop("Nenhuma amostra da tabela KO foi encontrada nos metadados.", call. = FALSE)
}

metadata_use <- metadata %>%
  filter(Sample %in% sample_cols) %>%
  mutate(
    GroupLabel = make_group_label(Treatment, Soil),
    SampleTypeClean = str_to_lower(Sample.Type)
  )

kegg_maps <- get_kegg_module_map(cache_dir)
ko_module_map <- kegg_maps$ko_module_map
module_names <- kegg_maps$module_names

ko_long <- ko_tbl %>%
  select(KO, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample",
    values_to = "Abundance"
  )

ko_module_long <- ko_long %>%
  inner_join(ko_module_map, by = "KO") %>%
  left_join(module_names, by = "Module") %>%
  left_join(metadata_use, by = "Sample")

write_csv(ko_module_long, file.path(tables_dir, "KO_to_module_long_with_metadata.csv"))

module_sample_tbl <- ko_module_long %>%
  group_by(Module, Module_Name, Sample) %>%
  summarise(ModuleAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

module_sample_mat <- module_sample_tbl %>%
  select(Module, Sample, ModuleAbundance) %>%
  pivot_wider(names_from = Sample, values_from = ModuleAbundance, values_fill = 0) %>%
  as.data.frame()

module_names_vec <- module_sample_mat$Module
rownames(module_sample_mat) <- module_names_vec
module_sample_mat$Module <- NULL
module_sample_mat <- as.matrix(module_sample_mat)

module_sample_clr <- clr_transform_matrix(module_sample_mat, pseudocount = pseudocount)

module_sample_clr_df <- as.data.frame(module_sample_clr) %>%
  rownames_to_column("Module") %>%
  left_join(module_names, by = "Module")

write_csv(module_sample_clr_df, file.path(tables_dir, "module_sample_clr_matrix.csv"))

module_clr_long <- as.data.frame(module_sample_clr) %>%
  rownames_to_column("Module") %>%
  pivot_longer(
    cols = -Module,
    names_to = "Sample",
    values_to = "CLR_Abundance"
  ) %>%
  left_join(module_names, by = "Module") %>%
  left_join(metadata_use, by = "Sample")

write_csv(module_clr_long, file.path(tables_dir, "module_clr_long_with_metadata.csv"))

module_group_means <- module_clr_long %>%
  group_by(Module, Module_Name, Sample.Type, SampleTypeClean, Treatment, Soil, GroupLabel) %>%
  summarise(Mean_CLR = mean(CLR_Abundance, na.rm = TRUE), .groups = "drop")

write_csv(module_group_means, file.path(tables_dir, "module_group_means_clr.csv"))

make_module_label <- function(module, module_name) {
  ifelse(is.na(module_name) | module_name == "", module, paste0(module, " | ", module_name))
}

make_heatmap_matrix <- function(df, compartment_pattern, top_n = 40) {
  df_sub <- df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern))

  if (nrow(df_sub) == 0) {
    return(NULL)
  }

  top_modules <- df_sub %>%
    group_by(Module, Module_Name) %>%
    summarise(VarAcrossGroups = var(Mean_CLR, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(VarAcrossGroups)) %>%
    slice_head(n = top_n)

  df_sub %>%
    filter(Module %in% top_modules$Module) %>%
    mutate(ModuleLabel = make_module_label(Module, Module_Name)) %>%
    select(ModuleLabel, GroupLabel, Mean_CLR) %>%
    distinct() %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("ModuleLabel") %>%
    as.matrix()
}

plot_module_heatmap <- function(mat, title_txt, file_name) {
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    return(invisible(NULL))
  }

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

make_significant_heatmap_matrix <- function(df_group_means, kruskal_df, compartment_pattern) {
  if (is.null(kruskal_df) || nrow(kruskal_df) == 0) {
    return(NULL)
  }

  sig_modules <- kruskal_df %>%
    filter(!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) %>%
    pull(Module)

  if (length(sig_modules) == 0) {
    sig_modules <- kruskal_df %>%
      filter(!is.na(Kruskal_p) & Kruskal_p <= 0.05) %>%
      pull(Module)
  }

  sig_modules <- unique(sig_modules)

  if (length(sig_modules) == 0) {
    return(NULL)
  }

  df_sub <- df_group_means %>%
    filter(
      str_detect(SampleTypeClean, compartment_pattern),
      Module %in% sig_modules
    ) %>%
    mutate(ModuleLabel = make_module_label(Module, Module_Name)) %>%
    select(ModuleLabel, GroupLabel, Mean_CLR) %>%
    distinct()

  if (nrow(df_sub) == 0) {
    return(NULL)
  }

  df_sub %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("ModuleLabel") %>%
    as.matrix()
}

roots_mat <- make_heatmap_matrix(module_group_means, "root", top_n = top_n_modules)
stalk_mat <- make_heatmap_matrix(module_group_means, "stalk", top_n = top_n_modules)

plot_module_heatmap(
  roots_mat,
  "KEGG modules (CLR) - roots",
  "heatmap_kegg_modules_roots.png"
)

plot_module_heatmap(
  stalk_mat,
  "KEGG modules (CLR) - stalk",
  "heatmap_kegg_modules_stalk.png"
)

kruskal_by_compartment <- function(df, compartment_pattern) {
  df_sub <- df %>%
    filter(str_detect(SampleTypeClean, compartment_pattern)) %>%
    mutate(Group = make_group_label(Treatment, Soil))

  if (nrow(df_sub) == 0) {
    return(tibble())
  }

  out <- lapply(split(df_sub, df_sub$Module), function(x) {
    n_groups <- dplyr::n_distinct(x$Group)
    if (n_groups < 2) {
      return(NULL)
    }

    test_res <- tryCatch(
      kruskal.test(CLR_Abundance ~ Group, data = x),
      error = function(e) NULL
    )

    if (is.null(test_res)) {
      return(NULL)
    }

    tibble(
      Module = x$Module[1],
      Module_Name = x$Module_Name[1],
      Compartimento = unique(x$Sample.Type)[1],
      Groups_Compared = paste(sort(unique(x$Group)), collapse = "; "),
      Kruskal_p = unname(test_res$p.value)
    )
  })

  bind_rows(out) %>%
    mutate(Kruskal_FDR = p.adjust(Kruskal_p, method = "fdr")) %>%
    arrange(Kruskal_p)
}

kruskal_roots <- kruskal_by_compartment(module_clr_long, "root")
kruskal_stalk <- kruskal_by_compartment(module_clr_long, "stalk")

write_csv(kruskal_roots, file.path(tables_dir, "kruskal_kegg_modules_roots.csv"))
write_csv(kruskal_stalk, file.path(tables_dir, "kruskal_kegg_modules_stalk.csv"))

roots_sig_mat <- make_significant_heatmap_matrix(module_group_means, kruskal_roots, "root")
stalk_sig_mat <- make_significant_heatmap_matrix(module_group_means, kruskal_stalk, "stalk")

plot_module_heatmap(
  roots_sig_mat,
  "KEGG modules significativos (Kruskal-Wallis) - roots",
  "heatmap_kegg_modules_roots_significant.png"
)

plot_module_heatmap(
  stalk_sig_mat,
  "KEGG modules significativos (Kruskal-Wallis) - stalk",
  "heatmap_kegg_modules_stalk_significant.png"
)

message("")
message("Analise de modulos KEGG concluida.")
message("Tabelas: ", tables_dir)
message("Graficos: ", plots_dir)
message("")
message("Heatmaps gerados:")
message("- roots: ", file.path(plots_dir, "heatmap_kegg_modules_roots.png"))
message("- stalk: ", file.path(plots_dir, "heatmap_kegg_modules_stalk.png"))
message("- roots significativos: ", file.path(plots_dir, "heatmap_kegg_modules_roots_significant.png"))
message("- stalk significativos: ", file.path(plots_dir, "heatmap_kegg_modules_stalk_significant.png"))
message("")
message("Testes Kruskal-Wallis:")
message("- roots: ", file.path(tables_dir, "kruskal_kegg_modules_roots.csv"))
message("- stalk: ", file.path(tables_dir, "kruskal_kegg_modules_stalk.csv"))
