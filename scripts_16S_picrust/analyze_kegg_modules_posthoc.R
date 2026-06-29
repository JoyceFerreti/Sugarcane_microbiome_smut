# ============================================================
# Pos-teste para modulos KEGG ja calculados
# Reutiliza resultados existentes de kegg_modules_analysis:
# - Dunn post-hoc para modulos significativos no Kruskal-Wallis
# - Heatmaps com simbolos nas celulas dos grupos envolvidos
#   em pelo menos uma comparacao significativa
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

install_if_missing(c("readr", "dplyr", "tidyr", "tibble", "stringr", "pheatmap", "rstatix", "multcompView"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(rstatix)
  library(multcompView)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

base_dir <- file.path(project_dir, "picrust2_results", "kegg_modules_analysis")
tables_dir <- file.path(base_dir, "tables")
plots_dir <- file.path(base_dir, "plots")
posthoc_dir <- file.path(base_dir, "posthoc")
dir.create(posthoc_dir, recursive = TRUE, showWarnings = FALSE)

module_clr_long_file <- file.path(tables_dir, "module_clr_long_with_metadata.csv")
module_group_means_file <- file.path(tables_dir, "module_group_means_clr.csv")
kruskal_roots_file <- file.path(tables_dir, "kruskal_kegg_modules_roots.csv")
kruskal_stalk_file <- file.path(tables_dir, "kruskal_kegg_modules_stalk.csv")

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Arquivo nao encontrado: ", path, call. = FALSE)
  }
}

stop_if_missing(module_clr_long_file)
stop_if_missing(module_group_means_file)
stop_if_missing(kruskal_roots_file)
stop_if_missing(kruskal_stalk_file)

module_clr_long <- read_csv(module_clr_long_file, show_col_types = FALSE)
module_group_means <- read_csv(module_group_means_file, show_col_types = FALSE)
kruskal_roots <- read_csv(kruskal_roots_file, show_col_types = FALSE)
kruskal_stalk <- read_csv(kruskal_stalk_file, show_col_types = FALSE)

pick_significant_modules <- function(kruskal_df) {
  mods <- kruskal_df %>%
    filter(!is.na(Kruskal_FDR) & Kruskal_FDR <= 0.05) %>%
    pull(Module)

  if (length(mods) == 0) {
    mods <- kruskal_df %>%
      filter(!is.na(Kruskal_p) & Kruskal_p <= 0.05) %>%
      pull(Module)
  }

  unique(mods)
}

run_dunn_for_compartment <- function(long_df, modules, compartment_pattern, compartment_name) {
  df_sub <- long_df %>%
    filter(
      str_detect(str_to_lower(Sample.Type), compartment_pattern),
      Module %in% modules
    ) %>%
    mutate(Group = paste(str_remove(Treatment, "^IAC-"), str_to_lower(Soil), sep = "_"))

  if (nrow(df_sub) == 0 || length(modules) == 0) {
    return(tibble())
  }

  out <- lapply(split(df_sub, df_sub$Module), function(x) {
    if (dplyr::n_distinct(x$Group) < 2) {
      return(NULL)
    }

    res <- tryCatch(
      rstatix::dunn_test(
        data = x,
        CLR_Abundance ~ Group,
        p.adjust.method = "BH"
      ),
      error = function(e) NULL
    )

    if (is.null(res) || nrow(res) == 0) {
      return(NULL)
    }

    res %>%
      mutate(
        Module = x$Module[1],
        Module_Name = x$Module_Name[1],
        Compartimento = compartment_name
      ) %>%
      select(Compartimento, Module, Module_Name, group1, group2, n1, n2, statistic, p, p.adj, p.adj.signif)
  })

  bind_rows(out)
}

make_module_label <- function(module, module_name) {
  ifelse(
    is.na(module_name) | module_name == "",
    module,
    paste0(module, " | ", module_name)
  )
}

build_heatmap_with_symbols <- function(group_means_df, dunn_df, compartment_pattern, title_txt, file_name) {
  if (is.null(dunn_df) || nrow(dunn_df) == 0) {
    return(invisible(NULL))
  }

  sig_modules <- dunn_df %>%
    filter(!is.na(p.adj) & p.adj <= 0.05) %>%
    pull(Module) %>%
    unique()

  if (length(sig_modules) == 0) {
    return(invisible(NULL))
  }

  heat_df <- group_means_df %>%
    filter(
      str_detect(SampleTypeClean, compartment_pattern),
      Module %in% sig_modules
    ) %>%
    mutate(
      GroupLabel = paste(str_remove(Treatment, "^IAC-"), str_to_lower(Soil), sep = "_"),
      ModuleLabel = make_module_label(Module, Module_Name)
    ) %>%
    select(Module, ModuleLabel, GroupLabel, Mean_CLR) %>%
    distinct()

  mat <- heat_df %>%
    select(ModuleLabel, GroupLabel, Mean_CLR) %>%
    pivot_wider(names_from = GroupLabel, values_from = Mean_CLR, values_fill = 0) %>%
    as.data.frame()

  rownames(mat) <- mat$ModuleLabel
  mat$ModuleLabel <- NULL
  mat <- as.matrix(mat)

  module_label_key <- heat_df %>%
    select(Module, ModuleLabel) %>%
    distinct()

  make_letters_for_module <- function(df_mod, group_levels) {
    group_levels <- unique(group_levels)
    if (length(group_levels) < 2) {
      return(tibble(GroupLabel = group_levels, Letters = "a"))
    }

    pmat <- matrix(1, nrow = length(group_levels), ncol = length(group_levels))
    rownames(pmat) <- group_levels
    colnames(pmat) <- group_levels

    for (i in seq_len(nrow(df_mod))) {
      g1 <- df_mod$group1[i]
      g2 <- df_mod$group2[i]
      pval <- df_mod$p.adj[i]

      if (is.na(g1) || is.na(g2) || is.na(pval)) {
        next
      }

      if (g1 %in% group_levels && g2 %in% group_levels) {
        pmat[g1, g2] <- pval
        pmat[g2, g1] <- pval
      }
    }

    letters_obj <- multcompView::multcompLetters(pmat, threshold = 0.05)

    tibble(
      GroupLabel = names(letters_obj$Letters),
      Letters = unname(letters_obj$Letters)
    )
  }

  display_df <- lapply(sig_modules, function(mod) {
    df_mod <- dunn_df %>% filter(Module == mod)
    groups_present <- heat_df %>%
      filter(Module == mod) %>%
      pull(GroupLabel) %>%
      unique()

    letters_df <- make_letters_for_module(df_mod, groups_present)
    module_label <- module_label_key %>%
      filter(Module == mod) %>%
      pull(ModuleLabel) %>%
      unique()

    tibble(
      Module = mod,
      ModuleLabel = module_label[1],
      GroupLabel = letters_df$GroupLabel,
      Symbol = letters_df$Letters
    )
  }) %>%
    bind_rows() %>%
    distinct()

  display_mat <- display_df %>%
    pivot_wider(names_from = GroupLabel, values_from = Symbol, values_fill = "") %>%
    as.data.frame()

  rownames(display_mat) <- display_mat$ModuleLabel
  display_mat$ModuleLabel <- NULL
  display_mat <- as.matrix(display_mat)

  display_mat <- display_mat[rownames(mat), colnames(mat), drop = FALSE]

  pheatmap(
    mat,
    scale = "row",
    clustering_method = "complete",
    display_numbers = display_mat,
    number_color = "black",
    fontsize_number = 10,
    fontsize_row = 6,
    fontsize_col = 10,
    main = title_txt,
    filename = file.path(plots_dir, file_name),
    width = 11,
    height = 14
  )
}

sig_roots <- pick_significant_modules(kruskal_roots)
sig_stalk <- pick_significant_modules(kruskal_stalk)

dunn_roots <- run_dunn_for_compartment(
  module_clr_long,
  sig_roots,
  "root",
  "Rhizosphere+Root"
)

dunn_stalk <- run_dunn_for_compartment(
  module_clr_long,
  sig_stalk,
  "stalk",
  "Stalk"
)

write_csv(dunn_roots, file.path(posthoc_dir, "dunn_kegg_modules_roots.csv"))
write_csv(dunn_stalk, file.path(posthoc_dir, "dunn_kegg_modules_stalk.csv"))

build_heatmap_with_symbols(
  module_group_means,
  dunn_roots,
  "root",
  "KEGG modules significativos com Dunn post-hoc - roots",
  "heatmap_kegg_modules_roots_dunn_symbols.png"
)

build_heatmap_with_symbols(
  module_group_means,
  dunn_stalk,
  "stalk",
  "KEGG modules significativos com Dunn post-hoc - stalk",
  "heatmap_kegg_modules_stalk_dunn_symbols.png"
)

message("")
message("Pos-teste de Dunn concluido.")
message("Tabelas: ", posthoc_dir)
message("Heatmaps com letras de agrupamento:")
message("- roots: ", file.path(plots_dir, "heatmap_kegg_modules_roots_dunn_symbols.png"))
message("- stalk: ", file.path(plots_dir, "heatmap_kegg_modules_stalk_dunn_symbols.png"))
message("")
message("Interpretacao das letras:")
message("Grupos que compartilham a mesma letra nao diferem entre si; letras diferentes indicam diferenca no pos-teste de Dunn (p.adj <= 0.05).")
