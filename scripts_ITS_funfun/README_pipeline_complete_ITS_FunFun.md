# ITS / FunFun Pipeline Scripts

This folder consolidates the main scripts used in the fungal functional prediction workflow.

## Recommended final pipeline

1. `analyze_funfun_functions.R`
   - Reads the FunFun `Results.tsv`, combines function-by-ASV predictions with ASV abundance, reconstructs sample-level functional tables, applies CLR, and performs general exploratory analyses.

2. `analyze_funfun_kegg_pathways.R`
   - Retains only `PATH:ko...` annotations, reconstructs KEGG pathway abundance per sample, applies CLR, and runs Kruskal-Wallis and Dunn tests.

3. `classify_funfun_pathways_exploratory.R`
   - Classifies fungal KEGG pathways into biological response categories, generates UpSet plots, and creates focused heatmaps with Dunn grouping letters.

## Additional scripts retained for record

- `picrust2_its_run_complete.R`
  - Initial ITS PICRUSt2 attempt.

- `analyze_picrust2_its_ec.R`
  - EC-based analysis from the initial ITS PICRUSt2 attempt.

- `analyze_picrust2_its_pathways.R`
  - Pathway-based analysis from the initial ITS PICRUSt2 attempt.

- `annotate_picrust2_its_ec_top50.R`
  - Annotation of top EC features from the initial ITS PICRUSt2 branch.

## Notes

- The final fungal functional analysis used in the thesis is the **FunFun-based KEGG pathway workflow**.
- The earlier ITS PICRUSt2 scripts were retained for documentation, but the final interpretation relied on FunFun outputs.
- Fungal pathway classification was treated as exploratory in the final figures and text.
