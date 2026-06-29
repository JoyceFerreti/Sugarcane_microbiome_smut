# 16S / PICRUSt2 Pipeline Scripts

This folder consolidates the main scripts used in the bacterial functional prediction workflow.

## Recommended execution order

1. `picrust2_run_complete.R`
   - Prepares the 16S inputs and runs PICRUSt2.

2. `analyze_picrust2_kegg.R`
   - General KEGG KO-level exploration.

3. `analyze_picrust2_kegg_modules.R`
   - Maps `KO -> KEGG modules`, aggregates module abundances, applies CLR, and generates the core module tables.

4. `analyze_kegg_modules_posthoc.R`
   - Runs Dunn post hoc tests for significant bacterial KEGG modules.

5. `classify_kegg_modules_upset.R`
   - Classifies modules into biological response categories and generates UpSet summaries.

## Visualization scripts

- `plot_kegg_modules_category_heatmaps.R`
  - Heatmaps by category.

- `plot_kegg_modules_focus_heatmaps.R`
  - Focused heatmaps using selected top modules.

- `plot_kegg_modules_focus_heatmaps_with_letters.R`
  - Focused heatmaps with Dunn grouping letters.

- `plot_kegg_modules_final_figure.R`
  - Final bacterial figure with prioritized categories.

- `plot_kegg_modules_final_heatmap.R`
  - Alternative final heatmap version.

- `plot_thesis_bacteria_fungi_panels.R`
  - Thesis-ready heatmap panels for bacteria and fungi.

## Additional analysis

- `analyze_picrust2_16s_ec.R`
  - Exploratory EC-based analysis for the bacterial dataset.

## Notes

- The main bacterial functional analysis used in the thesis is the **KEGG module-based workflow**.
- Category assignment for focused figures was used for visualization and does not replace the complete multi-category classification retained in the output tables.
