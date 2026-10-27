---
name: sc-figure-plan
description: Plan and render the standard figure set for a single-cell RNA-seq study (QC, UMAPs, composition, DEGs, enrichment, modules, TF activity, cell-cell communication) with replicate-aware statistics. Use for any request to plot, visualize or make publication figures from single-cell data.
metadata:
  authored_via: "claude-science"
  published_by: "user"
  published_at: "2026-10-01T09:02:40.969Z"
  authored_by: "user"
  authored_at: "2026-10-01T09:02:29.150Z"
  last_modified_by: "user"
  last_modified_at: "2026-10-01T09:02:40.951Z"
---
# Single-cell figure plan

Turn an annotated single-cell dataset and its analysis results into a consistent, publication-ready figure set. This skill covers what to plot, how to encode it, and what to report. It does not cover how to run the upstream analyses (integration, annotation, DE, NMF, TF inference, communication); use the matching analysis skill or the user's existing results for those.

## Workflow

1. **Inspect the data before planning.** List `obs` columns and map them to roles: `cell_type`, `sample` (the biological replicate), `group`/`condition`, `time`, `tissue`/`source`, `batch`, `donor`. Count samples per group. Check which result tables already exist (DE, enrichment, module scores, TF activity, communication). Ask the user only when a role is ambiguous, most often which column is the biological replicate.
2. **Choose the figures.** Go through the catalog below and keep only figures whose prerequisites are met. Tell the user which were skipped and why. A figure forced onto a design that cannot support it (for example a forest plot with one sample per group) is misleading rather than incomplete.
3. **Share the plan** as a short list (figure, panels, variables, statistic) and proceed unless the user objects.
4. **Render** with the global rules below. Do not recompute an analysis that already exists; plot the existing result so figures match the text.
5. **Look at every rendered figure** before delivering: overlapping labels, clipped legends, unreadable text, empty panels, wrong colour-scale limits. Fix and re-render.
6. **Deliver** per figure: vector PDF, PNG preview, source-data CSV, and a legend entry (see Legend checklist).

## Global rules

**Statistical units.** Cells are not replicates. Keep cells, samples and biological replicates distinct in every figure. Any group comparison uses the biological replicate as the unit (pseudobulk or per-sample summaries), because cell-level tests treat thousands of correlated cells as independent and inflate significance. With fewer than 3 replicates per group, show the data descriptively and say no inference was done.

**Colour.**
- Build one palette per categorical variable (cell type, group, tissue, batch) once, store it (for example `adata.uns['<key>_colors']` or a JSON file), and reuse it in every figure. A cell type that changes colour between figures costs the reader more than any other styling flaw.
- Use colourblind-safe categorical palettes; with more than about 20 categories, also label directly on the plot.
- Diverging blue-white-red only for centred quantities (z-scored expression, log2FC, activity scores), symmetric around 0. For non-negative quantities (raw or log-normalized expression, counts, -log10 p) use a sequential map; a diverging map there invents a midpoint that means nothing.
- Down = blue, up = red, not significant = grey, everywhere.

**Typography and layout.**
- Arial throughout. Check that it is actually installed (`matplotlib.font_manager.findfont('Arial', fallback_to_default=False)`); plotting libraries fall back silently. If absent, use Helvetica or Liberation Sans (metric-compatible) and tell the user.
- Text 5-8 pt at final print size. Design at final size (single column about 89 mm, double about 183 mm) rather than scaling afterwards.
- Adjust width or height to the content: scale dot plot and heatmap dimensions with the number of rows and columns so labels never overlap; keep UMAP panels square with equal aspect.
- Choose single-panel or multi-panel layouts by study design. Split panels share axes and colour scales unless stated.

**Embeddings.**
- Plot cells in random order, otherwise the last-drawn category hides the others.
- Rasterize the point layer (300+ dpi) and keep axes, text and legends as vectors, so PDFs stay small and editable.
- In split UMAPs, draw all other cells in light grey behind the highlighted subset for context.
- Scale point size to cell number.

**Export.** PDF with editable text (`pdf.fonttype = 42`) plus a 300 dpi PNG. Name files `Fig01_qc_overview.pdf`, `Fig02_celltypes_markers.pdf`, and so on. Save the plotted values as `FigNN_source_data.csv`.

Minimal matplotlib setup:

```python
import matplotlib as mpl
mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans'],
    'font.size': 7, 'axes.titlesize': 8, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.spines.top': False, 'axes.spines.right': False,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
})
```

The same rules apply in R (ggplot2/Seurat): `theme(text = element_text(family = 'Arial'))`, `ggrastr::rasterise()` for point layers, `cairo_pdf` for export.

## Figure catalog

Each entry: what it needs, what to draw, what to report.

### 1. Sample quality-control overview
- **Needs:** per-sample cell counts before and after QC.
- **Draw:** paired bars (before, after) per sample with the retained percentage labelled. Order samples by group, source or batch, and mark that variable with bar colour or an annotation strip under the x-axis. If per-cell metrics are available, add violins of genes per cell, UMIs per cell and mitochondrial fraction per sample.
- **Report:** QC thresholds and doublet method.

### 2. Cell types and marker genes
- **Needs:** embedding, cell-type labels, marker list.
- **Draw:** overall UMAP coloured by cell type, plus UMAPs split by sample source, tissue region or another key variable. Dot plot of markers by cell type: colour = mean expression scaled per gene (blue-to-red, centred), size = fraction of cells expressing. Order cell types by lineage and group markers by cell type so the diagonal reads cleanly. To check that markers hold across samples, repeat the dot plot per sample or stratum.
- **Report:** normalization, scaling, and how markers were chosen.

### 3. Sample distribution and cell composition
- **Needs:** sample labels; two or more samples.
- **Draw:** overall and per-sample UMAPs coloured by sample, to show how samples spread across the atlas (and whether integration left batch structure). Percentage-stacked bars of cell-type composition with one bar per sample, grouped and annotated by experimental group, source or batch. Prefer per-sample bars over pooled per-group bars, since pooling hides replicate variability and lets large samples dominate.
- **Report:** the denominator of each proportion (all cells in the sample, or cells within a compartment).

### 4. Composition changes across groups and conditions
- **Needs:** two or more groups. Forest plot needs 3 or more biological replicates per group.
- **Draw:** overall UMAP and UMAPs split by group or condition. Optional embedding density plots; subsample groups to equal cell numbers first, otherwise density reflects sequencing depth of the group rather than biology. Stacked bars of composition across groups or time points. With enough replicates, test proportion differences with a replicate-aware compositional method (for example propeller, sccomp, scCODA, or a mixed model on per-sample proportions) and show effect estimates with confidence intervals in a forest plot, one row per cell type, vertical line at no effect.
- **Report:** method, n per group, effect-size definition, interval level, multiple-testing correction.

### 5. Differential-expression overview
- **Needs:** DE results per cell type and prespecified contrast, ideally pseudobulk by sample.
- **Draw:** heatmap of DEG counts, rows = cell types, columns = contrasts, faceted by tissue or source; up and down shown separately (two panels or split cells), counts printed in cells. Mark untested combinations as missing, not zero. For selected contrasts, volcano plots: down blue, up red, not significant grey, threshold lines drawn, top five genes per direction labelled with non-overlapping labels.
- **Report:** DE method and unit, adjusted p threshold, log2FC threshold, expression filters, identical across all contrasts. State that DEG counts depend on cell and replicate numbers, so a count difference between cell types is not by itself a difference in response; annotate cell or sample numbers alongside when they vary widely.

### 6. Functional enrichment of DEGs
- **Needs:** DEG lists; a defined background.
- **Draw:** run up- and down-regulated genes separately per contrast. Dot plot of terms by cell type: one encoding for colour (for example -log10 adjusted p or enrichment score) and one for size (gene count or gene ratio), both defined in the legend. Cap the terms shown (for example top 5 per cell type by adjusted p), collapse near-duplicate terms, and make separate plots per tissue or source as needed.
- **Report:** gene-set database and version, background (the genes tested in that cell type, not the whole genome), correction method, term-selection rule.

### 7. Functional modules by cell type
- **Needs:** module decomposition (for example NMF) run within each cell type.
- **Draw:** clustered heatmap of cells (or samples) against module usage, with annotation bars for sample and group. Rank-selection plot (reconstruction error, stability or cophenetic coefficient against number of modules) as a supplementary panel that justifies the chosen number. Annotate modules by enrichment of top genes (bar plots) and gene-set activity scores such as AUCell (heatmap of modules or populations against gene sets).
- **Report:** input genes, rank-selection criterion, how top module genes were defined, scoring method.

### 8. Condition-dependent module trends
- **Needs:** module scores per sample; trend classes need 3 or more ordered conditions or time points. With two conditions, reduce to increase or decrease.
- **Draw:** per cell type and module, line plot across conditions with one point per biological replicate and mean with error bars, faceted by group. Summarize trend classes (sustained increase, sustained decrease, increase then decrease, decrease then increase, other) in a tile plot of cell type against module.
- **Report:** scoring method, error-bar definition, n, and the explicit rule that assigns each trend class, applied to replicate-level summaries.

### 9. Transcription-factor activity
- **Needs:** TF activity estimates per cell type or condition.
- **Draw:** dot plot of TFs against cell types (or conditions): colour = activity score (diverging, centred), size = significance or fraction of cells with active regulon. Keep score and significance as separate encodings so a large score is not read as a significant one. Show cell-type-specific or condition-associated TFs chosen by a stated rule.
- **Report:** inference method and regulon resource with version, TF selection rule.

### 10. TF, module and DE integration network
- **Needs:** figures 5, 7 and 9 results.
- **Draw:** network per cell type or contrast. TF nodes coloured by activity change; target nodes coloured by log2FC and grouped or outlined by module; distinct shapes for TFs and targets; edge style for activation or repression. Limit to a few key TFs and a bounded number of targets per panel so every label is legible; a hairball conveys nothing.
- **Report:** edge source, edge and node filtering criteria, how changes are encoded, and the prespecified ranking used to pick key TFs.

### 11. Cell-cell communication within a tissue or source
- **Needs:** communication inference run separately per tissue, source or prespecified subset.
- **Draw:** dot plot with sender-receiver pairs on one axis and ligand-receptor pairs on the other; colour = interaction strength, size = significance or specificity. DE results may filter or highlight pairs whose ligand or receptor changes; mark those visibly. Cap the pairs shown by a stated rule.
- **Report:** inference method and database, expression filter (minimum fraction of expressing cells), significance threshold, selection rule.

### 12. Cell-cell communication between tissues or sources
- **Needs:** two or more tissues or sources where cross-source signalling is biologically plausible.
- **Draw:** as figure 11, with senders and receivers labelled by source. Restrict to secreted or circulating ligands, since contact-dependent pairs cannot act between separate tissues. Use DE results to filter or highlight pairs.
- **Report:** how sources were paired, how cell populations were matched, expression filters, method. Label the result as inferred potential communication.

## Legend checklist

Write one legend entry per figure in `figure_legends.md`. Each states, where applicable:
- what each panel shows and what a point, bar or dot represents (cell, sample, replicate)
- n: cells, samples and biological replicates per group
- statistical test, unit of analysis, multiple-testing correction, thresholds
- error-bar and interval definitions
- every colour and size encoding, including scaling and clipping limits
- software, method and reference-resource versions

Anything too long for the legend goes into a Methods note delivered alongside.
