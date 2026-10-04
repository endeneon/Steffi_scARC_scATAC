# Handoff: GWAS-annotated Gviz SNP panels

**Date:** 2026-10-03
**Primary script:** `Steffi_works/generate_gviz_multiple_panels.R`
**Reference (adapted from):** `Steffi_works/gviz_plot_5_macrophage_candidates.R`
**Submission wrapper:** `Steffi_works/bsub_plot_gviz_general.sh`

## Goal

For each significant ASoC SNP in macrophages, render a multi-track Gviz panel
(per-cell-type ATAC coverage + gene model + optional HCC GWAS track) and tile
them into a multipage landscape PDF under `macrophage_SNP_by_cell_type/`.

## What the script does (end to end)

1. **Init block** — loads libraries (unused ones commented out), sets an
   LSF/IDE-aware `future` plan (`multicore` under Rscript, `multisession` in
   IDE), configures ArchR threads/genome, and registers `doFuture`.
2. **Inputs**
   - ArchR project: `ArchR_merged_ATAC_multiome_obj_multiple_cell_types`
     (grouping column `projected_barcodes_multiple_cell_types`, 6 cell types:
     B/Plasma cell, Endothelial, Fibroblast/HSC, Hepatocyte, Macrophage, T cell;
     323,803 cells, no NAs).
   - SNP table: `sig_ASoC_by_celltype/sig_ASoC_in_Macrophage_annotated_nopromoter.tsv`
     read with `quote = ""` (annotation fields contain apostrophes → otherwise
     "EOF within quoted string" drops rows). Subset to `df_2_plot`
     (seqnames, start, variantID, SYMBOL).
   - GWAS VCFs (GRCh38, chr-prefixed, bgzipped + tabix-indexed):
     - EAS: `HCC_GWAS/GRCh38/hcc_ea_011123.vcf.gz`
     - EUR: `HCC_GWAS/GRCh38/hcc_eur_200324.vcf.gz`
     - PVAL is in the INFO field.
3. **`plot_snp_tracks(snp_chr, snp_pos, gviz_window, gviz_bin_size, main)`** —
   builds, per SNP: ideogram, axis, optional GWAS dot track, one RPGC-normalized
   coverage `DataTrack` per cell type, optional peak track, squished gene track;
   wraps data tracks in a `HighlightTrack` at the SNP; calls `Gviz::plotTracks()`.
4. **Batch build** — `foreach %dopar%` over all rows of `df_2_plot`, each panel
   captured as a grid grob via `grid::grid.grabExpr(..., wrap.grobs = TRUE)`
   (critical: without `wrap.grobs = TRUE`, Gviz reuses grob names across tracks
   and all but one coverage track collapse). `.errorhandling = "pass"` + a
   post-loop `grid::is.grob` filter drop failed SNPs (warn + message which).
5. **PDF output** — 2×2 panels per page on a 14×10.5in landscape page, 5pt
   per-panel subtitle (`SYMBOL, variantID (chr:pos)`), per-page "Page N of M"
   footer, written to `macrophage_SNP_by_cell_type/macrophage_SNP_by_cell_type_panels.pdf`.

## GWAS track (`use_GWAS_track`)

- Toggle at line ~191: `use_GWAS_track <- TRUE`.
- Helper `.read_gwas_region(vcf_file, chr, from, to, ancestry)` uses
  `Rsamtools::scanTabix` to read ONLY the window (VCFs are 220–630 MB), writes
  header+lines to a tempfile, parses with `vcfR::read.vcfR`, extracts INFO/PVAL,
  returns CHROM/POS/`-log10(PVAL)`/ancestry (finite rows only).
- Track: a `Gviz::OverlayTrack` placed at the TOP of the track list, wrapping:
  1. a points `DataTrack`, two groups (EAS, EUR), `type = "p"`,
     `col = c("darkred","darkblue")`, `alpha = 0.8`, `legend = FALSE`,
     with scalar `baseline = 0` (solid grey70 lower-margin line); and
  2. a flat `type = "l"` significance-line `DataTrack` at y = -log10(5e-8) ≈ 7.3
     (dashed grey50), spanning `[from, to]`.
  Both share the same `ylim` so the overlay lines up.
- **Y-axis:** fixed to `c(0, 5)` unless any point in the window exceeds 5, then
  autoscale to `c(0, max)`.
- **Why an OverlayTrack:** a vector `baseline = c(0, -log10(5e-8))` on a single
  DataTrack ERRORS in Gviz (R 4.3: `'length = 2' in coercion to 'logical(1)'`),
  which silently dropped all panels in batch job 325705866. `baseline` must be
  scalar; the second line is therefore a separate overlaid track.
- When `use_GWAS_track == TRUE`, batch `gviz_window` widens 5000 → 50000.

## Environment / dependencies

- vcfR was NOT installed; installed via mamba into env `jwen_scRNA_singCellaR`
  (`mamba install -n jwen_scRNA_singCellaR -c conda-forge -c bioconda r-vcfr`,
  r-vcfr 1.15.0). Rsamtools + VariantAnnotation already present.
- foreach worker `.packages` include: Gviz, grid, gridExtra, ArchR,
  GenomicRanges, Rsamtools, vcfR.
- R 4.3.0 session; ArchR project path is relative to `setwd("Steffi_works")`.

## Status

- Script parses cleanly after the latest edits.
- Verified by single-panel renders: track order (GWAS on top → 6 coverage →
  genes), legend removal, and the 50 kb window.
- **GWAS y-axis + baselines VERIFIED (2026-10-03):**
  - Low-signal window (rs186811887, chr18:47749304, SMAD2): y-axis fixed at
    0–5, y=0 line drawn, 5e-8 line correctly off-screen (above 5).
  - Genome-wide-significant window (rs8100204, chr19:19282905, EUR
    PVAL=1.74e-09): autoscale branch triggered (y-axis ~0–50) and the dashed
    5e-8 line drawn across the point cloud. Lead variant found by awk-scanning
    the VCFs for the smallest PVAL.

## Next steps

1. Run the full parallel batch (submit via `bsub_plot_gviz_general.sh`); note
   ±50 kb @ 50 bp bins = ~2000 bins/track, ~10× more fragment reading than the
   5 kb window, so budget more time/cores.
2. Inspect page 1 of the regenerated PDF (expect all 30 SNPs to render now that
   the vector-baseline bug is fixed).

```json
{
  "project": "gviz_gwas_snp_panels",
  "primary_script": "Steffi_works/generate_gviz_multiple_panels.R",
  "reference_script": "Steffi_works/gviz_plot_5_macrophage_candidates.R",
  "submit_script": "Steffi_works/bsub_plot_gviz_general.sh",
  "handoff_file": "Steffi_works/HANDOFF_gviz_gwas_panels.md",
  "working_dir": "Steffi_works (script does setwd)",
  "r_version": "4.3.0",
  "conda_env": "jwen_scRNA_singCellaR",
  "inputs": {
    "archr_project": "ArchR_merged_ATAC_multiome_obj_multiple_cell_types",
    "group_col": "projected_barcodes_multiple_cell_types",
    "cell_types": ["B/Plasma cell", "Endothelial", "Fibroblast/HSC", "Hepatocyte", "Macrophage", "T cell"],
    "snp_table": "sig_ASoC_by_celltype/sig_ASoC_in_Macrophage_annotated_nopromoter.tsv",
    "snp_table_read_opts": {"quote": ""},
    "gwas_vcf_EAS": "HCC_GWAS/GRCh38/hcc_ea_011123.vcf.gz",
    "gwas_vcf_EUR": "HCC_GWAS/GRCh38/hcc_eur_200324.vcf.gz",
    "gwas_pval_field": "INFO/PVAL",
    "genome": "hg38_chr_prefixed"
  },
  "key_function": "plot_snp_tracks(snp_chr, snp_pos, gviz_window, gviz_bin_size, main)",
  "gwas_helper": ".read_gwas_region(vcf_file, chr, from, to, ancestry)",
  "gwas_track": {
    "toggle_var": "use_GWAS_track",
    "toggle_value": true,
    "position": "top of trackList",
    "structure": "Gviz::OverlayTrack(points DataTrack + flat significance-line DataTrack), shared ylim",
    "type": "p",
    "colors": {"EAS": "darkred", "EUR": "darkblue"},
    "alpha": 0.8,
    "legend": false,
    "yaxis": "fixed c(0,5); autoscale to c(0,max) only if any point > 5",
    "y0_line": {"y": 0, "col": "grey70", "lty": "solid", "role": "lower margin", "impl": "scalar baseline on points track"},
    "sig_line": {"y": "-log10(5e-8)", "col": "grey50", "lty": "dashed", "role": "genome-wide significance", "impl": "separate flat type='l' DataTrack overlaid"},
    "window_when_on": 50000,
    "window_when_off": 5000
  },
  "pipeline": {
    "parallel": "foreach %dopar% over doFuture",
    "worker_packages": ["Gviz", "grid", "gridExtra", "ArchR", "GenomicRanges", "Rsamtools", "vcfR"],
    "grob_capture": "grid::grid.grabExpr(..., wrap.grobs = TRUE)",
    "error_handling": "pass + grid::is.grob filter drops failed SNPs",
    "layout": "2x2 per page, 14x10.5in landscape",
    "subtitle": "5pt 'SYMBOL, variantID (chr:pos)'",
    "footer": "Page N of M",
    "output_pdf": "macrophage_SNP_by_cell_type/macrophage_SNP_by_cell_type_panels.pdf"
  },
  "critical_gotchas": [
    "wrap.grobs=TRUE is REQUIRED or Gviz tracks collapse to one when grabbed",
    "read.table needs quote='' for the annotation TSV",
    "each worker re-derives its own Arrow handles (HDF5 handles not shareable)",
    "+/-50kb @ 50bp bins is ~10x slower fragment reading than the 5kb window",
    "Gviz baseline must be SCALAR; vector baseline errors ('length = 2' in coercion to 'logical(1)') and dropped all panels in job 325705866 -> use an OverlayTrack for the second line"
  ],
  "status": "GWAS track fully verified on single panels (low-signal fixed 0-5 + high-signal autoscale/5e-8); OverlayTrack fix applied; script parses; batch not yet rerun",
  "verified": {
    "low_signal_window": {"variant": "rs186811887", "locus": "chr18:47749304", "gene": "SMAD2", "result": "y 0-5 fixed, y=0 line, 5e-8 off-screen"},
    "gwsig_window": {"variant": "rs8100204", "locus": "chr19:19282905", "ancestry": "EUR", "pval": 1.74e-09, "result": "autoscale y~0-50, dashed 5e-8 line drawn"}
  },
  "next_steps": [
    "run full parallel batch via bsub wrapper (budget extra time for 50kb windows)",
    "inspect page 1 of regenerated PDF; expect all 30 SNPs to render"
  ]
}
```
