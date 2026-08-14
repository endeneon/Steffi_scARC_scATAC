# HCC scRNA-seq Annotation Script Package

Scripts used to annotate the integrated 19-sample HCC object
(`merged_integrated_seurat_obj.qs2`; 223,521 cells → 165,115 after QC;
Harmony-integrated; clustering `SCT_snn_res.0.5`, 16 clusters), plus a new
reference-based alternative.

## Important: which method was actually used?

The **original annotation was marker-based (manual), NOT reference-based.**
No reference dataset was fetched or used in the original analysis. Cell types
were assigned by cross-referencing (a) per-cluster Wilcoxon marker DE
(`presto`), (b) `AddModuleScore` on 15 canonical HCC-TME lineage marker sets,
and (c) canonical single-gene expression — see `annotation_dictionary.md`
(delivered with the original analysis) for the full per-cluster reasoning.

This package therefore contains two pipelines:

| Package | Folder | Status | Reference dataset? |
|---|---|---|---|
| **A. Marker-based** | `marker_based/` | Faithful reproduction of the original analysis | No — markers from literature/domain knowledge |
| **B. Reference-based** | `reference_based/` | **New** alternative (SingleR) | Yes — Lu et al. 2022, GSE149614 |

## Package A — marker-based pipeline (`marker_based/`)

Run order (each script also works via `source()` + function call):

```bash
Rscript 01_load_qc.R merged_integrated_seurat_obj.qs2 step1.qs2 25
Rscript 02_marker_de.R step1.qs2 cluster_markers_res0.5
Rscript 03_module_scores.R step1.qs2 step3.qs2
Rscript 04_annotate_clusters.R step3.qs2 step4.qs2
Rscript 05_malignant_classifier.R step4.qs2 step5.qs2
Rscript 06_macrophage_subclustering.R step5.qs2 annotated_seurat.qs2
Rscript 07_tests.R annotated_seurat.qs2
```

| Script | What it does |
|---|---|
| `01_load_qc.R` | Load `.qs2`; `percent.mt ≤ 25` filter (223,521 → 165,115 cells) |
| `02_marker_de.R` | `JoinLayers` on RNA assay; `presto::wilcoxauc` per cluster |
| `03_module_scores.R` | `AddModuleScore` for the 15 lineage sets (`data/marker_sets_15.R`); per-cluster score + 69-gene expression summaries |
| `04_annotate_clusters.R` | Apply the 16-cluster label map (`data/cluster_label_map.csv`) → `celltype_fine` / `celltype_broad` |
| `05_malignant_classifier.R` | 9-gene malignant + 8-gene mature hepatocyte module scores; quartile classifier → `malignant_status` |
| `06_macrophage_subclustering.R` | Re-integrate 25,516 macrophages (Harmony per sample), cluster at res 0.3 → 6 subtypes incl. 4 doublet classes → `macrophage_subtype` |
| `07_tests.R` | Sanity suite: counts, label coverage, module-score peaks, classifier split, macrophage-subtype marker checks |

### Reuse on a new dataset

- Steps 01–03 are dataset-agnostic (only the mt-gene pattern `^MT-` assumes
  human data).
- Step 04 requires a label map matching **your** clustering: the shipped
  `cluster_label_map.csv` is specific to the original 16 clusters. For a new
  object, inspect the outputs of steps 02–03 (top markers + module scores per
  cluster) and write a new CSV with the same columns
  (`cluster, fine_label, broad_label, n_cells, key_evidence`).
- Steps 05–06 assume `celltype_broad` contains `"Hepatocyte"` and
  `"Macrophage"` labels.

## Package B — reference-based pipeline (`reference_based/`)

See `reference_based/README.md`. Downloads the Lu et al. 2022 HCC atlas
(GSE149614) from GEO, builds a SingleR reference, annotates the query, and
compares against the marker-based labels.

## Dependencies

R packages: `Seurat` (v5), `qs2`, `dplyr`, `Matrix`, `presto`, `harmony`,
`SingleR`, `SingleCellExperiment`.

## Provenance & limitations

- The 15 module-set **names** are exact (recovered from the original output
  `canonical_module_scores_per_cluster.csv`). The **gene membership** is a
  documented reconstruction from the original 69-gene canonical panel
  (`canonical_marker_avg_expression.csv`); the original defining code cell was
  truncated in the session record. See `marker_based/data/marker_sets_15.R`.
- The malignant classifier is marker-based and **not CNV-validated**; an
  inferCNV/CopyKAT confirmation run is recommended.
- Macrophage doublets are flagged by marker co-expression, not by a formal
  doublet detector.
- Full biological interpretation and per-cluster evidence:
  `annotation_dictionary.md` (delivered with the original analysis).
