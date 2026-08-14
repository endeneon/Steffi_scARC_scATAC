# Reference-based annotation (Package B)

An **alternative, new** annotation route that was **not** part of the original
analysis. It annotates the query object by mapping cells/clusters onto a
published HCC single-cell atlas with [SingleR](https://github.com/SingleR-inc/SingleR).

## Reference dataset

**Lu et al. 2022, *Nature Communications*** — "A single-cell atlas of the
multicellular ecosystem of primary and metastatic hepatocellular carcinoma"

- **GEO accession**: [GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614)
- **Size**: 71,915 single cells; 10 HCC patients; 4 tissue sites
  (non-tumor liver, primary tumor, portal vein tumor thrombus, metastatic
  lymph node)
- **Files used**:
  - `GSE149614_HCC.metadata.updated.txt.gz` (478 KB) — cell annotations
    (columns: `Cell, sample, res.3, site, patient, stage, virus, celltype`)
  - `GSE149614_HCC.scRNAseq.S71915.count.txt.gz` (158 MB) — genes × cells
    count matrix
- **Reference cell-type labels** (`celltype` column): T/NK (25,591),
  Hepatocyte (20,782), Myeloid (15,947), B (3,685), Endothelial (3,644),
  Fibroblast (2,266)

## Scripts

| Script | Purpose |
|---|---|
| `01_fetch_reference.R` | Download both files from GEO; build a log-normalized SingleR reference (optionally aggregated to per-celltype pseudobulk profiles for speed); save as `.rds` |
| `02_singler_annotate.R` | Run SingleR on the query object at cell level and cluster level; adds `singler_label_cell` / `singler_label_cluster` |
| `03_compare_methods.R` | Confusion matrix + per-cluster agreement between the marker-based labels (Package A) and SingleR labels, using a harmonization map between the two taxonomies |

## Quick start

```bash
Rscript 01_fetch_reference.R reference_data
Rscript 02_singler_annotate.R ../annotated_seurat.qs2 \
    reference_data/gse149614_singler_ref.rds annotated_singler.qs2
Rscript 03_compare_methods.R annotated_singler.qs2 marker_vs_singler
```

## Notes

- The Lu et al. taxonomy is **coarser** than the marker-based labels
  (6 lineages vs 15 fine labels) and has no "Cycling" state; the comparison
  script harmonizes labels before scoring agreement.
- `aggregateReference()` (default) collapses the 71,915-cell reference to 6
  pseudobulk profiles — much faster and appropriate for lineage-level
  mapping. Set `aggregate = FALSE` for full single-cell reference scoring.
- Dependencies: `SingleR`, `SingleCellExperiment`, `Seurat`, `qs2`, `Matrix`.
