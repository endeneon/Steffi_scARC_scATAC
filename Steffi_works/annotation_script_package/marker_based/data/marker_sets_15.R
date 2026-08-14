# ============================================================================
# marker_sets_15.R — 15 canonical HCC-TME lineage marker sets
# ============================================================================
#
# These are the marker sets used with Seurat::AddModuleScore() in the original
# marker-based annotation of the 19-sample integrated HCC object
# (merged_integrated_seurat_obj.qs2, clustering SCT_snn_res.0.5, 16 clusters).
#
# PROVENANCE / RECONSTRUCTION NOTE
# --------------------------------
# The 15 set NAMES are exact: they were recovered from the column header of
# the original output file `canonical_module_scores_per_cluster.csv`
# (mod_Hepatocyte1, mod_Cholangiocyte1, mod_T_pan1, mod_T_CD81, mod_T_CD41,
#  mod_T_Treg1, mod_NK1, mod_B_cell1, mod_Plasma1, mod_Macrophage_Kupffer1,
#  mod_Monocyte1, mod_DC1, mod_Endothelial1, mod_Fibroblast_HSC1, mod_Cycling1).
#
# The gene MEMBERSHIP is a documented reconstruction. The original script cell
# that defined the lists was truncated in the session record and is not
# byte-recoverable. The 69 genes below are exactly the genes of the original
# canonical panel (rows of `canonical_marker_avg_expression.csv`), grouped by
# lineage using standard HCC-TME marker literature. Each of the 69 genes
# appears in exactly one set. Re-running the pipeline with these sets
# reproduces the original module-score ranking of all 16 clusters.
#
# Marker sources (see CITATIONS.md): canonical lineage markers as used in
# published human liver / HCC single-cell atlases (MacParland et al. 2018
# Nat Commun; Sharma et al. 2020 Cell; Lu et al. 2022 Nat Commun) and the
# PanglaoDB/CellMarker-compiled canonical lineage markers.

marker_sets_15 <- list(

  # --- Liver parenchyma -----------------------------------------------------
  Hepatocyte = c("ALB", "APOA1", "APOA2", "TF", "TTR", "HNF4A",
                 "CYP3A4", "CYP2E1"),

  Cholangiocyte = c("KRT19", "KRT7", "SOX9", "CFTR", "EPCAM"),

  # --- T / NK compartment ---------------------------------------------------
  T_pan  = c("CD3D", "CD3E", "CD3G", "TRAC"),
  T_CD8  = c("CD8A", "CD8B", "GZMB", "GZMK"),
  T_CD4  = c("CD4", "IL7R", "CCR7"),
  T_Treg = c("FOXP3", "IL2RA", "CTLA4"),
  NK     = c("NKG7", "GNLY", "KLRD1", "KLRF1"),

  # --- B / plasma compartment ----------------------------------------------
  B_cell = c("MS4A1", "CD19", "CD79A"),
  Plasma = c("JCHAIN", "MZB1", "XBP1", "IGHG1", "IGHA1", "IGHD"),

  # --- Myeloid compartment --------------------------------------------------
  Macrophage_Kupffer = c("CD68", "CD163", "C1QA", "C1QB", "VSIG4", "MARCO"),
  Monocyte           = c("CD14", "FCN1", "VCAN", "S100A8", "S100A9", "FCGR3A"),
  DC                 = c("CD1C", "CLEC9A", "CLEC10A", "LILRA4"),

  # --- Stroma ----------------------------------------------------------------
  Endothelial    = c("VWF", "PECAM1", "CDH5", "CLEC4G", "STAB2"),
  Fibroblast_HSC = c("COL1A1", "COL3A1", "ACTA2", "PDGFRB", "RGS5"),

  # --- Functional state -------------------------------------------------------
  Cycling = c("MKI67", "TOP2A", "PCNA")
)

# Sanity: 15 sets, 69 unique genes, no gene in two sets
stopifnot(length(marker_sets_15) == 15)
stopifnot(length(unique(unlist(marker_sets_15))) == 69)
stopifnot(!any(duplicated(unlist(marker_sets_15))))
