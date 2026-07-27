# Chat Session Summary

## Session Metadata

- **Date:** 2026-07-26
- **Workspace:** `Multiome_main` (repo `endeneon/Steffi_scARC_scATAC`, branch `main`)
- **Primary input file:** `Steffi_works/sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv`
- **Primary deliverable:** `Steffi_works/analysis/ASoC_hepatocyte_liver_cancer_association_summary.md`
- **Goal of the session:** Interrogate an allele-specific open chromatin (ASoC) result table from human hepatocytes to determine which variants (`variantID`, dbSNP rsIDs) and which genes (`ENSEMBL`, `SYMBOL`) are associated with hepatocellular carcinoma (HCC), liver-related disease/traits more broadly, or gastric cancer; then narrow to genes whose SNPs strongly disrupt TF motifs and attach **measured** TCGA-LIHC tumor-vs-normal fold-changes. Progressively expand the evidence base from ClinVar to literature, GWAS Catalog, formal LD-proxy testing, and a scripted TCGA-LIHC differential-expression analysis, recording findings in a reusable markdown report.

## Data Context

- The TSV is an **ASoC analysis in hepatocytes**. Each row is a significant heterozygous SNP showing allelic imbalance in chromatin accessibility, annotated with genomic coordinates (GRCh38), overlapping/nearest gene (`geneId`, `ENSEMBL`, `SYMBOL`, `GENENAME`), peak coordinates, `annotation` (mostly `Promoter (<=1kb)`), TF motif matches (`motif_names`), and motifbreakR disruption columns (`mb_disrupted_TFs`, etc.).
- Coordinates are **GRCh38**. Variants are common regulatory SNPs (promoter/intronic/intergenic), not coding pathogenic alleles.

## Tasks Completed

1. **Initial read + first-pass annotation of the ASoC table.**
   - Read the full table and identified candidate HCC/liver-relevant genes from the `SYMBOL`/`ENSEMBL` columns.
   - Verified two representative SNPs against dbSNP (`rs2980223` near HNF4G; `rs2060982` near CPT1A): both "Not Reported in ClinVar."
   - **Outcome:** Established that variant-level disease relevance is essentially absent in ClinVar; relevance is gene-level. Created the summary markdown file (Section 1 + Section 2).

2. **Expanded the variant search online (literature + GWAS Catalog + LD-proxy reasoning).**
   - Europe PMC text-mining by exact rsID: no direct variant→HCC associations; only stray/false matches.
   - GWAS Catalog by rsID: the 6p21 MHC ASoC SNP `rs16899941` (GPANK1/BAG6) is genome-wide significant (p ≈ 3×10⁻⁵²) but for **hypothyroidism / rheumatoid arthritis** (autoimmune), not HCC.
   - Cross-referenced ASoC coordinates against known HCC GWAS loci; only the 6p21 MHC cluster sits near a documented HCC locus (MICA, `rs2596542`, ~0.13 Mb away).
   - Strengthened gene-level evidence with Europe PMC co-mention counts.
   - **Outcome:** Added **Section 1b** (expanded online search), added Europe PMC citation counts to gene tables, added `TOB1` (gastric-cancer SNP literature at `rs9898809`).

3. **Answered whether Ensembl LD REST API needs a token; ran formal LD-proxy test (EUR).**
   - Confirmed the **Ensembl LD REST API is public / no token required** — only a 1000 Genomes population is needed (LDlink is the one requiring a token).
   - Ran pairwise LD in `1000GENOMES:phase_3:EUR` between the four 6p21 MHC ASoC SNPs (`rs16899941`, `rs2517765`, `rs1632906`, `rs78905038`) and the HCC MICA lead `rs2596542`. All returned empty arrays → **not LD proxies**.
   - **Outcome:** Added **Section 1b(iv)** documenting the rejected proxy hypothesis; updated caveat/summary/footer.

4. **Repeated the LD-proxy test in EAS.**
   - Same four SNPs vs `rs2596542` in `1000GENOMES:phase_3:EAS` → all empty → not LD proxies in EAS either.
   - **Outcome:** Updated Section 1b(iv) table to include an EAS r² column; conclusion holds in both EUR and EAS.

5. **Broadened the search to any liver-related trait and to gastric cancer.**
   - Europe PMC grouped queries across ~140 rsIDs against liver terms (NAFLD, NASH, steatosis, cirrhosis, hepatitis, bilirubin, aminotransferase, cholestasis, "gamma glutamyl") and against gastric-cancer terms.
   - GWAS Catalog structured EFO traits for the common metabolic SNPs.
   - **Findings:** Exactly two variant-level signals emerged — `rs2060982` (CPT1A → hepatic lipid/fatty-acid metabolism) and `rs9898809` (TOB1 → gastric cancer). HNF4G `rs2980223` and KLF10 `rs2471847` have no GWAS Catalog associations.
   - **Outcome:** Added **Section 1c** (broadened search results) and updated the summary/footer.

6. **Expanded the investigation to gene-level HCC-vs-normal expression (literature first).**
   - Ran Europe PMC direction-informative queries for the top protein-coding ASoC genes (FZD7, USP22, NDRG1, TP53INP1, MTSS1, CPT1A, CES1, KLF10, NEK6, TRIB1, etc.).
   - **Outcome:** Added **Section 3** with an initial literature-based up/down/context classification, Europe PMC co-mention counts, and TCGA-LIHC verification-resource URLs (GEPIA2/UALCAN/TNMplot/HPA).

7. **Narrowed the gene panel to confirmed strong TF-motif disruption (`mb_n_strong` > 0).**
   - Filtered Section 3 to genes whose ASoC SNP(s) have `mb_n_strong` > 0 (max across a gene's SNPs). Removed four genes whose every SNP had `mb_n_strong = 0`: **TRIB1** (rs17663005), **TOB1/TOB1-AS1** (rs9898809), **ING1** (rs1441042), **TLR5** (rs2302597).
   - Kept NDRG1 (via rs58065091 = 15) and MAP3K2 (via rs13034125 = 3) despite their promoter SNPs scoring 0.
   - **Outcome:** Added an `mb_n_strong` column (with disrupted TFs) and a filter note.

8. **Replaced literature directions with MEASURED TCGA-LIHC fold-changes (scripted edgeR DE).**
   - Used the `r_45_python_312` conda env (R 4.5.3). Downloaded the UCSC Xena GDC hub `TCGA-LIHC.star_counts.tsv.gz` (log2(count+1) → raw counts), classified samples by TCGA barcode into **374 tumor (code 01)** vs **50 solid-normal (code 11)**, and ran **edgeR** (`filterByExpr` → TMM → `estimateDisp` → QL F-test) with thresholds **|log2FC| ≥ 0.2, FDR < 0.05**. No `TCGAbiolinks` install was needed (edgeR/limma/data.table/SummarizedExperiment already present).
   - Genome-wide ASoC intersection (230 genes with Ensembl IDs): **111 Up, 34 Down, 57 ns, 28 not-tested**.
   - Rewrote **Section 3** to report measured log2FC + FDR, organized by measured direction, noting concordance/discordance vs literature.
   - **Outcome:** Two literature-"Up" genes are measured **Down** in bulk TCGA-LIHC — **CPT1A** (−0.48) and **NEK6** (−0.40); literature tumor-suppressors **TP53INP1/MTSS1** were **not significant**. Confirmed Up: CHEK1 (+2.22), HDAC11 (+1.76), LAMC1 (+1.39), NDRG1 (+1.19), CAPN2, GSDME, YWHAZ, USP22, FZD7, TPD52. Confirmed Down: KLF10 (−1.44), DUSP10, HES1, MAP3K2, LRP5.

## Key Decisions & Rationale

- **Report disease relevance at the gene level, not the variant level.** ASoC SNPs are regulatory/eQTL-type variants; ClinVar/GWAS rarely annotate them directly. The biologically meaningful readout is the regulated gene.
- **Use Europe PMC + GWAS Catalog REST APIs over ad-hoc web pages** for reproducible, machine-readable evidence (hit counts, EFO traits) and to text-mine exact rsIDs.
- **Test the one plausible proxy hypothesis formally rather than asserting proximity.** The 6p21 MHC cluster is ~0.13 Mb from the HCC MICA lead, but genomic proximity in the gene-dense MHC does not imply LD. Pairwise LD in EUR and EAS rejected the proxy.
- **Ensembl LD REST API chosen over LDlink** because it is public/tokenless; population parameter (`EUR`, then `EAS`) supplied explicitly. Flagged that MHC LD is ancestry-specific and should be re-run in the donors' matched population if non-EUR/EAS.
- **Distinguish quantitative-trait GWAS signals from disease diagnoses.** `rs2060982`/CPT1A associates with fatty-acid/lipid measures (relevant to NAFLD physiology) but is not a liver-disease diagnosis.
- **Prefer measured over literature-inferred directions.** Section 3 directions were replaced with a scripted TCGA-LIHC edgeR analysis; discordances (CPT1A, NEK6) are surfaced explicitly rather than hidden.
- **Minimize new installs.** Reused base Bioconductor packages (edgeR/limma/data.table) and the UCSC Xena GDC pre-built matrix instead of installing/using `TCGAbiolinks` + full GDC download.
- **Restrict the expression panel to functionally credible variants.** Applied a `mb_n_strong` > 0 (motifbreakR strong TF disruption) filter so the expression analysis focuses on SNPs likely to alter TF binding.

## Code Changes

This was an analysis/research session. The concrete artifact is a markdown report plus a set of reproducible API queries.

### File: `Steffi_works/analysis/ASoC_hepatocyte_liver_cancer_association_summary.md` (created, then iteratively updated)

Structure of the report:
- **Section 1** — Variant-level (`variantID`) findings: no curated variant→HCC associations in dbSNP/ClinVar.
- **Section 1b** — Expanded online search: (i) Europe PMC literature by rsID; (ii) GWAS Catalog by rsID (MHC autoimmune signal); (iii) proximity to known HCC GWAS loci; (iv) formal Ensembl LD-proxy test (EUR + EAS), result: rejected.
- **Section 1c** — Broadened search: liver-related (any) and gastric-cancer associations. Two positive variant-level hits.
- **Section 2** — Gene-level (`ENSEMBL`/`SYMBOL`) findings with Europe PMC HCC co-mention counts.
- **Section 3** — Gene-level differential expression: strong-TF-disruption panel (`mb_n_strong` > 0) annotated with **measured** TCGA-LIHC tumor-vs-normal log2FC + FDR, organized by measured direction; data-source URLs; reproduce instructions.
- **Summary + Recommended follow-up + provenance footer.**

### New files: `Steffi_works/analysis/tcga_lihc_de/` (scripted TCGA-LIHC DE)

- `run_tcga_lihc_de.R` — downloads Xena GDC matrix + probemap, runs edgeR, writes full and ASoC-intersected DE tables.
- `intersect_asoc.R` — re-intersects the saved DE table with ASoC Ensembl IDs (fixes a merge key case mismatch) and prints the Section-3 panel.
- `TCGA-LIHC_tumor_vs_normal_DE_full.tsv` — genome-wide DE (22,759 tested genes).
- `ASoC_genes_TCGA-LIHC_DE.tsv` — 230 ASoC genes with measured log2FC/FDR/direction.
- `data/` — cached Xena matrix + probemap (~43 MB; regenerable).

Key gene-level findings table (strong/direct liver-cancer relevance):

| SYMBOL   | ENSEMBL         | Variant(s)               | Lit. (HCC) | Note                                        |
| -------- | --------------- | ------------------------ | ---------- | ------------------------------------------- |
| HNF4G    | ENSG00000164749 | rs2980223                | ~187       | Hepatocyte-identity TF                      |
| NDRG1    | ENSG00000104419 | rs1017105790, rs58065091 | ~1123      | HCC progression/metastasis/prognosis        |
| TP53INP1 | ENSG00000164938 | rs10097617               | ~653       | p53 stress-response tumor suppressor        |
| CPT1A    | ENSG00000110090 | rs2060982                | metabolic  | Fatty-acid β-oxidation; lipid reprogramming |
| INSIG1   | ENSG00000186480 | rs576875943              | metabolic  | SREBP lipogenesis; NAFLD→HCC                |
| KLF10    | ENSG00000155090 | rs2471847                | moderate   | TGF-β-inducible tumor-suppressor-type TF    |

Secondary candidates: FZD7 (~755), MTSS1 (~311), USP22 (~699), TOB1 (gastric), CES1, DUSP10, LAMC1, TLR5, CAPN2, GSDME.

Two variant-level signals from the broadened search:

| variantID | Gene            | Source       | Trait                                                                | Category                         |
| --------- | --------------- | ------------ | -------------------------------------------------------------------- | -------------------------------- |
| rs2060982 | CPT1A           | GWAS Catalog | saturated FA / total FA % (p=3×10⁻³⁰); IDL triglycerides (p=2×10⁻¹³) | Liver-related (lipid metabolism) |
| rs9898809 | TOB1 / TOB1-AS1 | Europe PMC   | gastric cancer risk (Chinese Han)                                    | Gastric cancer                   |

### Reproducible API queries used

Ensembl LD pairwise (public, no token; substitute `EAS` for the East Asian panel):

```bash
# Pairwise r2 between a 6p21 MHC ASoC SNP and the HCC MICA lead rs2596542
curl -s "https://rest.ensembl.org/ld/human/pairwise/rs16899941/rs2596542?population_name=1000GENOMES:phase_3:EUR;content-type=application/json"
# -> [] (empty array = not in appreciable LD)

# Proxies for a SNP at r2 >= 0.8 in EUR
curl -s "https://rest.ensembl.org/ld/human/rs16899941/1000GENOMES:phase_3:EUR?r2=0.8;content-type=application/json"
```

GWAS Catalog by rsID (structured EFO traits at the end of the JSON):

```bash
curl -s "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs2060982/associations?projection=associationBySnp"
# -> efoTraits include "triglycerides in IDL measurement" (p=2e-13) and
#    "saturated fatty acids to total fatty acids percentage" (p=3e-30)

curl -s "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/rs16899941/associations"
# -> efoTraits: "hypothyroidism", "rheumatoid arthritis" (p=3e-52), 6p21.33 MHC
```

Europe PMC literature by exact rsID (compact hit counts + titles):

```bash
# Gastric cancer across rsIDs -> only TOB1 (rs9898809) hit
curl -s 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=(rs9898809)%20AND%20(%22gastric%20cancer%22)&format=json&pageSize=5&resultType=lite'

# Gene x HCC co-mention counts (hitCount field)
curl -s 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%22NDRG1%22%20AND%20%22hepatocellular%20carcinoma%22&format=json&pageSize=1&resultType=lite'
```

### Scripted TCGA-LIHC differential expression (edgeR; no TCGAbiolinks)

Run in the `r_45_python_312` conda env (R 4.5.3). Core of `run_tcga_lihc_de.R`:

```r
suppressMessages({ library(data.table); library(edgeR) })

# 1) Xena GDC hub TCGA-LIHC STAR counts (log2(count+1)) + probemap
counts_url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_counts.tsv.gz"
download.file(counts_url, "TCGA-LIHC.star_counts.tsv.gz", mode = "wb")

# 2) ASoC Ensembl IDs (strip version)
asoc <- fread("sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv")
asoc[, ENSEMBL_base := sub("\\..*$", "", ENSEMBL)]

# 3) counts matrix; classify tumor(01)/normal(11) by TCGA barcode chars 14-15
mat <- fread(cmd = "zcat TCGA-LIHC.star_counts.tsv.gz")
setnames(mat, 1, "Ensembl"); mat[, Ensembl_base := sub("\\..*$", "", Ensembl)]
samples <- setdiff(colnames(mat), c("Ensembl", "Ensembl_base"))
stype   <- substr(samples, 14, 15)
grp     <- ifelse(stype == "11", "normal", ifelse(grepl("^0", stype), "tumor", NA))
keep_s  <- samples[!is.na(grp)]; grp <- factor(grp[!is.na(grp)], c("normal", "tumor"))

# 4) reverse log2(x+1) -> raw counts, edgeR QL F-test tumor vs normal
cnt <- round(2^as.matrix(mat[, ..keep_s]) - 1); cnt[cnt < 0] <- 0
rownames(cnt) <- mat$Ensembl_base; cnt <- cnt[!duplicated(rownames(cnt)), ]
y <- DGEList(cnt, group = grp); y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y); design <- model.matrix(~ grp); y <- estimateDisp(y, design)
fit <- glmQLFit(y, design); qlf <- glmQLFTest(fit, coef = 2)
res <- as.data.table(topTags(qlf, n = Inf)$table, keep.rownames = "ENSEMBL_base")

# 5) intersect + classify (|log2FC|>=0.2 & FDR<0.05)
m <- merge(unique(asoc[, .(ENSEMBL_base, SYMBOL)]), res, by = "ENSEMBL_base", all.x = TRUE)
m[, direction := fifelse(FDR < 0.05 & logFC >=  0.2, "Up",
                 fifelse(FDR < 0.05 & logFC <= -0.2, "Down", "ns"))]
```

Measured result for the strong-TF-disruption panel (log2FC, TCGA-LIHC tumor vs normal): CHEK1 +2.22, HDAC11 +1.76, LAMC1 +1.39, NDRG1 +1.19, CAPN2 +0.76, GSDME +0.70, YWHAZ +0.64, USP22 +0.59, FZD7 +0.58, TPD52 +0.45 (Up); KLF10 −1.44, DUSP10 −1.14, HES1 −0.98, MAP3K2 −0.51, **CPT1A −0.48**, **NEK6 −0.40**, LRP5 −0.29 (Down); CES1, INSIG1, HNF4G, MTSS1, TP53INP1, GPAM (ns).

## Outstanding Issues / Next Steps

- **Confirm ASoC donor ancestry.** LD was tested only in EUR and EAS. If donors are AFR/AMR/SAS, re-run the 6p21 pairwise LD in the matched 1000 Genomes population (MHC LD is population-specific).
- **Optionally test additional HCC lead SNPs** for proxy relationships (e.g., HLA-DQ/DP leads `rs9272105`, `rs3077`) against the 6p21 ASoC cluster.
- **Investigate the CPT1A and NEK6 discordances** (literature ↑ vs measured ↓ in bulk TCGA-LIHC): check GTEx-normal baseline (GEPIA2), protein-level data (HPA), and whether the effect is HCC-subtype-specific or masked by non-parenchymal cells.
- **Integrate allelic direction.** Combine ASoC allelic imbalance (`percRef`, ref/alt counts) and motifbreakR allele effect (`mb_max_abs_alleleDiff`) with the measured tumor-vs-normal direction to predict whether the accessibility-increasing allele is pro- or anti-tumor.
- **Cross-reference the gene list against curated HCC gene sets** (COSMIC Cancer Gene Census, DisGeNET "Carcinoma, Hepatocellular").
- **Motif-enrichment test:** assess whether ASoC-disrupted TF motifs (HNF4A/HNF4G, FOXA, RXR) are enriched among the significant variants.

## Context for LLM Handoff

This session analyzed an allele-specific open chromatin (ASoC) table from human hepatocytes (`Steffi_works/sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv`, GRCh38) to find variant- and gene-level associations with hepatocellular carcinoma (HCC), broader liver disease/traits, and gastric cancer, and then to quantify tumor-vs-normal expression of the implicated genes. The deliverable is `Steffi_works/analysis/ASoC_hepatocyte_liver_cancer_association_summary.md`. Variant level: ClinVar/dbSNP, Europe PMC (by exact rsID), and GWAS Catalog found **no direct variant→HCC association**; only two variant-level signals exist — `rs2060982` (CPT1A → hepatic lipid/fatty-acid GWAS traits) and `rs9898809` (TOB1/TOB1-AS1 → gastric cancer). The 6p21 MHC ASoC cluster (`rs16899941` etc.) has a strong **autoimmune** GWAS signal (hypothyroidism/RA), and a formal Ensembl LD-proxy test (public REST API, EUR and EAS) **rejected** it being an LD proxy of the HCC MICA locus `rs2596542`. Gene level: the analysis narrowed to genes whose ASoC SNPs show strong TF-motif disruption (`mb_n_strong` > 0; removed TRIB1, TOB1, ING1, TLR5), then attached **measured** TCGA-LIHC tumor-vs-normal fold-changes from a scripted edgeR workflow (UCSC Xena GDC STAR counts; 374 tumor vs 50 normal; |log2FC|≥0.2, FDR<0.05; scripts in `Steffi_works/analysis/tcga_lihc_de/`). Measured Up in HCC: CHEK1, HDAC11, LAMC1, NDRG1, CAPN2, GSDME, YWHAZ, USP22, FZD7, TPD52; measured Down: KLF10, DUSP10, HES1, MAP3K2, **CPT1A**, **NEK6**, LRP5; not significant: CES1, INSIG1, HNF4G, MTSS1, TP53INP1, GPAM. Notable discordances vs literature are **CPT1A** and **NEK6** (literature ↑ but measured ↓ in bulk tumor). Genome-wide, 111/34 of 230 ASoC genes are up/down. Open items: confirm donor ancestry for LD, probe the CPT1A/NEK6 discordances (GTEx/HPA/subtype), and integrate ASoC allelic direction with tumor-vs-normal direction. All external queries used reproducible Europe PMC, GWAS Catalog, Ensembl LD REST, and UCSC Xena GDC endpoints documented above; the R workflow requires only base Bioconductor packages (edgeR/limma/data.table), no `TCGAbiolinks`.
