# ASoC Variants in Hepatocytes — Association with Hepatocellular Carcinoma / Liver Neoplasms

**Input file:** `Steffi_works/sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated.tsv`
**Analysis type:** Allele-Specific Open Chromatin (ASoC) analysis in hepatocytes. Each row is a significant heterozygous SNP showing allelic imbalance in chromatin accessibility, annotated with its overlapping/nearest gene (ENSEMBL, SYMBOL) and predicted motif disruptions.
**Question addressed:** Which variants (`variantID`) and which genes (`ENSEMBL` / `SYMBOL`) have a known connection to **hepatocellular carcinoma (HCC)** or, more broadly, **liver-related neoplasms**?

> **Important interpretive caveat.** These ASoC SNPs are common regulatory variants (promoter / intronic / intergenic), not coding pathogenic mutations. Direct checks against dbSNP/ClinVar (e.g. `rs2980223`, `rs2060982`), an **expanded Europe PMC literature search by exact rsID**, and **GWAS Catalog** lookups all show **no direct variant→hepatocellular-carcinoma association** for any variant in this table — as expected for regulatory eQTL-type SNPs. Therefore the disease relevance below is established **at the gene level** (the gene each variant regulates), not by a curated variant-disease claim. Any HCC link should be read as "the regulated gene is implicated in liver cancer / hepatocyte biology," which is the biologically meaningful readout for an ASoC study. See **Section 1b** for the expanded online search (literature, GWAS Catalog, LD-proxy reasoning).

---

## 1) Variant-level (`variantID`) findings

No variant in this table carries a **direct, curated variant→hepatocellular-carcinoma association** in dbSNP/ClinVar. The variants are regulatory (mostly `Promoter (<=1kb)`), and their cancer relevance is inherited from the target gene (Section 2). Representative variants checked:

| variantID                 | Gene                 | Annotation        | Variant-level cancer annotation               |
| ------------------------- | -------------------- | ----------------- | --------------------------------------------- |
| rs2980223                 | HNF4G (upstream)     | Distal Intergenic | Not reported in ClinVar                       |
| rs2060982                 | CPT1A (2kb upstream) | Promoter (1–2kb)  | Not reported in ClinVar                       |
| rs1017105790 / rs58065091 | NDRG1                | Promoter / Distal | Not reported in ClinVar (gene-level HCC link) |
| rs10097617                | TP53INP1             | Promoter          | Not reported in ClinVar (gene-level p53 link) |

**Conclusion for (1):** the value of this dataset is regulatory — the SNPs modulate accessibility at cancer-relevant hepatocyte genes rather than being classified disease alleles themselves.

---

## 1b) Expanded online search — literature, GWAS Catalog, and LD-proxy reasoning

Because very few individual SNPs are formally documented for HCC (as the ASoC variants are common regulatory SNPs), the search was expanded beyond ClinVar to **(i)** full-text literature (Europe PMC), **(ii)** the **NHGRI-EBI GWAS Catalog** (queried by rsID), and **(iii)** proximity/LD reasoning against known HCC GWAS loci (coordinates are GRCh38).

### (i) Literature search by exact rsID (Europe PMC)

Direct text-mining of the actual `variantID`s returned **essentially no hepatocellular-carcinoma hits** — confirming the premise that these specific variants are undocumented for HCC:

| rsID(s) queried (representative)                                   | Literature hits for liver/HCC/cancer                                                        | Note                                                     |
| ------------------------------------------------------------------ | ------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| rs2471847, rs9935289, rs78255464, rs2270200, rs2385160, rs11020802 | 0 relevant (2 stray matches were unrelated lung-cancer/ATM papers)                          | No direct HCC association                                |
| rs9898809, rs2231516, rs10283129, rs62527000, rs7775082, rs2517765 | 1 relevant → **TOB1 gene** (near rs9898809 / `TOB1-AS1`) SNPs studied in **gastric cancer** | Gene-level cancer link only, not HCC, not this exact SNP |

→ **No ASoC variant in this table has a published, direct variant→HCC association.** The one signal that surfaced is at the **TOB1 / TOB1-AS1** locus (rs9898809), where *TOB1* is a documented tumor-suppressor whose polymorphisms have been studied in gastric cancer.

### (ii) GWAS Catalog lookups (by rsID)

Most variants have **no GWAS Catalog entry**. The notable exception is the cluster of **MHC / chromosome-6p21 ASoC SNPs**:

| rsID                   | Locus                                             | GWAS Catalog result                                                                           | HCC relevance                                         |
| ---------------------- | ------------------------------------------------- | --------------------------------------------------------------------------------------------- | ----------------------------------------------------- |
| **rs16899941**         | 6p21.33 (`GPANK1`/`BAG6`/`CSNK2B`, MHC class III) | Genome-wide significant (**p ≈ 3×10⁻⁵²**) for **hypothyroidism** and **rheumatoid arthritis** | Strong signal, but **autoimmune**, *not* liver cancer |
| rs2517765, rs1632906   | 6p22.1 (`HCG4B`, near HLA-A/HLA-F)                | No direct HCC association returned                                                            | MHC region (see proxy note)                           |
| rs78905038, rs30556713 | 6p21.33 (`GNL1`, MHC)                             | No direct HCC association                                                                     | MHC region                                            |

→ The 6p21 MHC ASoC hits are **real, functionally active immune-locus variants**, but their catalogued genome-wide associations are with **autoimmune/thyroid traits, not HCC**. This is an important *negative* clarification.

### (iii) Proxy-SNP / known-HCC-locus proximity (GRCh38)

Established HCC GWAS risk loci (esp. in HBV-related HCC) and their GRCh38 positions were compared against the ASoC coordinates:

| Known HCC GWAS locus                | Lead SNP                       | GRCh38 position     | Nearest ASoC SNP(s) in this table | Approx. distance                                                 |
| ----------------------------------- | ------------------------------ | ------------------- | --------------------------------- | ---------------------------------------------------------------- |
| MHC class I / MICA region (6p21.33) | rs2596542 (MICA)               | chr6 ≈ 31.4 Mb      | **rs16899941** (chr6:31,665,912)  | ~0.13 Mb (proximity only — LD-proxy tested & rejected, see (iv)) |
| HLA-DQ/DP (6p21.32)                 | rs9272105 / rs9275319 / rs3077 | chr6 ≈ 32.5–33.0 Mb | rs16899941 (31.67 Mb)             | ~0.8–1.3 Mb (likely outside a single LD block)                   |
| STAT4 (2q32.2)                      | rs7574865                      | chr2 ≈ 191.1 Mb     | none nearby                       | —                                                                |
| KIF1B (1p36.22)                     | rs17401966                     | chr1 ≈ 10.2 Mb      | none nearby                       | —                                                                |
| GRIK1 (21q21)                       | —                              | chr21               | none nearby                       | —                                                                |

→ Only the **6p21 MHC cluster** falls near a documented HCC GWAS locus (the MICA region). This proxy hypothesis was then **tested formally** (below).

### (iv) Formal LD-proxy test — Ensembl LD REST API (EUR + EAS), **result: rejected**

The Ensembl LD REST API is **public and requires no token** — only a 1000 Genomes population. Using **`1000GENOMES:phase_3:EUR`** and **`:EAS`**, pairwise LD was computed between each 6p21 MHC ASoC SNP and the HCC MICA lead SNP **rs2596542**:

| ASoC SNP   | vs HCC lead      | Pairwise r² (EUR) | Pairwise r² (EAS) | Proxy (r²≥0.8)? |
| ---------- | ---------------- | ----------------- | ----------------- | --------------- |
| rs16899941 | rs2596542 (MICA) | none returned     | none returned     | **No**          |
| rs2517765  | rs2596542 (MICA) | none returned     | none returned     | **No**          |
| rs1632906  | rs2596542 (MICA) | none returned     | none returned     | **No**          |
| rs78905038 | rs2596542 (MICA) | none returned     | none returned     | **No**          |

`GET https://rest.ensembl.org/ld/human/pairwise/{rsID}/rs2596542?population_name=1000GENOMES:phase_3:{EUR|EAS}` returned an **empty array `[]`** for all four SNPs in **both** populations, i.e. they are **not in appreciable LD** with the MICA HCC signal. Consistent with this, the EUR proxy sets of `rs16899941` (e.g. rs2077102, rs2242657/8, rs2295665 — all r²≈1) and of `rs2596542` (e.g. rs2428478, rs2596562, rs2596483 — all r²≈0.93–1) **do not overlap**.

→ **Conclusion: the 6p21 MHC ASoC SNPs are NOT LD proxies of the known HCC MICA locus in either EUR or EAS.** Their genome-wide catalogued signal remains autoimmune (hypothyroidism / rheumatoid arthritis). The ~0.13 Mb genomic proximity reflects the gene-dense extended MHC, not shared LD. This removes the last candidate variant-level HCC link and reinforces that the dataset's HCC relevance is **gene-level**, not variant-level.

> Note: LD was tested in **EUR and EAS** (the two most relevant reference panels). Both reject the proxy; if the ASoC donors are of another ancestry (e.g. AFR/AMR/SAS), re-run with that population, since MHC LD structure is population-specific.


---

## 1c) Broadened search — liver-related (any) and gastric-cancer associations

The search was widened beyond HCC to **any liver-related disease/trait** (NAFLD/NASH, steatosis, cirrhosis, hepatitis, bilirubin, liver enzymes/aminotransferases, cholestasis, lipid/fatty-acid metabolism) and to **gastric cancer**, using Europe PMC (text-mined by exact rsID) and the NHGRI-EBI GWAS Catalog (structured EFO traits by rsID).

### Positive hits (variant-level)

| variantID      | Gene              | Source           | Trait / disease                                 | Significance                     | Category                                            |
| -------------- | ----------------- | ---------------- | ----------------------------------------------- | -------------------------------- | --------------------------------------------------- |
| **rs2060982**  | CPT1A (upstream)  | **GWAS Catalog** | **saturated fatty acids / total fatty acids %** | p = 3×10⁻³⁰ (β decrease)         | Liver-related (hepatic lipid/fatty-acid metabolism) |
| **rs2060982**  | CPT1A (upstream)  | **GWAS Catalog** | **triglycerides in IDL**                        | p = 2×10⁻¹³ (β decrease)         | Liver-related (lipid metabolism)                    |
| **rs9898809**  | TOB1 / `TOB1-AS1` | Europe PMC       | **gastric cancer risk** (Chinese Han)           | case–control assoc. (PMC5929081) | Gastric cancer                                      |
| **rs16899941** | GPANK1/BAG6 (MHC) | GWAS Catalog     | hypothyroidism; rheumatoid arthritis            | p ≈ 3×10⁻⁵²                      | Autoimmune (not liver, not gastric) — see 1b        |

→ **CPT1A / rs2060982** is the clearest new signal: it is a bona-fide GWAS variant for **hepatic lipid / fatty-acid metabolism** quantitative traits. CPT1A is the rate-limiting enzyme of mitochondrial fatty-acid β-oxidation in the liver, and this locus is repeatedly linked to circulating fatty-acid/lipid phenotypes (relevant to NAFLD/steatosis physiology, though the catalogued traits are quantitative metabolite measures, not a liver-disease diagnosis).

### Negative results (important for scoping)

- **Liver disease (non-cancer):** Text-mining **all ~140 rsIDs** against liver terms returned only **2–3 incidental hits** (a hepatic glucocorticoid-receptor review, an X-chromosome eQTL paper, and the Xue et al. 2018 type-2-diabetes GWAS) — none is a specific liver-disease association for these variants.
- **GWAS Catalog:** Besides rs2060982 (lipids) and rs16899941 (autoimmune), the tested SNPs (incl. HNF4G rs2980223, KLF10 rs2471847) have **no catalogued liver-disease or gastric-cancer associations**.
- **Gastric cancer:** Across all rsIDs, the **only** hit is **rs9898809 / TOB1** (as in the HCC search).

**Net:** broadening to liver-related and gastric-cancer traits adds exactly **two** variant-level signals — **rs2060982 (CPT1A → lipid/fatty-acid metabolism, liver-relevant)** and **rs9898809 (TOB1 → gastric cancer)** — while confirming that the remaining variants have no documented liver or gastric disease associations. This is consistent with the earlier conclusion: the dataset's disease relevance is predominantly **gene-level**, and the few variant-level GWAS signals point to **metabolic/immune quantitative traits rather than neoplasia**.

---

## 2) Gene-level (`ENSEMBL` / `SYMBOL`) findings

The following genes overlapping/regulated by significant ASoC SNPs have documented roles in **hepatocellular carcinoma or liver / hepatocyte biology**. Ranked roughly by strength and directness of the liver-cancer link. The **Lit. (HCC)** column gives the approximate number of Europe PMC records co-mentioning the gene with hepatocellular carcinoma (queried 2026-07-26) — a proxy for depth of published evidence.

### Strong / direct liver-cancer relevance

| SYMBOL       | ENSEMBL         | Variant(s)               | Lit. (HCC)  | Relevance to HCC / liver neoplasm                                                                                                                                                                                        |
| ------------ | --------------- | ------------------------ | ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **HNF4G**    | ENSG00000164749 | rs2980223                | ~187        | Hepatocyte Nuclear Factor 4 gamma — master hepatocyte transcription factor of the HNF4 family; HNF4 signaling is central to hepatocyte identity and is dysregulated in HCC (loss of differentiation, tumor progression). |
| **NDRG1**    | ENSG00000104419 | rs1017105790, rs58065091 | ~1123       | N-myc Downstream Regulated 1 — well-studied in HCC; associated with tumor progression, metastasis, hypoxia response, and patient prognosis in hepatocellular carcinoma.                                                  |
| **TP53INP1** | ENSG00000164938 | rs10097617               | ~653        | Tumor Protein p53 Inducible Nuclear Protein 1 — p53 target/stress protein; tumor-suppressor context; downregulation reported in HCC and other GI cancers.                                                                |
| **CPT1A**    | ENSG00000110090 | rs2060982                | (metabolic) | Carnitine Palmitoyltransferase 1A — rate-limiting enzyme of fatty-acid β-oxidation; lipid-metabolism reprogramming via CPT1A is implicated in HCC growth and progression (metabolic rewiring in liver cancer).           |
| **INSIG1**   | ENSG00000186480 | rs576875943              | (metabolic) | Insulin-Induced Gene 1 — regulator of SREBP-mediated lipogenesis; lipid/cholesterol metabolism strongly tied to hepatocarcinogenesis and NAFLD→HCC.                                                                      |
| **KLF10**    | ENSG00000155090 | rs2471847                | moderate    | Krüppel-Like Factor 10 (TIEG1) — TGF-β–inducible tumor-suppressor-type TF; implicated in hepatocyte proliferation control and HCC.                                                                                       |

### Moderate / plausible liver-cancer or hepatocyte-cancer relevance

| SYMBOL      | ENSEMBL         | Variant(s)             | Lit. (HCC) | Relevance                                                                                                                                                                                                                         |
| ----------- | --------------- | ---------------------- | ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **FZD7**    | ENSG00000155760 | rs144271591            | ~755       | Frizzled-7 — Wnt/β-catenin receptor; FZD7 is one of the most consistently overexpressed Wnt receptors in HCC and drives β-catenin activation.                                                                                     |
| **MTSS1**   | ENSG00000170873 | rs934323               | ~311       | Metastasis Suppressor 1 — actin-regulatory tumor/metastasis suppressor; loss associated with HCC invasion/metastasis.                                                                                                             |
| **TOB1**    | ENSG00000141232 | rs9898809 (`TOB1-AS1`) | gastric    | Transducer of ERBB2 — antiproliferative tumor suppressor; the **one locus in this table with published cancer-SNP literature** (TOB1 polymorphisms studied in gastric cancer). Antisense RNA `TOB1-AS1` is the annotated feature. |
| **CES1**    | ENSG00000198848 | rs878890956            | liver-spec | Carboxylesterase 1 — highly liver-enriched drug/lipid-metabolizing enzyme; altered expression in HCC and hepatic metabolic disease.                                                                                               |
| **DUSP10**  | ENSG00000143507 | rs4317817              | context    | Dual-Specificity Phosphatase 10 (MKP5) — MAPK/JNK pathway regulator; context-dependent roles in liver and GI cancer.                                                                                                              |
| **CAPN2**   | ENSG00000162909 | rs2154114              | context    | Calpain-2 — protease implicated in tumor invasion/migration including HCC.                                                                                                                                                        |
| **GSDME**   | ENSG00000105928 | rs796754117            | context    | Gasdermin E (DFNA5) — pyroptosis executioner; tumor-suppressor/immune context studied across cancers including liver.                                                                                                             |
| **LAMC1**   | ENSG00000135862 | rs112409204            | context    | Laminin subunit gamma-1 — ECM component; upregulation linked to HCC invasion and poor prognosis.                                                                                                                                  |
| **TLR5**    | ENSG00000187554 | rs2302597              | context    | Toll-Like Receptor 5 — innate-immune/inflammation signaling relevant to inflammation-driven hepatocarcinogenesis.                                                                                                                 |
| **CYP20A1** | ENSG00000119004 | rs1037204785           | context    | Cytochrome P450 20A1 — CYP-family; hepatic xenobiotic-metabolism context.                                                                                                                                                         |
| **TPD52**   | ENSG00000076554 | rs10111451             | context    | Tumor Protein D52 — proliferation-associated oncogene amplified in several carcinomas.                                                                                                                                            |

### Notable liver-metabolism / hepatocyte-identity genes (indirect)

- **MAP3K2** (ENSG00000169967; rs889055704, rs1050847260, rs13034125, rs13034122) — MAPK cascade kinase; recurrent ASoC hits, general oncogenic signaling.
- **NDUFA4 / MTCH1 / MRPL36 / PITRM1 / PUS1 / TARS3** — mitochondrial/metabolic genes; metabolic-reprogramming context in HCC but no direct disease claim.
- **SLC22A23**, **UBE2V2**, **USP22** (ENSG00000124422; rs140762989) — USP22 is a well-known cancer stemness/oncogenic deubiquitinase reported in HCC (~699 Europe PMC HCC co-mentions).

---

## 3) Gene-level differential expression in HCC vs normal liver

This section reports **measured tumor-vs-normal differential expression** for the ASoC-linked genes, computed from **TCGA-LIHC** RNA-seq (not literature inference). **Up** = higher in tumor; **Down** = lower in tumor.

**Data & method (scripted, reproducible).** STAR gene-level counts for TCGA-LIHC were downloaded from the UCSC Xena GDC hub (`TCGA-LIHC.star_counts.tsv.gz`, log2(count+1); reversed to raw counts), samples classified by TCGA barcode into **tumor (n = 374, sample code 01)** vs **solid-tissue normal (n = 50, code 11)**, and differential expression computed with **edgeR** (`filterByExpr` → TMM normalization → `estimateDisp` → quasi-likelihood F-test, tumor vs normal). Thresholds: **|log2FC| ≥ 0.2 and FDR (BH) < 0.05**. Scripts and outputs: [analysis/tcga_lihc_de/run_tcga_lihc_de.R](tcga_lihc_de/run_tcga_lihc_de.R), full table [analysis/tcga_lihc_de/TCGA-LIHC_tumor_vs_normal_DE_full.tsv](tcga_lihc_de/TCGA-LIHC_tumor_vs_normal_DE_full.tsv), ASoC intersection [analysis/tcga_lihc_de/ASoC_genes_TCGA-LIHC_DE.tsv](tcga_lihc_de/ASoC_genes_TCGA-LIHC_DE.tsv).

**Genome-wide ASoC-gene summary.** Of **230** ASoC-annotated genes with an Ensembl ID, TCGA-LIHC DE classifies **111 Up, 34 Down, 57 not-significant, 28 not-tested** (filtered out by `filterByExpr`).

> **Filter applied (motifbreakR strong TF disruption).** The panel below is **narrowed to genes whose ASoC SNP(s) show confirmed strong TF-motif disruption** (`mb_n_strong` > 0). Genes with `mb_n_strong = 0` for every ASoC SNP were removed (**TRIB1, TOB1/TOB1-AS1, ING1, TLR5**). `mb_n_strong` = maximum across the gene's ASoC SNPs (e.g. NDRG1 via rs58065091 = 15; MAP3K2 via rs13034125 = 3).

### Measured UP in TCGA-LIHC (log2FC ≥ +0.2, FDR < 0.05)

| SYMBOL     | ENSEMBL         | ASoC SNP (strong) | mb_n_strong | **log2FC** | **FDR**   | Literature (Section 2) concordance |
| ---------- | --------------- | ----------------- | ----------- | ---------- | --------- | ---------------------------------- |
| **CHEK1**  | ENSG00000149554 | rs555752          | 1           | **+2.22**  | 2.3×10⁻²⁴ | concordant (oncogenic, ↑)          |
| **HDAC11** | ENSG00000163517 | rs2290193         | 3           | **+1.76**  | 9.3×10⁻²⁷ | resolves literature "context" → ↑  |
| **LAMC1**  | ENSG00000135862 | rs112409204       | 2           | **+1.39**  | 3.3×10⁻¹⁸ | concordant (↑ invasion)            |
| **NDRG1**  | ENSG00000104419 | rs58065091        | 15          | **+1.19**  | 2.3×10⁻⁷  | concordant (↑ progression)         |
| **CAPN2**  | ENSG00000162909 | rs2154114         | 2           | **+0.76**  | 1.4×10⁻⁸  | concordant (↑)                     |
| **GSDME**  | ENSG00000105928 | rs796754117       | 3           | **+0.70**  | 2.7×10⁻³  | resolves "context" → ↑             |
| **YWHAZ**  | ENSG00000164924 | rs56165483        | 5           | **+0.64**  | 1.0×10⁻⁷  | concordant (↑)                     |
| **USP22**  | ENSG00000124422 | rs140762989       | 6           | **+0.59**  | 1.4×10⁻⁵  | concordant (oncogenic, ↑)          |
| **FZD7**   | ENSG00000155760 | rs144271591       | 2           | **+0.58**  | 2.7×10⁻²  | concordant (Wnt, ↑)                |
| **TPD52**  | ENSG00000076554 | rs10111451        | 2           | **+0.45**  | 8.7×10⁻⁴  | concordant (↑)                     |

### Measured DOWN in TCGA-LIHC (log2FC ≤ −0.2, FDR < 0.05)

| SYMBOL     | ENSEMBL         | ASoC SNP (strong)   | mb_n_strong | **log2FC** | **FDR**   | Literature (Section 2) concordance      |
| ---------- | --------------- | ------------------- | ----------- | ---------- | --------- | --------------------------------------- |
| **KLF10**  | ENSG00000155090 | rs2471847           | 1           | **−1.44**  | 2.1×10⁻²⁴ | concordant (tumor suppressor, ↓)        |
| **DUSP10** | ENSG00000143507 | rs4317817           | 1           | **−1.14**  | 3.9×10⁻¹⁶ | resolves "context" → ↓                  |
| **HES1**   | ENSG00000114315 | rs6444772 (distal)  | 16          | **−0.98**  | 8.0×10⁻¹⁵ | resolves "context" → ↓                  |
| **MAP3K2** | ENSG00000169967 | rs13034125          | 3           | **−0.51**  | 3.3×10⁻¹² | resolves "context" → ↓                  |
| **CPT1A**  | ENSG00000110090 | rs2060982           | 1           | **−0.48**  | 1.5×10⁻³  | **discordant** (lit ↑/FAO → measured ↓) |
| **NEK6**   | ENSG00000119408 | rs72759286 (distal) | 2           | **−0.40**  | 1.2×10⁻⁵  | **discordant** (lit ↑ → measured ↓)     |
| **LRP5**   | ENSG00000162337 | rs312024 (intron)   | 2           | **−0.29**  | 1.5×10⁻²  | resolves "context" → ↓                  |

### Not significant in TCGA-LIHC (|log2FC| < 0.2 or FDR ≥ 0.05)

| SYMBOL       | ENSEMBL         | ASoC SNP (strong)   | mb_n_strong | log2FC | FDR   | Note                                                |
| ------------ | --------------- | ------------------- | ----------- | ------ | ----- | --------------------------------------------------- |
| **CES1**     | ENSG00000198848 | rs878890956         | 6           | −0.32  | 0.19  | trend ↓ (lit ↓) but not significant                 |
| **INSIG1**   | ENSG00000186480 | rs576875943         | 4           | −0.40  | 0.090 | trend ↓, not significant                            |
| **HNF4G**    | ENSG00000164749 | rs2980223           | 10          | −0.20  | 0.20  | lit ↑ but measured ns (slight ↓) — discordant trend |
| **MTSS1**    | ENSG00000170873 | rs934323            | 2           | +0.07  | 0.58  | lit ↓ but measured ns                               |
| **TP53INP1** | ENSG00000164938 | rs10097617          | 21          | +0.02  | 0.90  | lit ↓ but measured ns                               |
| **GPAM**     | ENSG00000119927 | rs75027818 (distal) | 5           | −0.01  | 0.98  | ns                                                  |


### Data sources & verification (recommended for quantitative confirmation)

The measured directions above come from a **scripted TCGA-LIHC edgeR analysis** (this repository; see script/output links at the top of Section 3). For independent visual confirmation of any gene (boxplots, GTEx normal comparison), the same TCGA-LIHC data can be queried by gene symbol at:

- **GEPIA2** (TCGA-LIHC + GTEx normal): http://gepia2.cancer-pku.cn — "Expression DIY → Boxplot", dataset **LIHC**.
- **UALCAN** (TCGA-LIHC): https://ualcan.path.uab.edu/analysis.html — TCGA → Liver hepatocellular carcinoma.
- **TNMplot** (tumor vs normal, RNA-seq): https://tnmplot.com/analysis/ — Liver.
- **Human Protein Atlas Pathology** (per gene, protein/RNA): `https://www.proteinatlas.org/{ENSEMBL}/pathology` (e.g. https://www.proteinatlas.org/ENSG00000155760/pathology for FZD7).

**Reproduce:** `Rscript analysis/tcga_lihc_de/run_tcga_lihc_de.R` (downloads the Xena GDC matrix, runs edgeR, writes the full and ASoC-intersected DE tables). Only base Bioconductor packages (edgeR, limma, data.table, SummarizedExperiment) are required — no `TCGAbiolinks` install needed.

---

## 4) Cross-reference against curated HCC gene sets (COSMIC CGC + DisGeNET)

The **230** ASoC-annotated genes (Ensembl IDs) were cross-referenced against two curated hepatocellular-carcinoma gene sets. Both primary resources are now access-gated (the COSMIC Cancer Gene Census download requires a Sanger login; the DisGeNET API requires an API key), so **tokenless, reproducible mirrors** were used:

- **COSMIC Cancer Gene Census (HCC-mapped)** — taken from the **Open Targets Platform** GraphQL API as the `cancer_gene_census` evidence datasource for **hepatocellular carcinoma (`MONDO_0007256`)**. Open Targets ingests the census and maps each gene to its census tumour types, so a nonzero `cancer_gene_census` score for HCC ≈ "CGC gene whose census annotation includes liver/HCC."
- **DisGeNET "Carcinoma, Hepatocellular" (`C2239176`)** — retrieved via the **Enrichr** DisGeNET library (term **"Liver carcinoma"** = `C2239176`; **3,593** genes, text-mining + curated evidence combined), which is tokenless and not rate-limited. Because this set is text-mining-inclusive it is broad/low-specificity: membership means the gene is co-reported with HCC in DisGeNET, not necessarily a curated causal driver. A second script (below) instead queries the **native DisGeNET API** for the per-gene GDA `score` + Evidence Index, but that endpoint needs an API key and rate-limits hard (TRIAL `retry-after` can be ~22 h), so the tokenless Enrichr numbers are reported here.

Two interchangeable scripts produce the DisGeNET half (COSMIC-CGC half is identical in both):
- Enrichr (tokenless, used for the numbers below): [analysis/hcc_geneset_crossref/crossref_hcc_genesets_enrichr.R](hcc_geneset_crossref/crossref_hcc_genesets_enrichr.R) → [ASoC_genes_HCC_curated_crossref_enrichr.tsv](hcc_geneset_crossref/ASoC_genes_HCC_curated_crossref_enrichr.tsv)
- Native DisGeNET GDA score/EI (API key, rate-limited; LSF wrapper `bsub_crossref_hcc_genesets.sh`): [analysis/hcc_geneset_crossref/crossref_hcc_genesets_disgenet.R](hcc_geneset_crossref/crossref_hcc_genesets_disgenet.R) → `ASoC_genes_HCC_curated_crossref_disgenet.tsv`

**Result.** Of 230 ASoC genes: **5** are in the (HCC-mapped) **COSMIC CGC**, **33** are in **DisGeNET C2239176**, **1** is in **both** (NDRG1), and **37** are in **at least one** curated HCC set. The remaining 193 ASoC genes are not in either curated HCC set — consistent with Sections 1–1c: these are regulatory variants near hepatocyte genes, most of which are not catalogued cancer drivers.

### 4a) COSMIC Cancer Gene Census — HCC-mapped hits (5)

| SYMBOL    | ENSEMBL         | CGC(HCC) score | DisGeNET C2239176 | TCGA-LIHC dir (log2FC) | OT HCC assoc | Strong-TF panel |
| --------- | --------------- | -------------- | ----------------- | ---------------------- | ------------ | --------------- |
| **NDRG1** | ENSG00000104419 | 0.53           | ✔ (both)          | **Up** (+1.19)         | 0.35         | ✔               |
| **PBX1**  | ENSG00000185630 | 0.53           | —                 | ns (−0.32)             | 0.33         | —               |
| **ELK4**  | ENSG00000158711 | 0.46           | —                 | Up (+0.34)             | 0.29         | —               |
| **PLEC**  | ENSG00000178209 | 0.30           | —                 | Up (+0.71)             | 0.21         | —               |
| **N4BP2** | ENSG00000078177 | 0.30           | —                 | Down (−0.68)           | 0.21         | —               |

→ **NDRG1** is the single strongest result: it is in **both** curated sets, is measured **Up** in TCGA-LIHC, and carries a strong TF-motif disruption (`mb_n_strong` = 15 via rs58065091) — the most credible variant→gene→HCC candidate in the dataset. **PBX1** and **ELK4** are census genes by virtue of documented oncogenic **fusions/translocations** (their census tumour-type mapping reaches HCC in Open Targets); **PLEC** and **N4BP2** round out the CGC-HCC set.

### 4b) DisGeNET "Carcinoma, Hepatocellular" (C2239176) hits (33)

Sorted by Open Targets HCC association score (`✔` in the *Strong-TF panel* column = ASoC SNP with `mb_n_strong` > 0, i.e. the Section-3 panel):

| SYMBOL       | TCGA-LIHC dir (log2FC) | OT HCC assoc | Strong-TF panel | Also COSMIC CGC |
| ------------ | ---------------------- | ------------ | --------------- | --------------- |
| **NDRG1**    | Up (+1.19)             | 0.35         | ✔               | ✔               |
| **CHEK1**    | Up (+2.22)             | 0.12         | ✔               | —               |
| **USP22**    | Up (+0.59)             | 0.12         | ✔               | —               |
| **HDAC11**   | Up (+1.76)             | 0.11         | ✔               | —               |
| **YWHAZ**    | Up (+0.64)             | 0.11         | ✔               | —               |
| **SDC1**     | Down (−0.29)           | 0.11         | —               | —               |
| **MTSS1**    | ns (+0.07)             | 0.11         | ✔               | —               |
| **LAMC1**    | Up (+1.39)             | 0.10         | ✔               | —               |
| **TP53INP1** | ns (+0.02)             | 0.10         | ✔               | —               |
| **NEK6**     | Down (−0.40)           | 0.10         | ✔               | —               |
| **FZD7**     | Up (+0.58)             | 0.09         | ✔               | —               |
| **TFAM**     | Down (−0.38)           | 0.09         | —               | —               |
| **DCAF13**   | Up (+1.05)             | 0.09         | —               | —               |
| **CES1**     | ns (−0.32)             | 0.09         | ✔               | —               |
| **GALNT2**   | Down (−0.34)           | 0.08         | —               | —               |
| **MAP3K2**   | Down (−0.51)           | 0.08         | ✔               | —               |
| **CD1D**     | Down (−1.59)           | 0.08         | —               | —               |
| **TPD52**    | Up (+0.45)             | 0.08         | ✔               | —               |
| **PEX5**     | ns (−0.04)             | 0.08         | —               | —               |
| **ING1**     | ns (−0.04)             | 0.07         | —               | —               |
| **KLF10**    | Down (−1.44)           | 0.06         | ✔               | —               |
| **INSIG1**   | ns (−0.40)             | 0.06         | ✔               | —               |
| **HES1**     | Down (−0.98)           | 0.06         | ✔               | —               |
| **MOK**      | Up (+0.41)             | 0.04         | —               | —               |
| **SPG7**     | Down (−0.28)           | 0.04         | —               | —               |
| **XPR1**     | Up (+0.64)             | 0.03         | —               | —               |
| **FSD1L**    | Up (+1.17)             | 0.03         | —               | —               |
| **HNF4G**    | ns (−0.20)             | 0.03         | ✔               | —               |
| **MTR**      | Up (+0.76)             | 0.03         | —               | —               |
| **SP3**      | ns (−0.04)             | 0.02         | —               | —               |
| **PADI1**    | Up (+1.95)             | 0.02         | —               | —               |
| **LAMTOR2**  | Up (+0.63)             | 0.00         | —               | —               |
| **SH3BP4**   | ns (−0.07)             | 0.00         | —               | —               |

→ **17** of the 33 DisGeNET hits are also in the **strong-TF-disruption panel** (`mb_n_strong` > 0), sharpening the priority list: **NDRG1, CHEK1, USP22, HDAC11, YWHAZ, MTSS1, LAMC1, TP53INP1, NEK6, FZD7, CES1, MAP3K2, TPD52, KLF10, INSIG1, HES1, HNF4G**. These are ASoC genes that are simultaneously (i) HCC-associated in a curated set, (ii) measured/known differentially expressed, and (iii) regulated by a variant predicted to strongly alter TF binding.

### 4c) Per-variant join and manually-curated priority shortlist (22 variants)

The per-gene curated-HCC cross-reference (4a/4b) was joined back onto the **per-variant** ASoC table (matched on base Ensembl ID) so that every SNP row carries the curated-HCC columns (`COSMIC_CGC_HCC`, `DisGeNET_HCC`, `curated_HCC_hit`, `OT_*`, alongside the existing `TCGA_LIHC_dir/log2FC` and `mb_n_strong`). Script: [analysis/hcc_geneset_crossref/join_crossref_to_snps.R](hcc_geneset_crossref/join_crossref_to_snps.R) → [sig_ASoC_in_Hepatocyte_annotated_HCC_crossref.tsv](../sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated_HCC_crossref.tsv). Of the **315** hepatocyte ASoC variant rows, **47** fall in ≥1 curated HCC set (43 DisGeNET C2239176, 8 COSMIC CGC-HCC) — more rows than genes because several genes carry multiple ASoC SNPs.

From that joined table, a **manually-curated priority shortlist of 22 variants** was selected (after all automated analysis finished) — variants that are both **in the DisGeNET C2239176 HCC set** and carry **strong TF-motif disruption** (`mb_n_strong` > 0), with lower-confidence/redundant rows removed by manual review. Output: [sig_ASoC_in_Hepatocyte_annotated_HCC_crossref_DisGeNET_yes.tsv](../sig_ASoC_by_celltype/sig_ASoC_in_Hepatocyte_annotated_HCC_crossref_DisGeNET_yes.tsv) — **22 variant rows across 17 genes**, sorted by ASoC significance (`adjPVal`).

| variantID    | SYMBOL | Annotation        | mb_n_strong | TCGA-LIHC dir (log2FC) | COSMIC CGC | ASoC adjP |
| ------------ | ------ | ----------------- | ----------- | ---------------------- | ---------- | --------- |
| rs2471847    | KLF10  | Promoter (≤1kb)   | 1           | Down (−1.44)           | —          | 4.0×10⁻¹⁰ |
| rs145317211  | XPR1   | Promoter (≤1kb)   | 7           | Up (+0.64)             | —          | 4.9×10⁻⁶  |
| rs10111451   | TPD52  | Promoter (≤1kb)   | 2           | Up (+0.45)             | —          | 8.3×10⁻⁶  |
| rs1050847260 | MAP3K2 | Promoter (≤1kb)   | 1           | Down (−0.51)           | —          | 1.2×10⁻⁵  |
| rs58065091   | NDRG1  | Distal Intergenic | 15          | Up (+1.19)             | ✔          | 1.4×10⁻⁴  |
| rs6444772    | HES1   | Distal Intergenic | 16          | Down (−0.98)           | —          | 1.1×10⁻³  |
| rs13034125   | MAP3K2 | Promoter (≤1kb)   | 3           | Down (−0.51)           | —          | 1.2×10⁻³  |
| rs9935289    | SPG7   | Promoter (≤1kb)   | 6           | Down (−0.28)           | —          | 2.0×10⁻³  |
| rs13034122   | MAP3K2 | Promoter (≤1kb)   | 1           | Down (−0.51)           | —          | 2.8×10⁻³  |
| rs140762989  | USP22  | Promoter (≤1kb)   | 6           | Up (+0.59)             | —          | 3.3×10⁻³  |
| rs1342691469 | MTR    | Promoter (≤1kb)   | 2           | Up (+0.76)             | —          | 4.1×10⁻³  |
| rs117165000  | MTR    | Promoter (≤1kb)   | 13          | Up (+0.76)             | —          | 4.3×10⁻³  |
| rs859002     | CD1D   | Promoter (≤1kb)   | 5           | Down (−1.59)           | —          | 8.0×10⁻³  |
| rs72759286   | NEK6   | Distal Intergenic | 2           | Down (−0.40)           | —          | 8.9×10⁻³  |
| rs2279339    | TFAM   | Promoter (≤1kb)   | 1           | Down (−0.38)           | —          | 1.3×10⁻²  |
| rs56165483   | YWHAZ  | Promoter (≤1kb)   | 5           | Up (+0.64)             | —          | 1.9×10⁻²  |
| rs2290193    | HDAC11 | Promoter (≤1kb)   | 3           | Up (+1.76)             | —          | 2.0×10⁻²  |
| rs555752     | CHEK1  | Promoter (≤1kb)   | 1           | Up (+2.22)             | —          | 2.4×10⁻²  |
| rs1428654600 | NDRG1  | Promoter (≤1kb)   | 11          | Up (+1.19)             | ✔          | 2.8×10⁻²  |
| rs567378743  | MOK    | Promoter (≤1kb)   | 1           | Up (+0.41)             | —          | 3.7×10⁻²  |
| rs73252164   | TPD52  | Promoter (≤1kb)   | 2           | Up (+0.45)             | —          | 4.0×10⁻²  |
| rs3818598    | FSD1L  | Promoter (≤1kb)   | 20          | Up (+1.17)             | —          | 4.1×10⁻²  |

> This 22-variant shortlist is a **manually reviewed** subset, not a purely rule-based cut (an automatic DisGeNET ∩ strong-TF filter yields 36 variant rows; 14 lower-priority/redundant rows were dropped by hand). **NDRG1** is the only member also in COSMIC CGC (two independent ASoC SNPs, rs58065091 and rs1428654600, both measured Up), keeping it the top convergent candidate. Most rows are promoter-proximal, and the panel is enriched for genes measured **Up** in TCGA-LIHC (CHEK1, HDAC11, NDRG1, FSD1L, USP22, YWHAZ, TPD52, MTR, XPR1, MOK).

### Takeaways

- **Curated overlap is modest but coherent.** Only 37/230 ASoC genes appear in a curated HCC set, and just **NDRG1** appears in both COSMIC CGC and DisGeNET — reinforcing that this ASoC dataset's value is *regulatory* (which hepatocyte genes are modulated), not a rediscovery of known HCC driver mutations.
- **NDRG1 is the top convergent candidate**: COSMIC CGC (HCC) + DisGeNET + measured Up in TCGA-LIHC + strong TF-motif disruption.
- **A manually-curated 22-variant shortlist** (Section 4c; DisGeNET C2239176 + strong TF disruption, hand-reviewed) across 17 genes is the recommended priority set for follow-up, led by NDRG1 (2 SNPs), CHEK1, HDAC11, USP22, and FSD1L.
- **DisGeNET breadth caveat.** The DisGeNET C2239176 set is text-mining-inclusive (3,593 genes), so its 33 hits reflect literature co-mention depth rather than curated causality; the Open Targets `cancer_gene_census` (5 genes) is the stricter, driver-oriented signal.
- **Provenance note.** Because COSMIC and DisGeNET are login/key-gated, faithful tokenless proxies were used (Open Targets `cancer_gene_census` datasource for HCC; Enrichr DisGeNET mirror). If Sanger COSMIC and a DisGeNET API key are available, re-running against the native `cancer_gene_census.csv` and the DisGeNET GDA `score`/`EI` fields would let hits be filtered by curated evidence level (and by CGC tier 1/2, germline/somatic, and role-in-cancer).

---

## Summary

1. **Variant level:** None of the `variantID` SNPs are curated disease/cancer alleles in dbSNP/ClinVar, and an **expanded search by exact rsID across Europe PMC literature and the GWAS Catalog also found no direct variant→HCC association**. They are regulatory (ASoC) variants; their relevance is mechanistic — they alter chromatin accessibility (and predicted TF motif binding, per `mb_disrupted_TFs`) at cancer-relevant loci.
   - The only variant with *any* published cancer-SNP literature is **rs9898809** at the **TOB1 / TOB1-AS1** tumor-suppressor locus (studied in **gastric cancer**).
   - Broadening to **liver-related (non-cancer) traits** surfaced one clear GWAS signal: **rs2060982 (CPT1A)** is genome-wide significant for **hepatic lipid / fatty-acid metabolism** quantitative traits (saturated-fatty-acid %, p≈3×10⁻³⁰; IDL triglycerides, p≈2×10⁻¹³) — metabolic/physiological, not a liver-disease diagnosis (see 1c).
   - The only variants with genome-wide GWAS signals are the **6p21 MHC cluster** (e.g. **rs16899941**, p≈3×10⁻⁵²), but those catalogued associations are **autoimmune (hypothyroidism, rheumatoid arthritis), not HCC**. A **formal Ensembl LD test in EUR and EAS** showed these SNPs are **not LD proxies** of the HCC MICA locus (rs2596542) — pairwise r² below threshold for all four in both populations — so the ~0.13 Mb proximity is incidental (extended-MHC gene density), not a shared HCC signal.
2. **Gene level:** Multiple regulated genes have established links to **hepatocellular carcinoma / liver biology**, most notably **HNF4G, NDRG1, TP53INP1, CPT1A, INSIG1, KLF10**, with secondary candidates **FZD7, MTSS1, CES1, USP22, TOB1, DUSP10, LAMC1, TLR5, CAPN2, GSDME**. Europe PMC co-mention counts confirm deep literature for several (NDRG1 ~1123, FZD7 ~755, USP22 ~699, TP53INP1 ~653, MTSS1 ~311, HNF4G ~187). Themes: hepatocyte-identity transcription factors (HNF4G), lipid-metabolism reprogramming (CPT1A, INSIG1, CES1), Wnt signaling (FZD7), p53 axis (TP53INP1), and metastasis suppression (NDRG1, MTSS1).
3. **Measured TCGA-LIHC differential expression (Section 3), strong TF-disruption panel (`mb_n_strong` > 0):** edgeR tumor-vs-normal DE (374 tumor / 50 normal; |log2FC|≥0.2, FDR<0.05) gives measured fold-changes. **Significantly UP in HCC:** CHEK1 (+2.22), HDAC11 (+1.76), LAMC1 (+1.39), NDRG1 (+1.19), CAPN2 (+0.76), GSDME (+0.70), YWHAZ (+0.64), USP22 (+0.59), FZD7 (+0.58), TPD52 (+0.45). **Significantly DOWN in HCC:** KLF10 (−1.44), DUSP10 (−1.14), HES1 (−0.98), MAP3K2 (−0.51), CPT1A (−0.48), NEK6 (−0.40), LRP5 (−0.29). **Not significant:** CES1, INSIG1, HNF4G, MTSS1, TP53INP1, GPAM. Most literature calls were confirmed, but two are **discordant** (measured ↓ despite literature ↑): **CPT1A** and **NEK6**; and TP53INP1/MTSS1 (literature ↓) were not significant in bulk TCGA-LIHC. Genome-wide, 111/34 ASoC genes are up/down. Genes removed by the strong-disruption filter: **TRIB1, TOB1, ING1, TLR5**.

**Recommended follow-up:**
- **LD-proxy analysis was run** (Ensembl LD REST API, `1000GENOMES:phase_3:EUR` and `:EAS`): the 6p21 MHC ASoC SNPs are **not** proxies of the HCC MICA lead rs2596542 in either population. If donors are AFR/AMR/SAS, re-run in the matched population and, optionally, also test the HLA-DQ/DP HCC leads (rs9272105/rs3077).
- **Confirm differential expression quantitatively:** intersect the ASoC `ENSEMBL` gene list with a TCGA-LIHC tumor-vs-normal DE table (GEPIA2 export or `TCGAbiolinks` in R; |log2FC| ≥ 1, FDR < 0.05) to attach measured fold-change/direction to each gene.
- **Curated HCC gene-set cross-reference was run** (Section 4): of 230 ASoC genes, **5** are in the HCC-mapped **COSMIC Cancer Gene Census** (NDRG1, PBX1, ELK4, PLEC, N4BP2), **33** are in **DisGeNET "Carcinoma, Hepatocellular" (C2239176)**, and only **NDRG1** is in both. 17 DisGeNET hits also carry strong TF-motif disruption. Tokenless proxies (Open Targets `cancer_gene_census` datasource; Enrichr DisGeNET mirror) were used since COSMIC/DisGeNET are gated; re-run natively if Sanger/DisGeNET credentials are available.
- Test whether ASoC-disrupted TF motifs (e.g. HNF4A/HNF4G, FOXA, RXR) are enriched — those are the hepatocyte lineage factors most relevant to liver-cancer regulatory rewiring.

---
*Generated on 2026-07-26 (Section 4 curated cross-reference added 2026-07-31). Variant-level status verified via dbSNP/ClinVar, Europe PMC (literature by exact rsID), and the NHGRI-EBI GWAS Catalog (by rsID) — searched for HCC, any liver-related disease/trait (NAFLD, cirrhosis, hepatitis, bilirubin, liver enzymes, lipid/fatty-acid metabolism), and gastric cancer. Section 3 tumor-vs-normal expression is **measured** from TCGA-LIHC (UCSC Xena GDC STAR counts; edgeR QL F-test, 374 tumor vs 50 normal; |log2FC|≥0.2, FDR<0.05); functional roles summarized from Europe PMC. LD-proxy hypothesis for the 6p21 MHC SNPs was formally tested and rejected via the Ensembl LD REST API (1000GENOMES:phase_3:EUR and :EAS). Section 4 cross-references the ASoC gene list against curated HCC gene sets — COSMIC Cancer Gene Census (HCC-mapped, via the Open Targets `cancer_gene_census` datasource for MONDO_0007256) and DisGeNET "Carcinoma, Hepatocellular" C2239176 (via the Enrichr DisGeNET mirror) — because the native COSMIC/DisGeNET downloads are login/API-key gated. This is a computational annotation summary, not a clinical determination.*
