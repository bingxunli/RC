# Canonical numbers

Every number that appears in the manuscript, where it comes from, and how to
recompute it. **If a number here disagrees with the manuscript, this file is
right** (it is generated from the actual result files).

Last reconciled: 2026-10-19, after the PoPS / MAGMA / FLAMES re-run.

---

## Meta-analysis

| Quantity | Value | Source |
|---|---|---|
| Total N | **570,662** | 136,016 (Karjalainen) + 434,646 (Zoodsma) |
| Variants tested | **11,071,938** | rows in `magma/RC_magma.sumstats` |
| Genome-wide significant (P < 5×10⁻⁸) | **19,414** | `RC.R` clumping block |
| Independent lead SNPs | **125** | `RC.R` |
| FUMA loci | **112** | `FUMA_job714838/GenomicRiskLoci.txt` |
| Novel loci | **13** | `RC.R` novelty block |
| Phenotypic variance explained | ~7.7 % | `RC.R` |

---

## Fine-mapping

| Quantity | Value | Source |
|---|---|---|
| PIP threshold | 0.99 | `finemap.R` |
| High-confidence variants | **497** | `finemap/` |
| Fraction at PIP = 1.00 | **96.2 %** | `finemap/` |
| Loci with a credible set | **110** of 112 | `cred1/` |
| Loci not fine-mapped | 32 (chr4:71.8–74.2 Mb), 101 (CYP2A6) | 1kg v3 EUR panel re-download pending |

> **497 vs 506.** An earlier draft carried 506 in one sentence (Discussion).
> 497 is correct: it is the number defined in the Results paragraph
> ("retaining a high-confidence set of 497 variants") and repeated in the
> abstract. The 506 figure came from the raw `high_pip` row count before the
> duplicate-SNP filter. Unified to **497** on 2026-10-19.

---

## MAGMA

| Quantity | Value | Source |
|---|---|---|
| Genes tested | **18,952** | `magma/magma.genes.out` |
| Bonferroni threshold | **2.638×10⁻⁶** | 0.05 / 18,952 |
| Significant genes | **400** (all positive Z) | same |
| Top gene | **ALDH1A2**, Z = 12.00, P = 1.76×10⁻³³ | same |
| GTEx v8 tissues tested | **54** | `magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out` |
| Tissues passing Bonferroni (9.26×10⁻⁴) | **7** | same |
| Top tissue | **Liver**, β = 0.0654, P = 1.91×10⁻¹⁶ | same |
| Gene sets tested | **17,005** | `magma.gsa.out` |
| Gene sets passing (2.94×10⁻⁶) | **63** | same |
| Top gene set | WP_CHOLESTEROL_METABOLISM, P = 3.96×10⁻¹⁹ | same |
| Statin pathway rank | **4th**, P = 1.02×10⁻¹⁶ | WP_STATIN_INHIBITION_OF_CHOLESTEROL_PRODUCTION |

---

## PoPS

| Quantity | Value | Source |
|---|---|---|
| Genes scored | **18,383** | `pops/run_full/RC_pops_full.preds` |
| Features selected | **13,768** | `.coefs` (published run: 13,608) |
| RidgeCV α | **31622.7766** | identical to the published run |
| Correlation with published | Pearson **0.99717**, Spearman 0.99595 | `pops_vs_published.csv` |
| top-100 overlap | **97/100** | same |
| top-500 overlap | **479/500** | same |
| PoPS top gene | **ALDH1A2** (4.856) | same |

> `UBC` (ubiquitin C) ranks 7th in every PoPS run and is a housekeeping
> artefact — remove it from any table shown in the paper.

---

## FLAMES

| Quantity | Value | Source |
|---|---|---|
| Loci annotated | **110** | `annots/` |
| Genes annotated | **2,185** | `flames/rerun/FLAMES_scores_newpops.raw` |
| Nominations (`FLAMES_causal = 1`) | **74** | same |
| Unique genes | **73** | same (`MAFB` nominated at loci 109 and 110) |
| Loci with ≥1 nomination | **73** | same |
| Median estimated cumulative precision | **0.854** | same |
| Top nomination | **TRPS1** (scaled = 1.00) | single-gene credible set — structural, not strength of evidence |

> The earlier 75 came from merging two separate runs (62 genes / 89 loci +
> 13 genes / 21 rebuilt loci). The current 74/73 comes from **one** re-scored
> run over all 110 loci, which is the cleaner provenance and is what the
> manuscript now reports. The single change relative to the merged list is that
> **ATXN2** drops out (0.262 → 0.246, just under the 0.248 gate).

### Why the FLAMES list is a high-precision subset, not the full set

`FLAMES_scaled = raw / sum(raw)` normalises within locus, so gene-dense loci are
penalised:

| Statistic | Value |
|---|---|
| Spearman(genes per locus, max scaled) | **−0.738**, P = 3.9×10⁻²⁰ |
| Median genes/locus, nominated | **12** |
| Median genes/locus, silent | **31** (Mann–Whitney P = 2.4×10⁻⁹) |
| Loci passing `raw > 0.134` but failing `scaled > 0.248` | **23** |
| Nominations if the `scaled` gate is dropped | **74 → 97** |

The loci lost this way are the textbook ones: APOA1/C3/A4/A5 (locus 73),
SORT1/CELSR2/PSRC1 (locus 8), GCKR (locus 15), FADS1/2 (locus 71).

---

## Cross-method convergence

| Set | n |
|---|---|
| MAGMA significant | 400 |
| PoPS top-100 | 100 |
| FLAMES nominated | 73 |
| **MAGMA ∩ PoPS ∩ FLAMES** | **29** |
| MAGMA ∩ FLAMES, PoPS outside top-100 | **23** |
| MAGMA ∩ PoPS, not nominated | 52 |
| PoPS ∩ FLAMES, not MAGMA-significant | 31 |
| FLAMES only | 19 |

FLAMES agrees with MAGMA far more than expected: **48/73 (66 %)** of its picks are
the top MAGMA gene in the same window, **69/73 (95 %)** are in the top three.
Where it deviates, it follows PoPS (Wilcoxon P = 1.74×10⁻⁴ on paired PoPS ranks).

---

## SMR / coloc

| Quantity | Value |
|---|---|
| SMR gene–RC associations | **36** |
| Colocalized | **21** |
| Unique colocalized genes | **18** |
| Tissues | 49 GTEx v8 + eQTLGen whole blood + deCODE plasma pQTL |

---

## PheWAS-MR

| Quantity | Value |
|---|---|
| FinnGen traits tested | 1,271 |
| Traits with a causal effect | **21** |
| Strongest | hypercholesterolaemia OR 1.65; statin use OR 1.76; MI OR 1.22 |

---

## Regenerating this table

```bash
python code/09_analysis/interpret_three_methods.py
python code/09_analysis/locus_arbitration.py
```

`results/three_methods_interpret.txt` prints every MAGMA / PoPS / FLAMES count
used above; `results/locus_arbitration.csv` has the per-locus breakdown.
