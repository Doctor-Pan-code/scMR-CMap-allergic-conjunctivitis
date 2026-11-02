# scMR-CMap: causal target discovery and compound reversal in allergic conjunctivitis

Analysis code and result tables for the manuscript

> **Causal target and reversal-compound screening in allergic conjunctivitis: a single-cell-informed Mendelian randomization and connectivity-map study**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23890289.svg)](https://doi.org/10.5281/zenodo.23890289)
[![Licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](requirements.txt)

This repository contains every self-written script used to produce the manuscript: the druggable-genome pool construction, the genome-wide cis-eQTL Mendelian randomization screen, colocalization with SMR/HEIDI, the lineage-restricted single-cell expression-state contrast in human bulbar conjunctiva, the L1000 connectivity-map reversal query with null calibration, the phenome-wide safety sweep, and the five-axis evidence integration that yields the final candidate ranking. It also contains the figure and manuscript build scripts, so the six published figures and the submission `.docx` regenerate from the tables in this repository without any hand-typed numbers.

Everything is plain Python 3 with a deliberately small dependency set — the statistics engine (`scripts/lib_stats.py`) implements IVW, MR-Egger, weighted median, weighted/simple mode, Cochran's Q, leave-one-out, the F statistic, the Steiger test, `coloc.abf` and a HEIDI approximation in pure Python, so no R installation is required.

---

## Key results reproduced here

| Axis | Result |
| --- | --- |
| Druggable pool | 19,255 protein-coding genes × DGIdb × OpenGWAS cis-eQTL → 10,983 druggable genes in three pharmacological tiers; 8,009 instrumentable |
| MR screen | 6,329 genes, median first-stage F = 347.7; 494 at P < 0.05, 3 at FDR < 0.05, 2 Bonferroni-significant |
| Colocalization | 47 genes tested, 6 with PP.H4 ≥ 0.80; TLR1 PP.H4 = 0.970, SMR P = 9.38×10⁻⁷, HEIDI P = 0.092 |
| Single cell | TLR1 detected in 11.6–20.6 % of cells across eight conjunctival lineages; 600-gene state signature; mean cross-lineage directional concordance 73.5 % vs. +0.010 in the stratified null |
| Connectivity map | 44/600 signature genes are L1000 landmarks; observed score 0.250 vs. null 0.155 ± 0.039 at 15 genes per arm (z = 2.47, permutation P = 0.0476); BI-2536 top reversal specificity 0.0949 |
| Safety | 14 genes × 1,400 FinnGen R13 endpoints; 106 significant associations; 10 safe, 3 monitor, 1 high-risk (PBXIP1 withdrawn) |
| Integration | TLR1 composite score 0.7631, the only Tier 1 candidate |

---

## Repository layout

```
config/config.json                 all analysis thresholds, endpoints, weights, seeds
config/opengwas_token.txt.example  template for the OpenGWAS JWT (not distributed)
data/l1000_landmark_genes.txt      L1000 landmark gene list used by the CMap query
data/README.md                     how to fetch the four public input files
scripts/day1_druggable_pool.py     druggable genome universe and pharmacological tiers
scripts/day2_mr_screen.py          cis-eQTL MR screen (Wald ratio) + confirmatory pass
scripts/day3_coloc_smr.py          coloc.abf, SMR Wald ratio, HEIDI heterogeneity
scripts/day4_scrna_signature.py    lineage-restricted expression-state signature (v7)
scripts/day4_probe_genes.py        per-cell probe of candidate genes across lineages
scripts/day5_cmap_reversal.py      L1000CDS2 reversal query, null calibration, Enrichr cross-check
scripts/day6_phewas_safety.py      FinnGen R13 phenome sweep and three-tier safety verdict
scripts/day7_integration.py        five-axis scoring, tiering, figure assembly
scripts/lib_io.py                  paths, config loading, CSV/JSON IO, logging, SVG/PNG saving
scripts/lib_api.py                 cached HTTP clients (OpenGWAS, Ensembl, L1000CDS2, Enrichr, DGIdb)
scripts/lib_stats.py               pure-python statistics engine (no scipy, no R)
scripts/run_pipeline.py            sequential runner for Days 3-7
scripts/report/make_figs.py        the six published figures, read straight from tables/
scripts/report/build_manuscript.py submission .docx builder (template + normalised text)
tables/*.csv                       every result table cited in the manuscript
figures/fig1..fig6 .png/.pdf       600 dpi PNG and vector PDF of the six figures
manuscript_src/                    normalised manuscript text, reference JSON, SR Word template
docs/ANALYSIS_MAP.md               method section -> script -> table -> figure cross-reference
docs/DATA_DICTIONARY.md            column-level description of every output table
docs/gen_data_dictionary.py        regenerates DATA_DICTIONARY.md from tables/*.csv
CITATION.cff                       machine-readable citation metadata (concept DOI 10.5281/zenodo.23890289)
```

---

## Requirements

Python 3.10 or newer (developed on 3.14).

```
pip install -r requirements.txt
```

`numpy`, `pandas`, `matplotlib`, `requests`, `python-docx` and `Pillow` are the only direct dependencies. `scipy`, `statsmodels` and R are **not** required — the statistical estimators are implemented in `scripts/lib_stats.py` and validated against the published formulas.

---

## Reproducing the analysis

### 1. Fetch the public inputs

The four upstream files are not redistributed here (they total ~55 MB and each has its own licence and update cycle). Download them into `data/` as described in [`data/README.md`](data/README.md):

```
data/hgnc_complete_set.txt      HGNC complete set
data/dgidb_genes.tsv            DGIdb gene claims
data/dgidb_interactions.tsv     DGIdb drug-gene interactions
data/opengwas_catalog.json      OpenGWAS dataset catalogue
data/geo/GSE242382/             GSE242382 matrix.mtx.gz, features, barcodes, clusters
```

### 2. Provide an OpenGWAS token

Register at <https://api.opengwas.io> and write the JWT to `config/opengwas_token.txt` (single line, no quotes):

```
cp config/opengwas_token.txt.example config/opengwas_token.txt
# then paste your token
```

The token is read by `lib_io.token()`. It is deliberately excluded by `.gitignore` and is **not** part of this repository.

### 3. Run the pipeline

```bash
python scripts/day1_druggable_pool.py      # -> tables/day1_*.csv
python scripts/day2_mr_screen.py           # -> tables/day2_*.csv   (long: genome-wide API sweep)
python scripts/run_pipeline.py             # Day 3-7 in order
```

`run_pipeline.py` runs Day 3 through Day 7 sequentially, writes one log per step to `logs/`, and continues past a failing step so a single pass shows every independent failure. HTTP responses are cached under `data/intermediate/cache/`, so an interrupted run resumes without re-querying.

Runtime on the original machine: Day 2 is the long pole (a genome-wide sweep over 6,329 genes against the OpenGWAS API at ~1.2 s per request); Days 3–7 complete in roughly an hour once the cache is warm.

### 4. Regenerate the figures and the manuscript

```bash
python scripts/report/make_figs.py         # -> figures/fig1..fig6 .png + .pdf
python scripts/report/build_manuscript.py  # -> manuscript_out/_build.docx
```

`make_figs.py` reads `tables/*.csv` at run time and prints a per-figure audit line (canvas size in mm, type sizes, out-of-canvas elements, text collisions). No number in the figures is typed by hand. `build_manuscript.py` needs `manuscript_src/tpl.docx` (the Scientific Reports Word template), `manuscript_src/_norm.txt` (run-level normalised text) and `manuscript_src/_refs_nature.json` (formatted reference list), all of which are included.

---

## Analysis design in one paragraph

The study orders its evidence causally before pharmacologically. Genetic instruments are restricted to **cis** regulation of expression, so the exposure sits proximal to the causal chain; the pool is restricted to genes with existing pharmacological annotation, so a genetic association reads as drug-target evidence rather than as a generic biomarker association. The strongest candidate is then projected into a tissue single-cell atlas, where the contrast of interest is deliberately **not** a cell-type comparison: the target gene is detected in only 11.6–20.6 % of cells in every lineage, so splitting within lineage by residualised expression yields a *state* contrast that reproduces across six lineages (73.5 % directional concordance) and collapses to r = +0.010 under a randomly permuted stratifier. That signature is what enters the connectivity map, calibrated on both ends — a dexamethasone positive control that must be recovered, and a matched permutation null at each query size. Safety is a fallback filter with a veto, not a score: a candidate that moves severe systemic domains across 1,400 FinnGen endpoints leaves the list regardless of how strong its genetic evidence is.

## Thresholds

Every threshold, endpoint identifier, weight and random seed lives in [`config/config.json`](config/config.json) and is read at run time — cis window size, PP.H4 cut-offs, the single-cell QC and stratification parameters, the CMap query sizes and null draws, the PheWAS P-value and domain keywords, and the five-axis weights (MR 0.30, coloc 0.20, single cell 0.20, CMap 0.20, safety 0.10). No threshold is hard-coded in an analysis script.

## Data availability

All inputs are public: FinnGen release 13 (allergic conjunctivitis endpoint `finn-b-H7_ALLERGICCONJUNCTIVITIS`), eQTLGen whole-blood cis-eQTL through IEU OpenGWAS, the human bulbar conjunctiva single-cell atlas **GSE242382**, HGNC, DGIdb and LINCS L1000. No new primary data were generated.

## Citation

If this code contributes to your work, please cite the archived release. The concept DOI always resolves to the most recent version:

> Pan Z. *scMR-CMap-allergic-conjunctivitis: reproducible code and result tables for a causal-target and reversal-compound screen in allergic conjunctivitis*. Zenodo. https://doi.org/10.5281/zenodo.23890289

Version-specific DOIs: `v1.0.0` → https://doi.org/10.5281/zenodo.23890290; `v1.0.1` → https://doi.org/10.5281/zenodo.23891048. Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

## Licence

MIT — see [`LICENSE`](LICENSE).
