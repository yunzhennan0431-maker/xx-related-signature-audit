# XX-Related Prognostic Signature Audit

Code and processed data for a study auditing whether "XX-related prognostic signature" papers
(candidate gene pool → clustering → LASSO-Cox risk model, keyed to a named biological process such
as ferroptosis, hypoxia, cuproptosis, angiogenesis, etc.) systematically encode tumor-microenvironment
cell-type composition rather than genuine pathway activity.

This repository accompanies the manuscript *"规模化审计一种诊断结论:基于LLM文献挖掘与跨主题单细胞归因,
审计'XX相关'预后基因签名中命名与机制的脱节现象"* (working title; manuscript in preparation).

## What's here

- `scripts/` — all analysis scripts (01–39), covering: literature search and download (PubMed
  E-utilities, PMC Cloud Service), candidate-table scanning and scoring, single-cell reference-atlas
  attribution (z-score method), cross-topic meta-analysis (exact binomial test + Benjamini-Hochberg
  FDR, permutation-based heterogeneity test), candidate-pool tautology screening, duplicate-publication
  screening (gene-set Jaccard + targeted coefficient-level checks), dataset lineage-coverage grading,
  and figure generation. Scripts are numbered roughly in the order they were developed/run; see
  in-file docstrings for what each one does and which manuscript section it supports.
- `data/` — processed/intermediate outputs: batch-by-batch manual review records
  (`verified_extractions_batch*.json`), meta-analysis result tables (`phase3_meta_analysis_*.csv`),
  duplicate/tautology screening outputs, dataset quality grading, and the reference table of the
  108 models included in the quantitative analysis set (`reference_table_allmodels.csv`/`.md`,
  PMCID + title + cancer type + theme + gene count only — no reproduction of the underlying papers'
  content).
- `figures/` — the manuscript's generated figures (workflow diagram, trajectory plots, H1a/H1b
  distribution comparison), PNG + PDF.

## What's deliberately NOT here

- Raw single-cell reference expression matrices: these are public data from GEO and the TISCH2
  database; see the manuscript's Methods section for the full list of accessions (GSE81861,
  GSE176078, GSE151530, etc.) and fetch them directly from GEO / http://tisch.compbio.cn.
- The 108 audited papers' own supplementary files (tables, figures, raw data) — these remain the
  copyright of their original publishers; only their PMCID and extracted, minimal factual data
  (gene symbols, coefficients, cancer type) needed for this audit are reproduced here, consistent
  with standard scholarly data-extraction practice for a meta-analysis.
- Internal project logs and planning documents.

## Reproducing the analysis

Scripts assume local copies of the public single-cell datasets listed above, laid out under paths
configured near the top of each script (adjust `DATA_DIR`/dataset-specific paths to your own layout).
Core dependencies: `pandas`, `numpy`, `scipy`, `h5py`, `matplotlib`. No GPU or deep-learning framework
is required — the cell-type attribution step is a classical z-score pipeline (log-normalized
expression → mean by annotated cell type → z-score across cell types → argmax), chosen after an
initial single-cell foundation model (SCimilarity) evaluation found its practical interface
incompatible with this study's bare-gene-list query pattern (see manuscript §2.3).

## Status

This is a working repository accompanying a manuscript still in preparation — expect the script
set and data files to be updated as the manuscript is revised. Questions/issues welcome via GitHub
Issues.
