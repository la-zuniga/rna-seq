A Nextflow DSL2 pipeline for single-cell RNA-seq analysis, from raw FASTQ files through clustering, differential expression, and pathway analysis. Built around STARsolo for alignment/quantification and Scanpy for downstream single-cell processing.

**Pipeline Overview:**

The pipeline chains together the following steps:

STAR Index — Builds a STAR genome index from a reference FASTA and GTF annotation.

FastQC — Per-sample read quality assessment.

MultiQC — Aggregated quality report across all samples.

STARsolo — Alignment and gene expression quantification for single-cell data.

Merge Counts — Combines per-sample STARsolo count matrices into a single object.

Scanpy QC — Cell-level quality control and filtering

Scanpy Preprocess — Normalization, log-transformation, and highly variable gene selection.

Scanpy Integrate — Batch correction / data integration across samples.

Scanpy Cluster — Dimensionality reduction (PCA, UMAP) and cell clustering with annotation.

Pseudobulk DE (optional) — Pseudobulk differential expression analysis.

Pathway Analysis (optional) — Functional enrichment on DE results.


**Requirements:**

Nextflow (DSL2-compatible, ≥ 21.04)
STAR
FastQC
MultiQC
Python 3 with Scanpy

**Output:**

Results are written to the directory specified by --outdir and include:

FastQC / MultiQC reports — HTML quality summaries
STARsolo outputs — BAM files, count matrices, and alignment logs
Scanpy objects — Filtered, preprocessed, and annotated .h5ad files at each stage
DE results — Pseudobulk differential expression tables (when --run_de is enabled)
Pathway analysis — Enrichment results from DE gene sets (when --run_de is enabled)
