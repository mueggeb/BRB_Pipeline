# BRB-seq Pipeline — CLAUDE.md

## Project Purpose
This pipeline processes raw BRB-seq RNA-seq data on a SLURM cluster and produces
a gene count matrix plus QC reports. BRB-seq is a bulk RNA-seq method using
sample barcodes for multiplexing.

## Pipeline Architecture
Two-stage design — preserve this structure during refactoring:

**Stage 1 (per-sample, SLURM array):** scripts/BRB_Seq_Muegge_241112_TwoStage_noMulti.sh
  - Optional demultiplexing by RT barcode (cutadapt)
  - FastQC on raw reads
  - Adapter trimming (cutadapt)
  - PolyA trimming (cutadapt)
  - FastQC on trimmed reads
  - STAR alignment
  - featureCounts → per-sample count matrix
  - RSeQC + Qualimap QC
  - Optional cleanup of intermediates

**Stage 2 (aggregation, runs after all samples complete):** scripts/BRB_MultiQC_FullRun_240116.sh
  - Merge per-sample count matrices → single count matrix
  - Merge per-sample log files → single summary table
  - MultiQC report generation
  - Cleanup

## Reference Genomes
Mouse and human both in use. Reference paths (STAR index, GTF, BED) are passed
via config YAML. A future goal is selecting genome by name (e.g. "mouse", "human")
rather than requiring manual path entry.

## Environment
- Cluster: HTCF (Washington University)
- Current dependency manager: SPACK (fragile, breaks after admin updates)
- Target dependency manager: Mamba with pinned environments
- Python target version: 3.10+
- All tools used are standard bioinformatics packages available via conda-forge
  or bioconda. If a tool substitution or addition would be a clear improvement,
  explain the rationale before implementing.

## Refactoring Goals (in order)
1. Translate bash logic to Python — no new features
2. Replace SPACK loads with mamba environment
3. Consolidate mapping file + parameters CSV into single YAML config
4. Add reference genome selection by name or build (mouse/human, or GENCODE/ENSEMBL)
5. Improve error handling and logging

## Out of Scope (future work)
- Merging Stage 1 and Stage 2 into a single orchestrated pipeline
- Automatic resubmission of failed SLURM jobs

## Coding Conventions
- Lab members are biologists, not software engineers — prioritize readable code
- Every pipeline step should be a named Python function with a clear docstring
- Log key metrics at each step (match existing log file format where possible)
- All paths and parameters via YAML config — no hardcoded paths, no argparse
- Prefer explicit variable names over terse abbreviations

## Key Hardcoded Paths to Parameterize
- /ref/rmlab/software/spack-0.22/ (fastqc spack path)
- /ref/bmlab/software/BRB-Seq/ (multiqc config template location)

## Inputs
- Raw FASTQ files (Read 1: barcodes, Read 2: cDNA)
- Mapping file: TSV with columns SampleName, Group, RT Barcode, Index1, Index2
- YAML config: STAR index path, GTF path, BED path, demultiplex flag, remove flag

## Outputs
- Per-sample featureCount matrices → merged count matrix (TSV)
- Per-sample log files → merged summary table
- MultiQC HTML report
- Intermediate files optionally removed

## Do Not
- Change pipeline functionality during the refactor phase
- Break the two-stage structure
- Remove existing log metrics (needed for QC comparisons across runs)
- Introduce tools not available on HTCF without explaining why first
