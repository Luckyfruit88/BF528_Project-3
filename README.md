# BF528 Project 3: ChIP-seq Genomic Data Analysis Pipeline

This repository contains a modular, containerized Nextflow DSL2 pipeline designed for the comprehensive analysis of ChIP-seq data, as developed for the Boston University Bioinformatics course [BF528 - Genomic Data Analysis](https://bu-bioinfo.github.io/bf528/projects/project_3_chipseq/).

## Project Overview

The goal of this project is to process and analyze a ChIP-seq experiment from human cell lines. The dataset includes two biological replicates, each consisting of paired IP (factor of interest) and INPUT (background control) samples. The pipeline takes raw sequencing reads (FASTQ) and performs a complete analysis to identify regions of enriched binding, annotate these regions to genomic features, and prepare data for downstream integration with RNA-seq differential expression results.

The analysis is divided into three core computational phases (Weeks 1-3) implemented in Nextflow, followed by a data interpretation and reproducibility phase (Week 4) conducted in a Jupyter Notebook.

## Pipeline Architecture

This workflow is built using Nextflow DSL2 and strictly follows the "One Process per Module" design pattern. Each tool runs within its own isolated computational container (e.g., Docker/Singularity/Conda) and requests specific cluster resources dynamically via profiles.

### Directory Structure

```text
project_3_chipseq/
├── main.nf                              # Main workflow script with advanced channel routing (branch/join)
├── nextflow.config                      # Configuration for params, SGE cluster profiles, and absolute paths
├── full_samplesheet.csv                 # CSV driving the workflow with full datasets
├── subsampled_samplesheet.csv           # CSV with subset data for rapid testing and development
├── environment.yml                      # Conda environment for downstream Jupyter Notebook analysis
├── project3_downstream_analysis.ipynb   # Week 4 Notebook: visualization, motif enrichment, and paper reproduction
│
├── modules/                             # Independent Nextflow modules for each tool
│   ├── fastqc/                          # Raw read quality control
│   ├── trimmomatic/                     # Adapter trimming and low-quality base removal
│   ├── bowtie2_build/                   # Reference genome index generation
│   ├── bowtie2_align/                   # Read alignment to the reference genome
│   ├── samtools_sort/                   # BAM sorting
│   ├── samtools_idx/                    # BAM indexing
│   ├── samtools_flagstat/               # Alignment statistics generation
│   ├── multiqc/                         # Aggregation of all QC metrics
│   ├── deeptools_bamcoverage/           # Generation of bigWig coverage tracks
│   ├── deeptools_multibwsummary/        # Correlation matrix computation across bigWigs
│   ├── deeptools_plotcorrelation/       # Sample correlation visualization (Pearson/Spearman)
│   ├── homer_maketagdir/                # HOMER tag directory creation
│   ├── homer_findpeaks/                 # Peak calling (pairing IP & INPUT per replicate)
│   ├── homer_pos2bed/                   # Conversion of HOMER peaks to standard BED format
│   ├── bedtools_intersect/              # Generation of reproducible consensus peaks
│   ├── bedtools_remove/                 # Filtering out ENCODE blacklist regions
│   ├── homer_annotatepeaks/             # Annotation of peaks to nearest genomic features
│   ├── deeptools_computematrix/         # Signal computation across scaled gene regions (IP only)
│   ├── deeptools_plotprofile/           # Metagene signal profiling
│   └── homer_findmotifsgenome/          # Motif enrichment analysis in filtered peaks
│
├── data/                                # Directory for user-downloaded data
│   └── rnaseq_results.csv               # Provided RNA-seq DE results for cross-omics validation
│
└── results/                             # Automatically generated output directory
    ├── qc/                              # MultiQC reports and individual logs
    ├── alignments/                      # Processed BAMs and indices
    ├── bigwigs/                         # Genome browser tracks
    ├── final_peaks/                     # High-quality, reproducible BED files
    └── annotations/                     # Genomic feature annotations
```

## Workflow Execution Steps

### Week 1: Pre-processing and Alignment
- **QC & Trimming**: Raw reads are evaluated with `FastQC` and cleaned using `Trimmomatic`.
- **Alignment**: Reads are aligned to the human reference genome (hg38) using `Bowtie2` (non-splice aware, suitable for ChIP-seq).
- **BAM Processing**: Alignments are sorted and indexed using `Samtools`, and `flagstat` calculates alignment rates.
- **Coverage Tracks**: `deepTools bamCoverage` converts BAMs to normalized bigWig files for visualization.
- **QC Aggregation**: `MultiQC` compiles all FastQC, Trimmomatic, and flagstat logs into a single HTML report.

### Week 2: Peak Calling and Annotation
- **Sample Correlation**: `deepTools multiBigwigSummary` and `plotCorrelation` evaluate the similarity between biological replicates vs. controls.
- **Peak Calling**: `HOMER` (makeTagDirectory & findPeaks) is used to identify enriched regions. *Crucially, the pipeline dynamically parses sample IDs to ensure the IP sample is correctly paired with its matching INPUT control for each replicate.*
- **Consensus & Filtering**: `bedtools intersect` extracts reproducible peaks found in both replicates. `bedtools intersect -v` is then used to rigorously remove artifactual signals overlapping with the ENCODE Blacklist.
- **Annotation**: `HOMER annotatePeaks` maps the filtered peaks to their nearest gene elements (Promoters, Introns, Exons, etc.) using the hg38 GTF.

### Week 3: Signal Profiling and Motif Discovery
- **Signal Profiling**: Using a curated BED file of hg38 genes, `deepTools computeMatrix` and `plotProfile` generate a metagene plot showing the average IP signal distribution around the Transcription Start Sites (TSS) and Transcription Termination Sites (TTS).
- **Motif Enrichment**: `HOMER findMotifsGenome` scans the reproducible peak sequences to discover enriched DNA-binding motifs, hinting at potential transcription factor complexes.

### Week 4: Downstream Interpretation (Notebook)
All downstream analysis is performed within `project3_downstream_analysis.ipynb`. This phase focuses on reproducibility and cross-omics integration:
1. Re-creating genomic tracks of key genes (Figures 2D/2E of the original paper).
2. Integrating the Nextflow pipeline's ChIP-seq peaks with external RNA-seq differential expression data to link binding events with transcriptional changes (Figure 2F).
3. Analyzing discrepancies in mapped reads, peak counts, and correlation metrics compared to the original authors' supplementary data.
4. Performing biological pathway enrichment (e.g., GREAT/DAVID) on the annotated peaks.

## Usage Instructions

### 1. Setup Environment
Ensure you are logged into the BU SCC. The pipeline utilizes the Sun Grid Engine (SGE) and requires no local software installation other than Nextflow, as all tools are containerized.

### 2. Run the Workflow
To execute the pipeline using the subset data for rapid testing:
```bash
nextflow run main.nf -profile cluster --samplesheet subsampled_samplesheet.csv -resume
```

To run the full analysis:
```bash
nextflow run main.nf -profile cluster --samplesheet full_samplesheet.csv -resume
```

### 3. Launch the Jupyter Notebook (Week 4)
Build the required computational environment and launch Jupyter:
```bash
conda env create -f environment.yml
conda activate bf528_project3_w4
jupyter notebook
```

## Configuration Details
The `nextflow.config` file relies on BU SCC absolute paths for reference genomes to save storage and computation time:
- Adapter sequences (`TruSeq3-SE.fa`)
- Human Reference Genome (`GRCh38.primary_assembly.genome.fa`)
- GENCODE Annotations (`gencode.v45.primary_assembly.annotation.gtf`)
- ENCODE Blacklist v2 (`hg38-blacklist.v2.bed`)

*Note: If running outside of the BU SCC environment, these paths must be updated in `nextflow.config`.*
