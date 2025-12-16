# Nanopore Direct RNA m6A Detection Pipeline

This repository contains the data processing pipeline for detecting N6-methyladenosine (m6A) modifications in *Arabidopsis thaliana* using Oxford Nanopore direct RNA sequencing.

## Overview

This pipeline processes raw nanopore sequencing data (POD5 format) to identify m6A modification sites across different developmental stages of *Arabidopsis thaliana*.

### Experimental Design

- **Species**: *Arabidopsis thaliana* (Col-PEK1.5)
- **Developmental stages**: Day 9, Day 21, Day 35
- **Replicates**: 3 biological replicates per stage
- **Sequencing**: Oxford Nanopore direct RNA sequencing (R10.4.1 chemistry)
- **Barcoding**: SQK-RBK114-24 kit

## Pipeline Workflow

```
POD5 files
    │
    ▼
┌─────────────────────────────────────┐
│  Step 1: Basecalling (GPU)          │
│  - Dorado basecaller                │
│  - SUP model + m6A modification     │
│  - Barcode classification           │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Step 2: Demultiplexing             │
│  - Separate reads by barcode        │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Step 3: Merge Technical Replicates │
│  - Combine runs for same sample     │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Step 4: Alignment (CPU)            │
│  - minimap2 to reference genome     │
│  - Preserve modification tags       │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  Step 5: m6A Site Extraction        │
│  - modbam2bed                       │
│  - Filter high-confidence sites     │
└─────────────────────────────────────┘
    │
    ▼
m6A sites (BED format)
```

## Requirements

### Software Dependencies

| Software | Version | Purpose |
|----------|---------|---------|
| [Dorado](https://github.com/nanoporetech/dorado) | ≥1.1.1 | Basecalling with modification detection |
| [samtools](http://www.htslib.org/) | ≥1.17 | BAM file manipulation |
| [minimap2](https://github.com/lh3/minimap2) | ≥2.26 | Read alignment |
| [modbam2bed](https://github.com/epi2me-labs/modbam2bed) | ≥0.9.0 | m6A site extraction |

### Dorado Models

- Basecalling: `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`
- m6A detection: `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1`

### Hardware Requirements

- **GPU**: NVIDIA GPU with ≥16GB VRAM (H100 recommended for optimal speed)
- **CPU**: ≥16 cores for alignment
- **RAM**: ≥64GB

## Installation

```bash
# Clone this repository
git clone https://github.com/JohnnyChen1113/nanopore-m6a-pipeline.git
cd nanopore-m6a-pipeline

# Install dependencies via conda (recommended)
conda create -n nanopore_m6a python=3.10
conda activate nanopore_m6a
conda install -c bioconda samtools minimap2 modbam2bed

# Download Dorado from ONT GitHub releases
# https://github.com/nanoporetech/dorado/releases

# Download models
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1
```

## Usage

### Step 1: Basecalling with m6A Detection

```bash
# For SLURM clusters with GPU
sbatch scripts/01_basecalling.sh \
    -i raw_data/pod5/ \
    -o results/run1/ \
    -b SQK-RBK114-24 \
    -d /path/to/dorado \
    -m /path/to/models/

# Or run directly
dorado basecaller \
    dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
    raw_data/pod5/ \
    --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1 \
    --kit-name SQK-RBK114-24 \
    > results/basecalled.bam

dorado demux \
    --output-dir results/demultiplexed/ \
    --no-classify \
    results/basecalled.bam
```

### Step 2: Merge Samples from Multiple Runs

```bash
# Edit config/sample_config.tsv with your barcode-to-sample mapping
# Then run:
bash scripts/03_merge_samples.sh \
    -c config/sample_config.tsv \
    -o results/merged/
```

### Step 3: Alignment and m6A Extraction

```bash
# For each sample
sbatch scripts/02_alignment_m6a.sh \
    -i results/merged/day9_rep1.bam \
    -o results/m6a/day9_rep1 \
    -r reference/genome.fasta \
    -t 16 \
    -f 0.66

# Or run directly
samtools bam2fq -T MM,ML input.bam > reads.fastq
minimap2 -ax map-ont -y reference.fasta reads.fastq | samtools sort -o aligned.bam
samtools index aligned.bam
modbam2bed -e -m 6mA -t 16 -f 0.66 reference.fasta aligned.bam > m6a_sites.bed
```

## Output Files

| File | Description |
|------|-------------|
| `*_aligned.bam` | Aligned reads with modification tags |
| `*_m6a_sites.bed` | All detected m6A sites |
| `*_m6a_high_conf.bed` | High-confidence sites (coverage ≥5, mod rate ≥20%) |
| `*_alignment_stats.txt` | Alignment statistics |

### BED File Format

The m6A BED files contain the following columns:

| Column | Description |
|--------|-------------|
| 1 | Chromosome |
| 2 | Start position (0-based) |
| 3 | End position |
| 4 | Modified base |
| 5 | Score |
| 6 | Strand |
| 7-9 | BED12 fields |
| 10 | Coverage (read depth) |
| 11 | Modification percentage |

## Parameters

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| m6A threshold | 0.66 | Minimum probability for calling m6A |
| Min coverage | 5 | Minimum read depth for high-confidence sites |
| Min mod rate | 20% | Minimum modification rate for high-confidence sites |

## Directory Structure

```
.
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_basecalling.sh      # GPU basecalling with m6A detection
│   ├── 02_alignment_m6a.sh    # Alignment and m6A extraction
│   └── 03_merge_samples.sh    # Merge technical replicates
├── config/
│   └── sample_config.tsv.example
└── docs/
    └── workflow.md
```

## Citation

If you use this pipeline, please cite:

```
[Your paper citation here]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- Author: Junhao Chen
- Email: junhao.chen.1@slu.edu
- Institution: Saint Louis University

## Acknowledgments

- Oxford Nanopore Technologies for Dorado
- The developers of samtools, minimap2, and modbam2bed
