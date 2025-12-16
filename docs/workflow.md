# Detailed Workflow Documentation

## Overview

This document provides detailed information about each step of the m6A detection pipeline.

## Step 1: Basecalling with Modification Detection

### Description

Dorado basecaller converts raw electrical signals (POD5) to nucleotide sequences while simultaneously detecting m6A modifications. The SUP (super accuracy) model provides the highest accuracy for both basecalling and modification detection.

### Command

```bash
dorado basecaller \
    dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
    <pod5_directory> \
    --modified-bases-models dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1 \
    --kit-name SQK-RBK114-24 \
    > output.bam
```

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Model | `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` | SUP accuracy basecalling model |
| Mod model | `dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1` | m6A modification detection model |
| Kit | `SQK-RBK114-24` | Rapid barcoding kit (24 barcodes) |

### Output

- BAM file with:
  - Basecalled sequences
  - Quality scores
  - MM/ML tags containing modification information
  - BC tag containing barcode classification

## Step 2: Demultiplexing

### Description

Separates reads by barcode into individual BAM files.

### Command

```bash
dorado demux \
    --output-dir demultiplexed/ \
    --no-classify \
    basecalled.bam
```

### Parameters

| Parameter | Description |
|-----------|-------------|
| `--no-classify` | Use existing barcode classifications from basecalling |
| `--output-dir` | Directory for demultiplexed BAM files |

### Output

- One BAM file per barcode: `barcode01.bam`, `barcode02.bam`, etc.
- `unclassified.bam` for reads without clear barcode assignment

## Step 3: Merging Technical Replicates

### Description

When samples are sequenced across multiple runs, this step combines BAM files from the same biological sample.

### Command

```bash
samtools merge -o merged_sample.bam run1/barcode01.bam run2/barcode10.bam
```

### Configuration

Create a TSV file mapping sample names to input BAM files:

```
sample_name    bam_file1    bam_file2
day9_rep1      run1/barcode01.bam    run2/barcode10.bam
day9_rep2      run1/barcode02.bam    run2/barcode11.bam
```

## Step 4: Alignment

### Description

Aligns reads to the reference genome while preserving modification tags.

### Commands

```bash
# Convert BAM to FASTQ, preserving MM/ML tags
samtools bam2fq -T MM,ML -@ 16 input.bam > reads.fastq

# Align with minimap2
# -ax map-ont: Oxford Nanopore preset
# -y: Copy input FASTQ comments to output SAM
minimap2 -ax map-ont -y -t 16 reference.fasta reads.fastq \
    | samtools sort -@ 16 -o aligned.bam -

# Index the BAM file
samtools index -@ 16 aligned.bam
```

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-ax map-ont` | - | ONT read alignment preset |
| `-y` | - | Preserve modification tags from FASTQ comments |
| `-t 16` | 16 threads | Parallel processing |

## Step 5: m6A Site Extraction

### Description

Extracts m6A modification sites from aligned BAM files using modbam2bed.

### Command

```bash
modbam2bed -e -m 6mA -t 16 -f 0.66 reference.fasta aligned.bam > m6a_sites.bed
```

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `-e` | - | Extended BED output with coverage and mod % |
| `-m 6mA` | - | Target modification type |
| `-t 16` | 16 | Number of threads |
| `-f 0.66` | 0.66 | Minimum modification probability threshold |

### Output Format

Extended BED format with columns:

1. `chrom` - Chromosome name
2. `start` - Start position (0-based)
3. `end` - End position
4. `name` - Modification name (a = 6mA)
5. `score` - Score (unused)
6. `strand` - + or -
7. `thickStart` - Same as start
8. `thickEnd` - Same as end
9. `color` - RGB color
10. `coverage` - Number of reads at this position
11. `mod_percent` - Percentage of reads with modification

### High-Confidence Filtering

```bash
# Filter sites with coverage ≥ 5 and modification rate ≥ 20%
awk '$10 >= 5 && $11 >= 20' m6a_sites.bed > m6a_high_conf.bed
```

## Quality Control Metrics

### Alignment Statistics

```bash
samtools flagstat aligned.bam
```

Key metrics to check:
- Total reads
- Mapped reads percentage
- Properly paired reads (if applicable)

### m6A Detection Summary

| Metric | Description |
|--------|-------------|
| Total sites | All positions with potential m6A |
| Sites with coverage | Positions with ≥1 read |
| Sites with m6A detected | Positions where m6A was called |
| High-confidence sites | Filtered sites meeting thresholds |
| Average modification rate | Mean % across detected sites |

## Computational Resources

### Typical Resource Usage

| Step | Time | Memory | Notes |
|------|------|--------|-------|
| Basecalling | ~2-4 hours/100k reads | 16GB GPU | GPU required |
| Demultiplexing | ~10 min | 8GB | CPU only |
| Merging | ~5 min/sample | 4GB | CPU only |
| Alignment | ~30 min/sample | 32GB | CPU intensive |
| m6A extraction | ~20 min/sample | 16GB | CPU only |

### SLURM Job Templates

See `scripts/` directory for ready-to-use SLURM submission scripts.

## Troubleshooting

### Common Issues

1. **No m6A sites detected**
   - Check that modification tags (MM/ML) are preserved through pipeline
   - Verify `-y` flag is used with minimap2
   - Confirm `-T MM,ML` is used with samtools bam2fq

2. **Low alignment rate**
   - Verify correct reference genome
   - Check read quality with `samtools stats`

3. **Out of memory**
   - Reduce thread count
   - Process samples sequentially
