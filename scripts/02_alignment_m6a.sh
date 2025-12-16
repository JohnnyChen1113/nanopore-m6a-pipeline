#!/bin/bash
#SBATCH --job-name=m6a_analysis
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/m6a_%j.out
#SBATCH --error=logs/m6a_%j.err

# =============================================================================
# Alignment and m6A Site Extraction
# =============================================================================
# This script aligns basecalled reads to a reference genome and extracts
# m6A modification sites using modbam2bed.
#
# Usage:
#   sbatch 02_alignment_m6a.sh -i <input.bam> -o <output_prefix> -r <reference.fasta>
#
# Requirements:
#   - samtools
#   - minimap2
#   - modbam2bed
# =============================================================================

set -e

# -----------------------------------------------------------------------------
# Default parameters
# -----------------------------------------------------------------------------
INPUT_BAM=""
OUTPUT_PREFIX=""
REFERENCE=""
THREADS=16
M6A_THRESHOLD=0.66
MIN_COVERAGE=5
MIN_MOD_RATE=20

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: sbatch $(basename $0) [OPTIONS]

Options:
    -i, --input BAM         Input BAM file with modification tags (required)
    -o, --output PREFIX     Output file prefix (required)
    -r, --reference FASTA   Reference genome FASTA file (required)
    -t, --threads NUM       Number of threads (default: 16)
    -f, --threshold NUM     m6A detection threshold (default: 0.66)
    --min-cov NUM           Minimum coverage for high-confidence sites (default: 5)
    --min-mod NUM           Minimum modification rate % for high-confidence (default: 20)
    -h, --help              Show this help message

Example:
    sbatch $(basename $0) -i demux/barcode01.bam -o results/sample1 -r reference.fasta
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_BAM="$2"; shift 2 ;;
        -o|--output) OUTPUT_PREFIX="$2"; shift 2 ;;
        -r|--reference) REFERENCE="$2"; shift 2 ;;
        -t|--threads) THREADS="$2"; shift 2 ;;
        -f|--threshold) M6A_THRESHOLD="$2"; shift 2 ;;
        --min-cov) MIN_COVERAGE="$2"; shift 2 ;;
        --min-mod) MIN_MOD_RATE="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------
if [ -z "$INPUT_BAM" ] || [ -z "$OUTPUT_PREFIX" ] || [ -z "$REFERENCE" ]; then
    echo "ERROR: Input BAM, output prefix, and reference are required"
    usage
fi

if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: Input BAM not found: $INPUT_BAM"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE"
    exit 1
fi

# Create output directory
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# -----------------------------------------------------------------------------
# Print configuration
# -----------------------------------------------------------------------------
echo "============================================"
echo "Alignment and m6A Extraction Pipeline"
echo "============================================"
echo "Started at: $(date)"
echo ""
echo "Configuration:"
echo "  Input BAM: $INPUT_BAM"
echo "  Output prefix: $OUTPUT_PREFIX"
echo "  Reference: $REFERENCE"
echo "  Threads: $THREADS"
echo "  m6A threshold: $M6A_THRESHOLD"
echo "  Min coverage: $MIN_COVERAGE"
echo "  Min mod rate: $MIN_MOD_RATE%"
echo "============================================"
echo ""

# Define output files
TEMP_FASTQ="${OUTPUT_PREFIX}_temp.fastq"
ALIGNED_BAM="${OUTPUT_PREFIX}_aligned.bam"
M6A_BED="${OUTPUT_PREFIX}_m6a_sites.bed"
HIGH_CONF_BED="${OUTPUT_PREFIX}_m6a_high_conf.bed"
STATS_FILE="${OUTPUT_PREFIX}_alignment_stats.txt"

# -----------------------------------------------------------------------------
# Step 1: Convert BAM to FASTQ (preserving modification tags)
# -----------------------------------------------------------------------------
echo "[$(date)] Step 1/5: Converting BAM to FASTQ..."
samtools bam2fq -T MM,ML -@ "$THREADS" "$INPUT_BAM" > "$TEMP_FASTQ"
echo "  Generated: $TEMP_FASTQ"

# -----------------------------------------------------------------------------
# Step 2: Align reads with minimap2
# -----------------------------------------------------------------------------
echo "[$(date)] Step 2/5: Aligning reads with minimap2..."
minimap2 -ax map-ont -y -t "$THREADS" "$REFERENCE" "$TEMP_FASTQ" \
    | samtools sort -@ "$THREADS" -o "$ALIGNED_BAM" -
echo "  Generated: $ALIGNED_BAM ($(du -h "$ALIGNED_BAM" | cut -f1))"

# Clean up temporary FASTQ
rm -f "$TEMP_FASTQ"

# -----------------------------------------------------------------------------
# Step 3: Index BAM file
# -----------------------------------------------------------------------------
echo "[$(date)] Step 3/5: Indexing BAM..."
samtools index -@ "$THREADS" "$ALIGNED_BAM"
echo "  Generated: ${ALIGNED_BAM}.bai"

# -----------------------------------------------------------------------------
# Step 4: Generate alignment statistics
# -----------------------------------------------------------------------------
echo "[$(date)] Step 4/5: Generating alignment statistics..."
samtools flagstat "$ALIGNED_BAM" > "$STATS_FILE"
echo "  Generated: $STATS_FILE"
echo ""
echo "Alignment summary:"
cat "$STATS_FILE"
echo ""

# -----------------------------------------------------------------------------
# Step 5: Extract m6A modification sites
# -----------------------------------------------------------------------------
echo "[$(date)] Step 5/5: Extracting m6A sites with modbam2bed..."
modbam2bed -e -m 6mA -t "$THREADS" -f "$M6A_THRESHOLD" \
    "$REFERENCE" \
    "$ALIGNED_BAM" \
    > "$M6A_BED"

TOTAL_SITES=$(wc -l < "$M6A_BED")
echo "  Generated: $M6A_BED ($TOTAL_SITES sites)"

# Filter high-confidence sites
# BED columns: chrom, start, end, name, score, strand, start, end, color, coverage, mod_percent
awk -v cov="$MIN_COVERAGE" -v mod="$MIN_MOD_RATE" \
    '$10 >= cov && $11 >= mod' "$M6A_BED" > "$HIGH_CONF_BED"

HIGH_CONF_SITES=$(wc -l < "$HIGH_CONF_BED")
echo "  Generated: $HIGH_CONF_BED ($HIGH_CONF_SITES high-confidence sites)"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "============================================"
echo "Pipeline completed successfully!"
echo "Finished at: $(date)"
echo "============================================"
echo ""
echo "Output files:"
echo "  - Aligned BAM: $ALIGNED_BAM"
echo "  - BAM index: ${ALIGNED_BAM}.bai"
echo "  - Alignment stats: $STATS_FILE"
echo "  - All m6A sites: $M6A_BED ($TOTAL_SITES sites)"
echo "  - High-confidence m6A: $HIGH_CONF_BED ($HIGH_CONF_SITES sites)"
echo "============================================"
