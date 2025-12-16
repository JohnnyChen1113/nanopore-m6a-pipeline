#!/bin/bash
# =============================================================================
# Merge BAM Files from Multiple Sequencing Runs
# =============================================================================
# This script merges BAM files from multiple sequencing runs for the same
# biological sample, enabling technical replicate combination.
#
# Usage:
#   bash 03_merge_samples.sh -c <config.tsv> -o <output_dir>
#
# Requirements:
#   - samtools
#
# Config file format (TSV):
#   sample_name<TAB>bam_file1<TAB>bam_file2<TAB>...
# =============================================================================

set -e

# -----------------------------------------------------------------------------
# Default parameters
# -----------------------------------------------------------------------------
CONFIG_FILE=""
OUTPUT_DIR=""

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: bash $(basename $0) [OPTIONS]

Options:
    -c, --config FILE     Config file with sample-to-BAM mapping (required)
    -o, --output DIR      Output directory for merged BAM files (required)
    -h, --help            Show this help message

Config file format (TSV, no header):
    sample_name<TAB>bam_path1<TAB>bam_path2<TAB>...

Example config (merge_config.tsv):
    day9_rep1    run1/barcode01.bam    run2/barcode10.bam
    day9_rep2    run1/barcode02.bam    run2/barcode11.bam
    day21_rep1   run1/barcode04.bam    run2/barcode13.bam

Example:
    bash $(basename $0) -c merge_config.tsv -o merged_samples/
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config) CONFIG_FILE="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------
if [ -z "$CONFIG_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Config file and output directory are required"
    usage
fi

if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# -----------------------------------------------------------------------------
# Process each sample
# -----------------------------------------------------------------------------
echo "============================================"
echo "Merge BAM Files Pipeline"
echo "============================================"
echo "Started at: $(date)"
echo "Config file: $CONFIG_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "============================================"
echo ""

SAMPLE_COUNT=0
while IFS=$'\t' read -r sample_name bam_files; do
    # Skip empty lines and comments
    [[ -z "$sample_name" || "$sample_name" =~ ^# ]] && continue

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    echo "----------------------------------------"
    echo "Processing sample: $sample_name"
    echo "----------------------------------------"

    # Convert tab-separated BAM files to array
    IFS=$'\t' read -ra BAM_ARRAY <<< "$bam_files"

    # Check if all BAM files exist
    VALID=true
    for bam in "${BAM_ARRAY[@]}"; do
        if [ ! -f "$bam" ]; then
            echo "  WARNING: BAM file not found: $bam"
            VALID=false
        fi
    done

    if [ "$VALID" = false ]; then
        echo "  SKIPPED: Missing input files"
        continue
    fi

    OUTPUT_BAM="$OUTPUT_DIR/${sample_name}.bam"

    # Merge BAM files
    if [ ${#BAM_ARRAY[@]} -eq 1 ]; then
        # Single file, just copy
        echo "  Single input file, copying..."
        cp "${BAM_ARRAY[0]}" "$OUTPUT_BAM"
    else
        # Multiple files, merge
        echo "  Merging ${#BAM_ARRAY[@]} BAM files..."
        samtools merge -f "$OUTPUT_BAM" "${BAM_ARRAY[@]}"
    fi

    # Count reads
    READ_COUNT=$(samtools view -c "$OUTPUT_BAM")
    echo "  Output: $OUTPUT_BAM ($READ_COUNT reads)"
    echo ""

done < "$CONFIG_FILE"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo "============================================"
echo "Merge completed!"
echo "Processed $SAMPLE_COUNT samples"
echo "Finished at: $(date)"
echo "============================================"
