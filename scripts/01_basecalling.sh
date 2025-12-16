#!/bin/bash
#SBATCH --job-name=nanopore_basecall
#SBATCH --partition=machinelearning
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:nvidia_h100_nvl:1
#SBATCH --time=48:00:00
#SBATCH --output=logs/basecall_%j.out
#SBATCH --error=logs/basecall_%j.err

# =============================================================================
# Nanopore Basecalling with m6A Modification Detection
# =============================================================================
# This script performs basecalling on POD5 files using Dorado with m6A
# modification detection, followed by barcode demultiplexing.
#
# Usage:
#   sbatch 01_basecalling.sh -i <pod5_dir> -o <output_dir> [-b <barcode_kit>]
#
# Requirements:
#   - Dorado (tested with v1.1.1)
#   - NVIDIA GPU (H100 recommended)
#   - samtools
# =============================================================================

set -e

# -----------------------------------------------------------------------------
# Default parameters
# -----------------------------------------------------------------------------
POD5_DIR=""
OUTPUT_DIR=""
BARCODE_KIT="SQK-RBK114-24"
DORADO_PATH=""  # Set your dorado path
MODEL_DIR=""    # Set your model directory
MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
MOD_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.2.0_6mA@v1"
SKIP_DEMUX=false

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
usage() {
    cat << EOF
Usage: sbatch $(basename $0) [OPTIONS]

Options:
    -i, --input DIR       Input directory containing POD5 files (required)
    -o, --output DIR      Output directory (required)
    -b, --barcode KIT     Barcode kit name (default: SQK-RBK114-24)
    -d, --dorado PATH     Path to dorado executable
    -m, --model-dir DIR   Directory containing model files
    --skip-demux          Skip demultiplexing step
    -h, --help            Show this help message

Example:
    sbatch $(basename $0) -i raw_data/pod5 -o results/run1 -b SQK-RBK114-24
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input) POD5_DIR="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        -b|--barcode) BARCODE_KIT="$2"; shift 2 ;;
        -d|--dorado) DORADO_PATH="$2"; shift 2 ;;
        -m|--model-dir) MODEL_DIR="$2"; shift 2 ;;
        --skip-demux) SKIP_DEMUX=true; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------
if [ -z "$POD5_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Input and output directories are required"
    usage
fi

if [ ! -d "$POD5_DIR" ]; then
    echo "ERROR: Input directory not found: $POD5_DIR"
    exit 1
fi

POD5_COUNT=$(ls -1 "$POD5_DIR"/*.pod5 2>/dev/null | wc -l)
if [ "$POD5_COUNT" -eq 0 ]; then
    echo "ERROR: No POD5 files found in $POD5_DIR"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# -----------------------------------------------------------------------------
# Print configuration
# -----------------------------------------------------------------------------
echo "============================================"
echo "Nanopore Basecalling Pipeline"
echo "============================================"
echo "Started at: $(date)"
echo ""
echo "Configuration:"
echo "  Input directory: $POD5_DIR"
echo "  POD5 files: $POD5_COUNT"
echo "  Output directory: $OUTPUT_DIR"
echo "  Barcode kit: $BARCODE_KIT"
echo "  Model: $MODEL"
echo "  m6A model: $MOD_MODEL"
echo "  Skip demux: $SKIP_DEMUX"
echo "============================================"
echo ""

# Check GPU
echo "GPU information:"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo ""

# -----------------------------------------------------------------------------
# Step 1: Basecalling with m6A detection
# -----------------------------------------------------------------------------
echo "[$(date)] Step 1: Basecalling with m6A detection..."

MODEL_PATH="$MODEL_DIR/$MODEL"
MOD_MODEL_PATH="$MODEL_DIR/$MOD_MODEL"
BASECALL_OUTPUT="$OUTPUT_DIR/basecalled.bam"

if [ "$SKIP_DEMUX" = true ]; then
    $DORADO_PATH basecaller \
        "$MODEL_PATH" \
        "$POD5_DIR" \
        --modified-bases-models "$MOD_MODEL_PATH" \
        > "$BASECALL_OUTPUT"
else
    $DORADO_PATH basecaller \
        "$MODEL_PATH" \
        "$POD5_DIR" \
        --modified-bases-models "$MOD_MODEL_PATH" \
        --kit-name "$BARCODE_KIT" \
        > "$BASECALL_OUTPUT"
fi

echo "  Generated: $BASECALL_OUTPUT ($(du -h "$BASECALL_OUTPUT" | cut -f1))"

# -----------------------------------------------------------------------------
# Step 2: Demultiplexing
# -----------------------------------------------------------------------------
if [ "$SKIP_DEMUX" = false ]; then
    echo ""
    echo "[$(date)] Step 2: Demultiplexing..."

    DEMUX_DIR="$OUTPUT_DIR/demultiplexed"
    mkdir -p "$DEMUX_DIR"

    $DORADO_PATH demux \
        --output-dir "$DEMUX_DIR" \
        --no-classify \
        "$BASECALL_OUTPUT"

    echo ""
    echo "Demultiplexing summary:"
    for bam in "$DEMUX_DIR"/*.bam; do
        if [ -f "$bam" ]; then
            bc=$(basename "$bam" .bam)
            count=$(samtools view -c "$bam")
            echo "  $bc: $count reads"
        fi
    done
fi

# -----------------------------------------------------------------------------
# Complete
# -----------------------------------------------------------------------------
echo ""
echo "============================================"
echo "Basecalling completed successfully!"
echo "Finished at: $(date)"
echo "============================================"
