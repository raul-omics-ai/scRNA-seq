#!/bin/bash

# ===============================
# Script: run_cellranger_count_v3.sh
# Author: Raúl
# Description: Runs Cell Ranger 'count' on selected FASTQ samples,
#              moves only the 'outs' directory to the output folder,
#              and supports automatic resume if interrupted.
# ===============================

# -------------------------------
# Help / usage function
# -------------------------------
usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -t <threads> -r <reference_path> -s <samples.txt> [-m <memory_gb>]"
    echo
    echo "Options:"
    echo "  -i   Directory containing compressed FASTQ files (.fastq.gz)"
    echo "  -o   Destination directory where the 'outs' folders will be moved"
    echo "  -t   Number of CPU threads (cores) to use for Cell Ranger"
    echo "  -r   Path to the reference transcriptome (e.g., /path/to/refdata-gex-GRCh38-2020-A)"
    echo "  -s   Text file with sample names to process (one per line)"
    echo "  -m   Local memory in GB (default: 64)"
    echo
    echo "Example:"
    echo "  $0 -i ./fastqs -o /mnt/results -t 8 -r /refs/refdata-gex-GRCh38-2020-A -s samples.txt -m 96"
    echo
    exit 1
}

# -------------------------------
# Default values
# -------------------------------
MEMORY=64

# -------------------------------
# Parse command-line arguments
# -------------------------------
while getopts ":i:o:t:r:m:s:" opt; do
  case ${opt} in
    i ) INPUT_DIR="$OPTARG" ;;
    o ) OUTPUT_DIR="$OPTARG" ;;
    t ) THREADS="$OPTARG" ;;
    r ) REFERENCE="$OPTARG" ;;
    m ) MEMORY="$OPTARG" ;;
    s ) SAMPLE_FILE="$OPTARG" ;;
    \? ) echo "Invalid option: -$OPTARG" >&2; usage ;;
    : ) echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# -------------------------------
# Check required arguments
# -------------------------------
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$THREADS" || -z "$REFERENCE" || -z "$SAMPLE_FILE" ]]; then
    usage
fi

# -------------------------------
# Environment setup
# -------------------------------
echo "📂 Input directory: $INPUT_DIR"
echo "📦 Output directory: $OUTPUT_DIR"
echo "🧬 Reference transcriptome: $REFERENCE"
echo "📋 Sample list file: $SAMPLE_FILE"
echo "🧵 Threads: $THREADS"
echo "💾 Local memory: ${MEMORY}GB"
echo

mkdir -p "$OUTPUT_DIR"

# -------------------------------
# Read sample list
# -------------------------------
if [[ ! -f "$SAMPLE_FILE" ]]; then
    echo "❌ Sample list file not found: $SAMPLE_FILE"
    exit 1
fi

mapfile -t SAMPLES < "$SAMPLE_FILE"
TOTAL=${#SAMPLES[@]}

if [[ $TOTAL -eq 0 ]]; then
    echo "❌ No sample names found in $SAMPLE_FILE"
    exit 1
fi

echo "🔍 Found $TOTAL samples to process."
echo

# -------------------------------
# Main loop
# -------------------------------
count=0
for SAMPLE in "${SAMPLES[@]}"; do
    count=$((count+1))
    SAMPLE=$(echo "$SAMPLE" | xargs)  # trim whitespace
    [[ -z "$SAMPLE" ]] && continue

    echo "==========================================="
    echo "➡️  Processing sample: $SAMPLE ($count / $TOTAL)"
    echo "==========================================="

    # If outs directory already exists in output, skip
    if [[ -d "$OUTPUT_DIR/$SAMPLE/outs" ]]; then
        echo "⚠️  Sample $SAMPLE already processed — skipping."
        continue
    fi

    # Run Cell Ranger
    cellranger count \
    	--no-bam \
        --id="$SAMPLE" \
        --transcriptome="$REFERENCE" \
        --fastqs="$INPUT_DIR" \
        --sample="$SAMPLE" \
        --localcores="$THREADS" \
        --localmem="$MEMORY" \
        &> "${SAMPLE}.log"

    # Check exit status
    if [[ $? -ne 0 ]]; then
        echo "❌ Error processing sample $SAMPLE. See ${SAMPLE}.log for details."
        continue
    fi

    # Ensure output directories exist
    mkdir -p "$OUTPUT_DIR/$SAMPLE"

    # Move only the 'outs' folder to the final output directory
    if [[ -d "$SAMPLE/outs" ]]; then
        mv "$SAMPLE/outs" "$OUTPUT_DIR/$SAMPLE/"
        echo "📁 Moved 'outs' folder for sample $SAMPLE → $OUTPUT_DIR/$SAMPLE/"
    else
        echo "⚠️  No 'outs' folder found for sample $SAMPLE. Skipping move."
        continue
    fi

    # Optionally clean up temporary files
    rm -rf "$SAMPLE"

    # Progress bar
    PERCENT=$(( count * 100 / TOTAL ))
    echo -ne "Progress: ["
    for ((i=0; i<=$PERCENT/5; i++)); do echo -ne "#"; done
    for ((i=$PERCENT/5; i<20; i++)); do echo -ne " "; done
    echo -ne "] $PERCENT% \r"
    sleep 0.3
done

echo
echo "✅ All selected samples processed successfully."
echo "Results (outs folders) are available in: $OUTPUT_DIR"

