#!/bin/bash
#
# ==============================================================
# Script: make_aggr_csv.sh
# Author: Raúl Fernández
# Version: 1.0
#
# Description:
#   Automatically generates the CSV file required by Cell Ranger
#   for the `cellranger aggr` command, based on the existing
#   `molecule_info.h5` files in multiple sample directories.
#
# Features:
#     - Recursively searches for .h5 files in a given base directory.
#     - Extracts the sample ID (the folder containing "outs").
#     - Optionally replaces technical IDs with real sample names
#       using a user-provided mapping file.
#     - Allows choosing the output file path and name.
#
# Usage:
#   ./make_aggr_csv.sh -i <base_directory> -o <output_csv> [-m <mapping_file>]
#
# Example:
#   ./make_aggr_csv.sh \
#     -i /your/path/here \
#     -o /your/path/here/aggregated_samples.csv \
#     -m /your/path/here/map_samples.csv
#
# Optional mapping file format (CSV with header):
#   sample_id_real,sample_id_cellranger
#   Sample1,SAMN123456777
#   Sample2,SAMN123456788
#   Sample3,SAMN123456789
#
# Expected output:
#   sample_id,molecule_h5
#   Sample1,/path/SAMN123456777/outs/molecule_info.h5
#   Sample2,/path/SAMN123456788/outs/molecule_info.h5
#   Sample3,/path/SAMN123456789/outs/molecule_info.h5
#
# Requirements:
#   bash >= 4 (for associative arrays)
# ==============================================================

# -------------------------------
# 1️⃣ Parse command-line arguments
# -------------------------------
while getopts "i:o:m:" opt; do
  case $opt in
    i) BASE_DIR="$OPTARG" ;;   # Base directory containing the .h5 files
    o) OUTPUT_FILE="$OPTARG" ;; # Output CSV file
    m) MAP_FILE="$OPTARG" ;;    # Optional mapping file
    *)
      echo "Usage: $0 -i <input_dir> -o <output_csv> [-m <map_file>]"
      exit 1
      ;;
  esac
done

# Validate required parameters
if [ -z "$BASE_DIR" ] || [ -z "$OUTPUT_FILE" ]; then
  echo "❌ Error: You must specify at least -i (input dir) and -o (output csv)."
  echo "Usage: $0 -i <input_dir> -o <output_csv> [-m <map_file>]"
  exit 1
fi

# Check that the base directory exists
if [ ! -d "$BASE_DIR" ]; then
  echo "❌ Error: Base directory '$BASE_DIR' does not exist."
  exit 1
fi

# Create output directory if it doesn’t exist
OUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUT_DIR"

# -------------------------------
# 2️⃣ Load sample ID mapping (optional)
# -------------------------------
declare -A MAP  # Associative array to store mappings

if [ -n "$MAP_FILE" ] && [ -f "$MAP_FILE" ]; then
  echo "Reading mapping file: $MAP_FILE"
  # Skip the header line and read key-value pairs
  while IFS=, read -r real_id tech_id; do
    if [[ -n "$real_id" && -n "$tech_id" ]]; then
      MAP["$tech_id"]="$real_id"
    fi
  done < <(tail -n +2 "$MAP_FILE")
else
  if [ -n "$MAP_FILE" ]; then
    echo "⚠️  Warning: Mapping file '$MAP_FILE' not found. Original IDs will be used."
  fi
fi

# -------------------------------
# 3️⃣ Generate the CSV file
# -------------------------------
echo "sample_id,molecule_h5" > "$OUTPUT_FILE"

# Find all molecule_info.h5 files recursively
find "$BASE_DIR" -type f -name "molecule_info.h5" | while read -r FILE; do
  # Folder that directly contains 'outs'
  SAMPLE_ID=$(basename "$(dirname "$FILE")")

  # Parent folder of that folder (e.g., the true sample name)
  PARENT_DIR=$(basename "$(dirname "$(dirname "$FILE")")")

  # Determine final sample ID:
  # If mapping exists, use the mapped name; otherwise use the parent directory name.
  if [ -n "${MAP[$PARENT_DIR]}" ]; then
    FINAL_ID="${MAP[$PARENT_DIR]}"
  else
    FINAL_ID="$PARENT_DIR"
  fi

  # Append line to CSV
  echo "$FINAL_ID,$FILE" >> "$OUTPUT_FILE"
done

# -------------------------------
# 4️⃣ Final message
# -------------------------------
echo "✅ File '$OUTPUT_FILE' successfully generated with $(($(wc -l < "$OUTPUT_FILE") - 1)) samples."

