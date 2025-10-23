#!/bin/bash
###############################################################################
# Script: scrnaseq_data_download_sra_en.sh
# Author: Raúl [Bioinformatician]
# Description:
#   Automated pipeline to download and process scRNA-seq data from the
#   NCBI Sequence Read Archive (SRA) using prefetch and fasterq-dump.
#   Includes parameter validation, timestamped logging, parallel compression,
#   and optional read statistics using seqkit.
#   Update 23-10-2025: Now this function detects the average length of each
#   file and rename each file with the CellRanger convenction based on sequencing length.
#
# Important:
#   ⚠️ Run this script *from the directory where you want all files to be saved*.
#   The SRA Toolkit commands (prefetch, fasterq-dump) create subdirectories
#   relative to the current working directory.
#
# Usage:
#   ./sra_downloader.sh -i accession_list.txt -o output_dir [-t threads]
#
# Dependencies:
#   - SRA Toolkit (prefetch, fasterq-dump)
#   - pigz (for parallel compression)
#   - seqkit (optional, for read statistics)
#
# Parameters:
#   -i  File containing SRA accessions (one per line)
#   -o  Output directory for FASTQ files
#   -t  Number of threads (default: 8, max: 32)
#   -h  Show this help message
#
# Example:
#   ./sra_downloader.sh -i samples.txt -o data/fastq -t 16
#
# Change log:
#   v1.0 - Initial release with error handling, logging, and compression
#   v2.0 - Rename automatically all files with Cell Ranger convention
###############################################################################

# =========================
# Logging function with timestamp
# =========================
log() {
    local type=$1
    local message=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$type] $message"
}

# =========================
# Format time (seconds → Xm Ys)
# =========================
format_time() {
    local seconds=$1
    local minutes=$((seconds / 60))
    local remainder=$((seconds % 60))
    echo "${minutes}m ${remainder}s"
}

# =========================
# Display help message
# =========================
show_help(){
    echo "Usage: $0 -i accession_list.txt -o output_dir [-t threads]"
    echo ""
    echo "Options:"
    echo " -i ACCESSION_LIST   File with SRR accessions"
    echo " -o OUTPUT_DIR       Output directory"
    echo " -t THREADS          Number of threads (default: 8, max: 32)"
    echo " -h                  Show this help"
    echo ""
}

# =========================
# Main download function
# =========================
download_sra_files(){
    local accession_file=$1
    local output_dir=$2
    local threads=$3

    mkdir -p "$output_dir"
    mkdir -p "$output_dir/tmp"

    log "INFO" "Starting download of samples from: $accession_file"

    # Read each line (sample accession)
    while IFS= read -r line; do
        sample=$(echo "$line" | tr -d '\r\n')
        [[ -z "$sample" ]] && { log "WARN" "Empty line detected, skipping..."; continue; }

        log "INFO" "=============================================================================="
        log "INFO" "Processing sample: $sample"

        # --- PREFETCH ---
        log "INFO" "Running prefetch for $sample"
        local start_time=$(date +%s)
        prefetch "$sample" --progress
        local exit_code=$?
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        local time_spent=$(format_time "$duration")

        if [[ $exit_code -ne 0 ]]; then
            log "ERROR" "prefetch failed for $sample (exit code: $exit_code) after $time_spent"
            continue
        else
            log "SUCCESS" "prefetch completed in $time_spent"
        fi

        # --- FASTERQ-DUMP ---
        log "INFO" "Running fasterq-dump for $sample"
        start_time=$(date +%s)
        fasterq-dump "$sample" --split-files --include-technical -O "$output_dir" --threads "$threads" --temp "$output_dir/tmp"
        exit_code=$?
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        time_spent=$(format_time "$duration")

        if [[ $exit_code -ne 0 ]]; then
            log "ERROR" "fasterq-dump failed for $sample (exit code: $exit_code) after $time_spent"
            continue
        else
            log "SUCCESS" "fasterq-dump completed in $time_spent"
        fi

        # --- CLEANUP ---
        if [[ -d "$sample" ]]; then
            rm -rf "$sample"
            log "INFO" "Temporary folder '$sample' removed."
        fi

        # --- COMPRESSION ---
        log "INFO" "Compressing FASTQ files for $sample"
        start_time=$(date +%s)
        if compress_fastq "$sample" "$output_dir"; then
            end_time=$(date +%s)
            duration=$((end_time - start_time))
            time_spent=$(format_time "$duration")
            log "SUCCESS" "Compression completed in $time_spent"
        else
            log "ERROR" "Compression failed for $sample"
        fi

    done < "$accession_file"

    # --- REMOVE TEMP DIR ---
    if [[ -d "$output_dir/tmp" ]]; then
        rm -rf "$output_dir/tmp"
        log "INFO" "Temporary directory '$output_dir/tmp' removed."
    fi
}

# =========================
# Compress FASTQ files using pigz
# =========================
compress_fastq() {
    local sample=$1
    local output_dir=$2
    
    for fastq_file in "${output_dir}/${sample}"_*.fastq; do
        if [[ -f "$fastq_file" ]]; then
            log "INFO" "Compressing: $(basename "$fastq_file")"
            pigz "$fastq_file"
            if [[ $? -ne 0 ]]; then
                log "ERROR" "Error compressing $fastq_file"
                return 1
            fi
        fi
    done
    return 0
}

# --------------------------------------------------
# GLOBAL RENAME FUNCTION (AFTER ALL SAMPLES)
# --------------------------------------------------
renombrar_todos_cellranger() {
    local stats_file=$1
    local output_dir=$2

    log "INFO" "Starting global FASTQ renaming to Cell Ranger convention"

    if [[ ! -f "$stats_file" ]]; then
        log "ERROR" "Statistics file not found: $stats_file"
        return 1
    fi

    # Normalize line endings (remueve \r si vienen de Windows)
    # y procesar con awk; awk maneja bien números decimales.
    # Campos esperados por seqkit stats: file format type num_seqs sum_len min_len avg_len max_len

    # Usamos awk para:
    #  - saltar la cabecera (NR==1)
    #  - extraer basefile (basename)
    #  - decidir tipo (primero por patrón en el nombre, si no, por avg_len)
    #  - generar y ejecutar mv solo si el archivo existe
    awk -v outdir="$output_dir" '
    BEGIN { OFS = ""; }
    NR==1 { next } # saltar header
    {
        # eliminar \r que puede venir en sistemas con CRLF
        sub(/\r$/, "", $1)

        # file puede venir con path absoluto/relativo; extraer basename
        file = $1
        n = split(file, parts, "/")
        base = parts[n]

        # avg_len está en la 8ª columna (si seqkit stats clásico)
        avg = $8 + 0

        # decidir tipo por patrón en el nombre
        tipo = "UNK"
        if (base ~ /(_1\.fastq|_R1_|\_R1\.|_R1\.fastq)/) tipo = "R1"
        else if (base ~ /(_2\.fastq|_R2_|\_R2\.|_R2\.fastq)/) tipo = "R2"
        else if (base ~ /\.I1\.|_I1\./) tipo = "I1"  # por si hay indicadores técnicos
        else {
            # fallback a avg_len si no se puede por patrón
            if (avg > 60) tipo = "R2"
            else if (avg > 20 && avg < 40) tipo = "R1"
            else if (avg <= 20) tipo = "I1"
            else tipo = "UNK"
        }

        if (tipo == "UNK") {
            printf("CMD: echo \"WARN Unknown read type for %s (avg_len=%s); skipping\\n\"", base, avg) | "cat"
            next
        }

        # Construir sample: quitar sufijo de read y extensión
        # ejemplo: sample_S1_L001_R1_001.fastq.gz -> sample_S1_L001
        sample = base
        # quitar .fastq(.gz)
        sub(/\.fastq(\.gz)?$/, "", sample)
        # quitar read suffixes típicos
        # posibles sufijos: _R1_001, _R2_001, _1, _2, _S1_L001_R1_001, etc.
        sub(/(_R?[12](_[0-9]+)?$)|(_[12]$)|(_R[12]_[0-9]+$)/, "", sample)

        nuevo = outdir "/" sample "_S1_L001_" tipo "_001.fastq.gz"
        origen = outdir "/" base

        # imprimir y ejecutar solo si el archivo existe
        cmd = "test -f \"" origen "\""
        if (system(cmd) == 0) {
            # mv con -v para logging; usar sh -c para evitar problemas con comillas en awk print
            printf("mv -v \"%s\" \"%s\"\n", origen, nuevo) | "sh"
            printf("CMD: echo \"SUCCESS Renamed: %s -> %s\\n\"", base, sample "_S1_L001_" tipo "_001.fastq.gz") | "cat"
        } else {
            printf("CMD: echo \"WARN Source not found: %s (skipping)\\n\"", origen) | "cat"
        }
    }
    ' <(tr -d '\r' < "$stats_file")
    
    log "INFO" "Global renaming completed"
}

# =========================
# Default variables and argument parsing
# =========================
ACCESSION_FILE=""
OUTPUT_DIR=""
THREADS=8  # default value

while getopts "i:o:t:h" opt; do
    case $opt in 
        i) ACCESSION_FILE="$OPTARG";;
        o) OUTPUT_DIR="$OPTARG";;
        t) THREADS="$OPTARG";;
        h) show_help; exit 0;;
        *) log "ERROR" "Invalid option"; show_help; exit 1;;
    esac
done

# =========================
# Parameter validation
# =========================
if [[ -z "$ACCESSION_FILE" ]]; then
    log "ERROR" "You must specify the accession list file with -i"
    show_help; exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    log "ERROR" "You must specify an output directory with -o"
    show_help; exit 1
fi

if [[ ! -f "$ACCESSION_FILE" ]]; then
    log "ERROR" "The file '$ACCESSION_FILE' does not exist"
    exit 1
fi

if [[ ! "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]] || [[ "$THREADS" -gt 32 ]]; then
    log "ERROR" "Threads must be an integer between 1 and 32"
    exit 1
fi

# =========================
# Main execution
# =========================
log "INFO" "Starting SRA download pipeline"
download_sra_files "$ACCESSION_FILE" "$OUTPUT_DIR" "$THREADS"

# --- Generate summary report using seqkit ---
log "INFO" "Generating read statistics"
start_time=$(date +%s)
seqkit stats "$OUTPUT_DIR"/*.fastq.gz > "$OUTPUT_DIR"/reads_stats.txt 2>/dev/null
end_time=$(date +%s)
duration=$((end_time - start_time))
time_spent=$(format_time "$duration")

if [[ $? -eq 0 ]]; then
    log "SUCCESS" "Report generated in ${time_spent}: $OUTPUT_DIR/reads_stats.txt"
else
    log "WARN" "Statistics report failed (is seqkit installed?)"
fi

log "INFO" "Pipeline completed successfully"

