#!/bin/bash
# ================================================================
# Title: scRNA-seq SRA Downloader
# Description: Automates downloading, decompressing, compressing, 
#              and renaming single-cell RNA-seq FASTQ files 
#              from SRA, preparing them for Cell Ranger.
#
# Author: Raúl
# Version: 1.2
# ================================================================

# --------------------------------------------------
# LOGGING FUNCTION (with timestamp)
# --------------------------------------------------
log() {
    local tipo=$1
    local mensaje=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$tipo] $mensaje"
}

# --------------------------------------------------
# TIME FORMATTER (seconds → Xm Ys)
# --------------------------------------------------
formatear_tiempo() {
    local segundos=$1
    local minutos=$((segundos / 60))
    local resto=$((segundos % 60))
    echo "${minutos}m ${resto}s"
}

# --------------------------------------------------
# HELP MESSAGE
# --------------------------------------------------
mostrar_ayuda(){
    echo "Usage: $0 -i accession_list.txt -o output_dir [-t threads]"
    echo ""
    echo "Options:"
    echo " -i ACCESSION_LIST   File containing SRR accessions"
    echo " -o OUTPUT_DIR       Output directory"
    echo " -t THREADS          Number of threads (default: 8, max: 32)"
    echo " -h                  Show this help message"
    echo ""
    echo "NOTE: Move and run this script from the directory where you"
    echo "      want to store your FASTQ files."
    echo ""
}

# --------------------------------------------------
# MAIN DOWNLOAD FUNCTION
# --------------------------------------------------
download_sra_files(){
    local archivo=$1
    local output_dir=$2
    local threads=$3

    mkdir -p "$output_dir"
    mkdir -p "$output_dir/tmp"

    log "INFO" "Starting download from accession list: $archivo"

    while IFS= read -r linea; do
        sample=$(echo "$linea" | tr -d '\r\n')
        [[ -z "$sample" ]] && { log "WARN" "Empty line detected, skipping..."; continue; }

        log "INFO" "=============================================================================="
        log "INFO" "Processing sample: $sample"

        # -------------------
        # PREFETCH
        # -------------------
        log "INFO" "Running prefetch for $sample"
        local start_time=$(date +%s)
        prefetch "$sample" --progress
        local exit_code=$?
        local end_time=$(date +%s)
        local duration=$((end_time - start_time))
        local tiempo=$(formatear_tiempo "$duration")

        if [[ $exit_code -ne 0 ]]; then
            log "ERROR" "prefetch failed for $sample (code: $exit_code) after $tiempo"
            continue
        else
            log "SUCCESS" "prefetch completed in $tiempo"
        fi

        # -------------------
        # FASTERQ-DUMP
        # -------------------
        log "INFO" "Running fasterq-dump for $sample"
        start_time=$(date +%s)
        fasterq-dump "$sample" --split-files --include-technical -O "$output_dir" --threads "$threads" --temp "$output_dir/tmp"
        exit_code=$?
        end_time=$(date +%s)
        duration=$((end_time - start_time))
        tiempo=$(formatear_tiempo "$duration")

        if [[ $exit_code -ne 0 ]]; then
            log "ERROR" "fasterq-dump failed for $sample (code: $exit_code) after $tiempo"
            continue
        else
            log "SUCCESS" "fasterq-dump completed in $tiempo"
        fi

        # -------------------
        # REMOVE .sra DIRECTORY
        # -------------------
        if [[ -d "$sample" ]]; then
            rm -rf "$sample"
            log "INFO" "Temporary folder '$sample' removed."
        fi

        # -------------------
        # COMPRESSION
        # -------------------
        log "INFO" "Compressing FASTQ files for $sample"
        start_time=$(date +%s)
        if comprimir_fastq "$sample" "$output_dir"; then
            end_time=$(date +%s)
            duration=$((end_time - start_time))
            tiempo=$(formatear_tiempo "$duration")
            log "SUCCESS" "Compression completed in $tiempo"
        else
            log "ERROR" "Error compressing $sample"
        fi

    done < "$archivo"

    # -------------------
    # REMOVE TEMPORARY DIRECTORY
    # -------------------
    if [[ -d "$output_dir/tmp" ]]; then
        rm -rf "$output_dir/tmp"
        log "INFO" "Temporary folder '$output_dir/tmp' removed."
    fi
}

# --------------------------------------------------
# FASTQ COMPRESSION FUNCTION
# --------------------------------------------------
comprimir_fastq() {
    local sample=$1
    local output_dir=$2
    
    for fastq_file in "${output_dir}/${sample}"_*.fastq; do
        if [[ -f "$fastq_file" ]]; then
            log "INFO" "Compressing: $(basename "$fastq_file")"
            pigz "$fastq_file"
            if [[ $? -ne 0 ]]; then
                log "ERROR" "Compression failed for $fastq_file"
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

# --------------------------------------------------
# ARGUMENT PARSING AND VALIDATION
# --------------------------------------------------
ARCHIVO=""
OUTPUTDIR=""
THREADS=8

while getopts "i:o:t:h" opcion; do
    case $opcion in 
        i) ARCHIVO="$OPTARG";;
        o) OUTPUTDIR="$OPTARG";;
        t) THREADS="$OPTARG";;
        h) mostrar_ayuda; exit 0;;
        *) log "ERROR" "Invalid option"; mostrar_ayuda; exit 1;;
    esac
done

if [[ -z "$ARCHIVO" ]]; then
    log "ERROR" "You must specify the accession list with -i"
    mostrar_ayuda
    exit 1
fi

if [[ -z "$OUTPUTDIR" ]]; then
    log "ERROR" "You must specify an output directory with -o"
    mostrar_ayuda
    exit 1
fi

if [[ ! -f "$ARCHIVO" ]]; then
    log "ERROR" "Accession list '$ARCHIVO' does not exist"
    exit 1
fi

if [[ ! "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]] || [[ "$THREADS" -gt 32 ]]; then
    log "ERROR" "Threads must be a number between 1 and 32"
    exit 1
fi

# --------------------------------------------------
# MAIN EXECUTION
# --------------------------------------------------
log "INFO" "Starting SRA scRNA-seq download pipeline"
download_sra_files "$ARCHIVO" "$OUTPUTDIR" "$THREADS"

# --------------------------------------------------
# GENERATE FINAL STATS REPORT
# --------------------------------------------------
log "INFO" "Generating summary statistics"
start_time=$(date +%s)
seqkit stats "$OUTPUTDIR"/*.fastq.gz > "$OUTPUTDIR"/reads_stats.txt
end_time=$(date +%s)
duration=$((end_time - start_time))
tiempo=$(formatear_tiempo "$duration")

if [[ $? -eq 0 ]]; then
    log "SUCCESS" "Statistics generated in ${tiempo}s: $OUTPUTDIR/reads_stats.txt"
else
    log "WARN" "Could not generate statistics (is seqkit installed?)"
fi

# --------------------------------------------------
# RENAME FILES TO CELL RANGER FORMAT
# --------------------------------------------------
renombrar_todos_cellranger "$OUTPUTDIR/reads_stats.txt" "$OUTPUTDIR"

log "INFO" "Pipeline completed successfully"

