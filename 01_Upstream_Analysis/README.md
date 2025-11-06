# 🧬 scRNA-seq Upstream Analysis Scripts

This repository contains two Bash scripts designed to streamline the processing of single-cell RNA sequencing (scRNA-seq) data, from downloading raw SRA files to running the **Cell Ranger** pipeline.

---

## ⚠️ Important Notes

- ⚙️ **Local Execution Required:**  
  These scripts **must be executed from a local drive**, not an external hard drive or network-mounted volume.  
  External drives can cause I/O errors during Cell Ranger execution or file movement operations.
  
- 🧾 **Manual Sample List Creation:**  
  The file `samples.txt` **must be created manually** by the user.  
  It should contain one sample name per line, matching the FASTQ filenames (without extensions).

---

## 📁 Overview

### 1. `scRNAseq_SRA_Downloader.sh`
Automates downloading, extracting, compressing, and renaming **FASTQ** files from **SRA** accessions, preparing them for use with **Cell Ranger**.

### 2. `run_cellranger_count.sh`
Executes **Cell Ranger `count`** on selected FASTQ samples, automatically handling resuming, logging, and output organization.

---

## 📘 Understanding FASTQ Files in scRNA-seq Downloads

When downloading scRNA-seq data from the SRA (e.g., 10x Genomics datasets), each sample typically produces **three FASTQ files**: one for the cDNA reads, one for the cell barcodes + UMIs, and one for the index reads.

### File conventions

| File | Typical name | Contains | Description |
|------|----------------|-----------|--------------|
| **R1** | `_1.fastq.gz` | **Read 1** | Usually contains the **cell barcode** and **UMI** (*Unique Molecular Identifier*). These reads are short (≈ 26–28 bases). |
| **R2** | `_2.fastq.gz` | **Read 2** | Contains the **transcript (cDNA)** sequence. This is the longer read (≈ 90–100 bases). |
| **I1** | `_3.fastq.gz` (or sometimes `_I1.fastq.gz`) | **Index 1** | Contains the **sample index** sequence used for demultiplexing by the sequencer. Typically ~8 bases long. |

### Quick summary (example)

| Type | Represents | Typical length | Example in your case |
|------|-------------|----------------|----------------------|
| **R1** | Cell barcode + UMI | ~28 bp | `SRR123456789_1.fastq.gz` |
| **R2** | cDNA / transcript | ~90 bp | `SRR123456789_2.fastq.gz` |
| **I1** | Sample index | ~8 bp | `SRR123456789_3.fastq.gz` |

## 🚀 Script 1: `scrnaseq_data_download_sra_en.sh`

### **Description**
This script automates the download of single-cell RNA-seq datasets from the **SRA** database using `prefetch` and `fasterq-dump`, followed by compression (`pigz`) and renaming of files to the **Cell Ranger** format.


After running this script, the files are automatically renamed to match the **Cell Ranger** naming convention:

```bash
SRR123456789_S1_L001_R1_001.fastq.gz
SRR123456789_S1_L001_R2_001.fastq.gz
SRR123456789_S1_L001_I1_001.fastq.gz
```
This format allows the files to be used directly with the `cellranger count` command without any manual renaming.

---

## ⚠️ Important Notes

- You **must move and run this script from the directory** where you want your downloaded files to be saved.  
  The **SRA Toolkit** (`prefetch` and `fasterq-dump`) creates intermediate folders relative to the current working directory.  
  Example:

  ```bash
  cd /path/to/my/project/
  ./scrnaseq_data_download_sra_en.sh -i accession_list.txt -o fastq_data
  ```

- Make sure you have sufficient disk space — single-cell datasets can be large (10–100 GB per sample).

---

## 🚀 Usage

```bash
./sra_downloader.sh -i accession_list.txt -o output_dir [-t threads]
```

### Parameters

| Option | Description | Required | Default |
|---------|--------------|-----------|----------|
| `-i` | Text file containing SRA accessions (one per line) | ✅ | — |
| `-o` | Output directory for FASTQ files | ✅ | — |
| `-t` | Number of threads for `fasterq-dump` | ❌ | 8 |
| `-h` | Show help message | ❌ | — |

---

## 📦 Example

**Example accession list (`samples.txt`):**
```
SRR12345678
SRR12345679
```

**Run the script from the desired directory:**
```bash
cd /data/scRNAseq/
./sra_downloader.sh -i samples.txt -o fastq_data -t 16
```

**Expected output structure:**
```
fastq_data/
├── SRR12345678_S1_R1_L001_001.fastq.gz
├── SRR12345678_S1_R2_L001_001.fastq.gz
├── SRR12345678_S1_I1_L001_001.fastq.gz
├── SRR12345679_S1_R1_L001_001.fastq.gz
├── SRR12345679_S1_R2_L001_001.fastq.gz
├── SRR12345679_S1_I1_L001_001.fastq.gz
└── reads_stats.txt
```

---

## 🧠 Features

- ✅ Timestamped logging with severity levels (INFO, SUCCESS, ERROR)
- ✅ Per-sample error handling
- ✅ Parallel compression using `pigz`
- ✅ Automatic statistics report via `seqkit stats`
- ✅ Input and thread parameter validation
- ✅ Clean temporary folder management
- ✅ Directory-aware: run it from where you want your data stored

---

## 🧫 Script 2: `run_cellranger_count.sh`

### **Description**
Automates running the **Cell Ranger** `count` command on multiple samples, moving only the `outs` directory to a specified output folder. Includes automatic resume support if interrupted.

### **Usage**
```bash
bash run_cellranger_count.sh -i ./fastqs -o ./results -t 8 -r /path/to/reference -s samples.txt -m 64
```

### **Options**
| Flag | Argument | Description |
|------|-----------|-------------|
| `-i` | `INPUT_DIR` | Directory containing compressed FASTQ files. |
| `-o` | `OUTPUT_DIR` | Destination directory for final results. |
| `-t` | `THREADS` | Number of CPU threads to use. |
| `-r` | `REFERENCE_PATH` | Path to the Cell Ranger reference transcriptome. |
| `-s` | `SAMPLES_FILE` | Text file with sample names (one per line). |
| `-m` | `MEMORY_GB` | Local memory in GB (default: 64). |

### **Features**
- ✅ Automatic resume: skips already processed samples.
- 📜 Logs: creates a log file per sample (`sample.log`).
- 📦 Moves only the `outs` folder to the output directory.
- 🔁 Simple progress bar for tracking execution.

### **Dependencies**
- [Cell Ranger v6](https://www.10xgenomics.com/)
- Bash ≥ 4.0

---

## 📋 Example Workflow

```bash
# 1️⃣ Download and prepare FASTQ files
bash scrnaseq_data_download_sra_en.sh -i sra_list.txt -o ./fastqs -t 16

# 2️⃣ Create the sample list manually
echo -e "SampleA\nSampleB\nSampleC" > samples.txt

# 3️⃣ Run Cell Ranger count for selected samples
bash run_cellranger_count_v3.sh -i ./fastqs -o ./results -t 8 -r /refs/refdata-gex-GRCh38-2020-A -s samples.txt -m 64
```

---

## 🧩 Output Structure

After successful execution:

```
results/
├── SampleA/
│   └── outs/
├── SampleB/
│   └── outs/
└── ...
```

---

## ⚙️ Recommended Environment

- Linux (Ubuntu ≥ 20.04)
- 16+ CPU cores
- ≥ 64 GB RAM
- Stable internet connection (for SRA downloads)
- Local filesystem (not external or network drives)

---

## 🧠 License

This project is released under the **MIT License**.
Feel free to use, modify, and distribute with attribution.
