# 🧬 scRNA-seq Automation Scripts

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

### 2. `run_cellranger_count_v3.sh`
Executes **Cell Ranger `count`** on selected FASTQ samples, automatically handling resuming, logging, and output organization.

---

## 🚀 Script 1: `scRNAseq_SRA_Downloader.sh`

### **Description**
This script automates the download of single-cell RNA-seq datasets from the **SRA** database using `prefetch` and `fasterq-dump`, followed by compression (`pigz`) and renaming of files to the **Cell Ranger** format.

### **Usage**
```bash
bash scRNAseq_SRA_Downloader.sh -i accession_list.txt -o ./fastqs -t 16
```

### **Options**
| Flag | Argument | Description |
|------|-----------|-------------|
| `-i` | `ACCESSION_LIST` | File containing SRR accessions (one per line). |
| `-o` | `OUTPUT_DIR` | Output directory for FASTQ files. |
| `-t` | `THREADS` | Number of threads (default: 8, max: 32). |
| `-h` | – | Show help message. |

### **Main Steps**
1. **Download**: Uses `prefetch` to fetch `.sra` files.
2. **FASTQ Extraction**: Converts `.sra` to `.fastq` with `fasterq-dump`.
3. **Compression**: Compresses output FASTQ files using `pigz`.
4. **Statistics**: Generates a summary using `seqkit stats`.
5. **Renaming**: Converts filenames to Cell Ranger-compatible convention:  
   `sample_S1_L001_R1_001.fastq.gz`, `sample_S1_L001_R2_001.fastq.gz`, `sample_S1_I1_L001_001.fastq.gz`.

### **Dependencies**
- [`sra-tools`](https://github.com/ncbi/sra-tools) (`prefetch`, `fasterq-dump`)
- [`pigz`](https://zlib.net/pigz/)
- [`seqkit`](https://bioinf.shenwei.me/seqkit/)

---

## 🧫 Script 2: `run_cellranger_count_v3.sh`

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
bash scRNAseq_SRA_Downloader.sh -i sra_list.txt -o ./fastqs -t 16

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
---

## 🧠 License

This project is released under the **MIT License**.
Feel free to use, modify, and distribute with attribution.
