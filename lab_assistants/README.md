# 🧬 SRA Downloader for scRNA-seq

**Author:** Raúl  
**Version:** 2.0  
**Language:** Bash  

This Bash script automates the download and preprocessing of **single-cell RNA-seq** data from the **NCBI Sequence Read Archive (SRA)** using the **SRA Toolkit**.  
It performs validation, timestamped logging, parallel compression, and optional summary statistics generation with `seqkit`.

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

## 🧩 Dependencies

Make sure the following tools are installed:

```bash
conda install -c bioconda sra-tools pigz seqkit
```

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

## 🧪 Recommended Environment

- Linux (Ubuntu 22.04 / Debian 12)
- Bash ≥ 4.0
- SRA Toolkit --> 3.1.1
- pigz --> v2.8
- seqkit --> v2.10.1

---

## 📜 License

MIT License © 2025 Raúl
