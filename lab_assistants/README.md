# 🧬 SRA Downloader for scRNA-seq

**Author:** Raúl  
**Version:** 2.0  
**Language:** Bash  

This Bash script automates the download and preprocessing of **single-cell RNA-seq** data from the **NCBI Sequence Read Archive (SRA)** using the **SRA Toolkit**.  
It performs validation, timestamped logging, parallel compression, and optional summary statistics generation with `seqkit`.

---

## ⚠️ Important Notes

- You **must move and run this script from the directory** where you want your downloaded files to be saved.  
  The **SRA Toolkit** (`prefetch` and `fasterq-dump`) creates intermediate folders relative to the current working directory.  
  Example:

  ```bash
  cd /path/to/my/project/
  ./sra_downloader.sh -i accession_list.txt -o fastq_data
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
