# ONT Variant Calling Pipeline

## Overview

This pipeline is a modification of the Oxford Nanopore Technologies (ONT) variant calling pipeline and has been adapted for the study: Assessing the readiness of Oxford Nanopore sequencing for clinical genomics applications. It automates the process of ONT variant calling and structural variant detection from  sequencing data. It performs read alignment, quality control, small variant calling (Clair3), structural variant calling (Sniffles), copy number variation (CNV) detection, and short tandem repeat (STR) genotyping.

## Features

- Converts BAM files into FASTQ for alignment with **minimap2**.
- Aligns reads to a reference genome and outputs CRAM format.
- Performs **Quality Control (QC)** using **Alfred** and **Mosdepth**.
- Calls **single nucleotide variants (SNVs) and small insertions/deletions (INDELs)** with **Clair3**.
- Calls **structural variants (SVs)** using **Sniffles**.
- Generates **copy number variation (CNV)** statistics with **Spectre**.
- Infers **biological sex** using **goleft indexcov**.
- Calls **short tandem repeats (STRs)** using **Straglr**.

## Requirements
This pipeline requires the following software installed:
- **Bash** (Linux/macOS compatible)
- [**samtools**](http://www.htslib.org/)
- [**minimap2**](https://github.com/lh3/minimap2)
- [**sentieon**](https://www.sentieon.com/)
- [**Alfred**](https://github.com/tobiasrausch/alfred)
- [**Mosdepth**](https://github.com/brentp/mosdepth)
- [**Sniffles**](https://github.com/fritzsedlazeck/Sniffles)
- [**SURVIVOR**](https://github.com/fritzsedlazeck/SURVIVOR)
- [**Clair3**](https://github.com/HKU-BAL/Clair3)
- [**Spectre**](https://github.com/fritzsedlazeck/Spectre)
- [**goleft**](https://github.com/brentp/goleft)
- [**Straglr**](https://github.com/abishpi/straglr)


## Installation
Clone this repository to get started:
```bash
git clone https://github.com/juditharres/Nanopore4Clinics.git
chmod +x bin/ont_vc_pipeline.sh
```

## Usage
Run the pipeline with the following command:
```bash
./bin/ont_vc_pipeline.sh -s SAMPLE_ID -R RUN_ID -i IN_DIR -o OUT_DIR -t THREADS -r REFERENCE -m CLAIR3_MODEL -b BED_FILE
```

### Required Arguments:
| Argument | Description |
|----------|-------------|
| `-s SAMPLE_ID` | Unique sample identifier |
| `-R RUN_ID` | Sequencing run ID |
| `-i IN_DIR` | Input directory containing BAM files |
| `-o OUT_DIR` | Output directory for results |
| `-t THREADS` | Number of CPU threads to use |
| `-r REFERENCE` | Path to the reference genome (FASTA) |
| `-m CLAIR3_MODEL` | Clair3 model |
| `-b BED_FILE` | BED file with regions |

Example:
```bash
./bin/ont_vc_pipeline.sh -s SAMPLE123 -R RUN001 -i /path/to/input -o /path/to/output -t 16   -r /path/to/reference.fasta -m /path/to/clair3_model -b /path/to/regions.bed
```

## Output Files
The pipeline generates the following key output files:

| File | Description |
|-------------------------------|--------------------------------|
| `OUT_DIR/SAMPLE_ID.sorted.cram` | Aligned reads in CRAM format |
| `OUT_DIR/sample_alfred.json.gz` | Quality control JSON report |
| `OUT_DIR/SAMPLE_ID_sniffles2.vcf.gz` | Structural variant calls |
| `OUT_DIR/SAMPLE_ID_spectre/` | CNV results from Spectre |
| `OUT_DIR/goleft_straglr/gender.txt` | Predicted gender of the sample |
| `OUT_DIR/goleft_straglr/SAMPLE_ID.straglr.vcf` | STR genotyping results |

## Troubleshooting
- **Permission denied?** Run `chmod +x bin/ont_vc_pipeline.sh` before executing the script.
- **Missing dependencies?** Ensure all tools are installed
- **No BAM files found?** Verify that `IN_DIR` contains BAM files.

---
