# RNA_Seq_Reference_Based_Sequencing
# Drosophila melanogaster RNA-seq Analysis Pipeline

A complete RNA-seq analysis pipeline for studying tumor-free (w1118) vs tumor-bearing (MARCM RasV12, scrib1) Drosophila larvae, from raw data download to differential expression analysis.

## Dataset Information
- **SRA Accession**: [SRP324433](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP324433)
- **Study**: RNAseq from cuticles of tumor-free and tumor-bearing Drosophila larvae
- **Samples**:
  - Control (w1118): SRR14850651, SRR14850663
  - Treated (MARCM RasV12, scrib1): SRR14850697, SRR14850698

## Pipeline Overview

<img width="800" height="690" alt="image" src="https://github.com/user-attachments/assets/a672688e-2db4-4e0e-bec8-204f61b61282" />

1. Data acquisition from SRA/ENA
2. Quality control with FastQC
3. Adapter trimming with fastp
4. Alignment using:
   - STAR (splice-aware)
   - HISAT2 (alternative)
5. Quantification and differential expression with Cufflinks/Cuffdiff

## Requirements

### Software
- [SRA Toolkit](https://github.com/ncbi/sra-tools) (for direct SRA downloads)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [fastp](https://github.com/OpenGene/fastp)
- [STAR](https://github.com/alexdobin/STAR)
- [HISAT2](http://daehwankimlab.github.io/hisat2/)
- [samtools](http://www.htslib.org/)
- [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)

### Reference Data
- Drosophila melanogaster (dm6) genome:
  ```bash
  wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
  gunzip dm6.fa.gz
