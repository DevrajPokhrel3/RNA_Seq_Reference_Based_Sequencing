#Downloading Reference genome from SRA database
#Aim = *RNAseq from cuticles of tumour-free (w1118) and tumour-bearing Drosophila larvae (MARCM RasV12, scrib1) tumour-bearing Drosophila larvae*

#control1
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/051/SRR14850651/SRR14850651_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/051/SRR14850651/SRR14850651_2.fastq.gz

#control2
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/063/SRR14850663/SRR14850663_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/063/SRR14850663/SRR14850663_2.fastq.gz

#Treated1
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/097/SRR14850697/SRR14850697_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/097/SRR14850697/SRR14850697_2.fastq.gz

#Treated2
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/098/SRR14850698/SRR14850698_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/098/SRR14850698/SRR14850698_2.fastq.gz

#Extraction of samples
gunzip *.gz

#Taking small dataset
head -1000000 SRR14850651_1.fastq > control1_1.fastq
head -1000000 SRR14850651_2.fastq > control1_2.fastq

head -1000000 SRR14850663_1.fastq > control2_1.fastq
head -1000000 SRR14850663_2.fastq > control2_2.fastq

head -1000000 SRR14850697_1.fastq > treated1_1.fastq
head -1000000 SRR14850697_2.fastq > treated1_2.fastq

head -1000000 SRR14850698_1.fastq > treated2_1.fastq
head -1000000 SRR14850698_2.fastq > treated2_2.fastq

##Quality control
fastqc *.fastq

#make adapter file using adapter sequence from fastq report
touch adapter_c1.fasta    # adapters from (control1_1.fastq & control1_2.fastq)
touch adapter_c2.fasta    # adapters from (control2_1.fastq & control2_2.fastq)

touch adapter_t1.fasta    # adapters from (treated1_1.fastq & treated1_2.fastq)
touch adapter_t2.fasta    # adapters from (treated2_1.fastq & treated2_2.fastq)

##Removing adapter sequences from DNA sequencing data USING fastp (install)
sudo apt-get install fastp

##Trimmming and cleaning using fastp
fastp -i control1_1.fastq -o Trim_control1_1.fastq -I control1_2.fastq -O Trim_control1_2.fastq --adapter_fasta adapter_c1.fasta
fastp -i control2_1.fastq -o Trim_control2_1.fastq -I control2_2.fastq -O Trim_control2_2.fastq --adapter_fasta adapter_c2.fasta

fastp -i treated1_1.fastq -o Trim_treated1_1.fastq -I treated1_2.fastq -O Trim_treated1_2.fastq --adapter_fasta adapter_t1.fasta
fastp -i treated2_1.fastq -o Trim_treated2_1.fastq -I treated2_2.fastq -O Trim_treated2_2.fastq --adapter_fasta adapter_t2.fasta

#again check the quality
fastqc Trim_control1_1.fastq Trim_control1_2.fastq
fastqc Trim_control2_1.fastq Trim_control2_2.fastq

fastqc Trim_treated1_1.fastq Trim_treated1_2.fastq
fastqc Trim_treated2_1.fastq Trim_treated2_2.fastq

# Downloading Reference genome and GTF file
1. create a folder (dm6) and 2 files (dm6.fa.gz & dm6.gtf )
    #Refrence genome from UCSC 
    --> wget -c https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz

    #Unzip the fa.gz file in the terminal
    gunzip dm6.fa.gz
    
# download GTF file 
    UCSC browser --> Tools --> table browser -->
    clade insect
    genome D. melanogaster
    assembly current assembly
    groups genes and gene pridiction 
    track NCBI RefSeq
    region genome
    output-format GTF
    output-filename dm6.gtf
    file-return-type plain text --> get output

#Mapping using STAR tool  -- Run this in RNAseq folder result will generate in dm6 folder
    a. Building genome index using GTF file
      /mnt/d/DrOmicsClass/TOOLS/STAR/source/STAR --runThreadN 3 --runMode genomeGenerate --genomeDir dm6 --genomeFastaFiles /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.fa --sjdbGTFfile /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.gtf --sjdbOverhang 35 --genomeSAindexNbases 12

    b. Map paired-end reads to genome
      /mnt/d/DrOmicsClass/TOOLS/STAR/source/STAR --runThreadN 3 --readFilesIn Trim_control1_1.fastq Trim_control1_2.fastq --genomeDir dm6 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sorted_control1  --outSAMunmapped Within
      /mnt/d/DrOmicsClass/TOOLS/STAR/source/STAR --runThreadN 3 --readFilesIn Trim_control2_1.fastq Trim_control2_2.fastq --genomeDir dm6 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sorted_control1  --outSAMunmapped Within

      /mnt/d/DrOmicsClass/TOOLS/STAR/source/STAR --runThreadN 3 --readFilesIn Trim_treated1_1.fastq Trim_treated1_2.fastq --genomeDir dm6 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sorted_control1  --outSAMunmapped Within
      /mnt/d/DrOmicsClass/TOOLS/STAR/source/STAR --runThreadN 3 --readFilesIn Trim_treated2_1.fastq Trim_treated2_2.fastq --genomeDir dm6 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix sorted_control1  --outSAMunmapped Within

#Mapping using Hisat2 tool
    sudo apt-get intall hisat2
    
    1. creation of indexed genome
    hisat2-build dm6.fa dm6

    2. mapping of read to indexed genome
    hisat2 --dta-cufflinks -x /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6 -1 Trim_control1_1.fastq -2 Trim_control1_2.fastq > mapped_control1.bam
    hisat2 --dta-cufflinks -x /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6 -1 Trim_control2_1.fastq -2 Trim_control2_2.fastq > mapped_control2.bam

    hisat2 --dta-cufflinks -x /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6 -1 Trim_treated1_1.fastq -2 Trim_treated1_2.fastq > mapped_treated1.bam
    hisat2 --dta-cufflinks -x /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6 -1 Trim_treated2_1.fastq -2 Trim_treated2_2.fastq > mapped_treated2.bam

    3. samtool = conversion of BAM into sorted BAM
    samtool sort mapped_control1.bam > sorted_control1.bam
    samtool sort mapped_control2.bam > sorted_control2.bam

    samtool sort mapped_treated1.bam > sorted_treated1.bam
    samtool sort mapped_treated2.bam > sorted_treated2.bam

#To normalize the data and count the number of readsusing cufflinks
#Normalization = remove noise, biases and variability in the data, technical factor = depth and gene length 

    cufflinks -o cufflinks_c1 -p 2 -G /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.gtf sorted_control1.bam
    cufflinks -o cufflinks_c2 -p 2 -G /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.gtf sorted_control2.bam

    cufflinks -o cufflinks_c3 -p 2 -G /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.gtf sorted_treated1.bam
    cufflinks -o cufflinks_c4 -p 2 -G /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.gtf sorted_treated2.bam

#Using the cuffmerge 
    1. create a assembly file (assembly.txt) containt path and folder name of 4 cufflinks files (cufflinks_c1/c2/t1/t2)
        /mnt/d/DrOmicsClass/7RNAseq/demo/cufflinks_c1/transcripts.gtf
        /mnt/d/DrOmicsClass/7RNAseq/demo/cufflinks_c2/transcripts.gtf
        /mnt/d/DrOmicsClass/7RNAseq/demo/cufflinks_t1/transcripts.gtf
        /mnt/d/DrOmicsClass/7RNAseq/demo/cufflinks_t2/transcripts.gtf
    
    2. Merging using the cuffmerge
        cuffmerge -o cuffmerge_output -g /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.gtf -s /mnt/d/DrOmicsClass/7RNAseq/demo/dm6/dm6.fa assembly.txt

#Quantifying the differential gene expressoin using cuffdiff
    # 1st - Decide the analysis plan   - - - treated vs control
        cuffdiff -o cuffdiff_output -p 3 -L control,Treated -u /mnt/d/DrOmicsClass/7RNAseq/demo/cuffmerge_output/merged.gtf sorted_control1.bam,sorted_control2.bam sorted_treated1.bam,sorted_treated2.bam
