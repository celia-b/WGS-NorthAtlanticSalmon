#!/bin/bash

## These commands are to be run sequentially on a Unix terminal
## to reproduce the genomic analysis on NGS whole genome data 
## of the project "NORTH ATLANTIC SALMON WHOLE GENOME SEQUENCING ANALYSIS"
## by Celia Burgos Sequeros, Amaya Zaratiegui Pedrosa and Bence Daniel Kovaks.

# Module loads ---------------------------- #
module load plink2/1.90beta6.18
module load gcc/9.3.0
module load intel/perflibs/64/2019_update5
module load proj/7.0.0
module load R/3.6.1
module load sratoolkit/2.10.7
module load perl/5.30.2
module load java/1.8.0
module load fastqc/0.11.2
module load samtools/1.10
module load trimmomatic/0.38
module load bwa/0.7.15
module load bcftools
module load tabix
module load bedtools
module load vcftools

# File downloads ---------------------------- #
# St. John River (SJR)
prefetch –max-size 30000000 SRR9925703 && fastq-dump --split-files SRR9925703.sra 
prefetch –max-size 30000000 SRR9925704 && fastq-dump --split-files SRR9925704.sra 
prefetch –max-size 30000000 SRR9925706 && fastq-dump --split-files SRR9925706.sra 
prefetch –max-size 30000000 SRR9925707 && fastq-dump --split-files SRR9925707.sra 
# Penobscot River (PR)
prefetch –max-size 30000000 SRR9925722 && fastq-dump --split-files SRR9925722.sra 
prefetch –max-size 30000000 SRR9925723 && fastq-dump --split-files SRR9925723.sra 
prefetch –max-size 30000000 SRR9925724 && fastq-dump --split-files SRR9925724.sra 

# Quality Control ---------------------------- #
fastqc SRR9925703_1.fastq -o .
fastqc SRR9925703_2.fastq -o .

fastqc SRR9925704_1.fastq -o .
fastqc SRR9925704_2.fastq -o .

fastqc SRR9925706_1.fastq -o .
fastqc SRR9925706_2.fastq -o .

fastqc SRR9925707_1.fastq -o .
fastqc SRR9925707_2.fastq -o .

fastqc SRR9925722_1.fastq -o .
fastqc SRR9925722_2.fastq -o .

fastqc SRR9925723_1.fastq -o .
fastqc SRR9925723_2.fastq -o .

fastqc SRR9925724_1.fastq -o .
fastqc SRR9925724_2.fastq -o .

# Read trimming w/trimommatic ------------------ #
trimommatic -phred33 SRR9925703_1.fastq SRR9925703_2.fastq SRR9925703_reduced_1_paired.fastq SRR9925703_reduced_1_unpaired.fastq SRR9925703_reduced_2_paired.fastq SRR9925703_reduced_2_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36
trimommatic -phred33 SRR9925704_1.fastq SRR9925704_2.fastq SRR9925704_reduced_1_paired.fastq SRR9925704_reduced_1_unpaired.fastq SRR9925704_reduced_2_paired.fastq SRR9925704_reduced_2_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36
trimommatic -phred33 SRR9925706_1.fastq SRR9925706_2.fastq SRR9925706_reduced_1_paired.fastq SRR9925706_reduced_1_unpaired.fastq SRR9925706_reduced_2_paired.fastq SRR9925706_reduced_2_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36
trimommatic -phred33 SRR9925707_1.fastq SRR9925707_2.fastq SRR9925707_reduced_1_paired.fastq SRR9925707_reduced_1_unpaired.fastq SRR9925707_reduced_2_paired.fastq SRR9925707_reduced_2_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36
trimommatic -phred33 SRR9925722_1.fastq SRR9925722_2.fastq SRR9925722_reduced_1_paired.fastq SRR9925722_reduced_1_unpaired.fastq SRR9925722_reduced_2_paired.fastq SRR9925722_reduced_2_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36
trimommatic -phred33 SRR9925723_1.fastq SRR9925723_2.fastq SRR9925723_reduced_1_paired.fastq SRR9925723_reduced_1_unpaired.fastq SRR9925723_reduced_2_paired.fastq SRR9925723_reduced_2_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36
trimommatic -phred33 SRR9925724_1.fastq SRR9925724_2.fastq SRR9925724_reduced_1_paired.fastq SRR9925724_reduced_1_unpaired.fastq SRR9925724_reduced_2_paired.fastq SRR9925724_reduced_2_unpaired.fastq SLIDINGWINDOW:4:15 MINLEN:36

# Reference Genome Download -------------------- #
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.fna.gz

# Keep only chromosomes, remove scaffolds
grep -n '>' GCF_000233375.1_ICSASG_v2_genomic.fna | head -30
head -28002605 GCF_000233375.1_ICSASG_v2_genomic.fna > Salmo_salar_RFG.fa

# Index genome
bwa index Salmo_salar_RFG.fa


# Read Alignment -------------------------------- #
bwa mem Salmo_salar_RFG.fa SRR9925703_reduced_1_paired.fastq SRR9925703_reduced_2_paired.fastq > SRR9925703.sam
bwa mem Salmo_salar_RFG.fa SRR9925704_reduced_1_paired.fastq SRR9925704_reduced_2_paired.fastq > SRR9925704.sam
bwa mem Salmo_salar_RFG.fa SRR9925706_reduced_1_paired.fastq SRR9925706_reduced_2_paired.fastq > SRR9925706.sam
bwa mem Salmo_salar_RFG.fa SRR9925707_reduced_1_paired.fastq SRR9925707_reduced_2_paired.fastq > SRR9925707.sam
bwa mem Salmo_salar_RFG.fa SRR9925722_reduced_1_paired.fastq SRR9925722_reduced_2_paired.fastq > SRR9925722.sam
bwa mem Salmo_salar_RFG.fa SRR9925723_reduced_1_paired.fastq SRR9925723_reduced_2_paired.fastq > SRR9925723.sam
bwa mem Salmo_salar_RFG.fa SRR9925724_reduced_1_paired.fastq SRR9925724_reduced_2_paired.fastq > SRR9925724.sam

# Convert .sam to .bam --------------------------- #
samtools view -b SRR9925703.sam > SRR9925703.bam
samtools view -b SRR9925704.sam > SRR9925704.bam
samtools view -b SRR9925706.sam > SRR9925706.bam
samtools view -b SRR9925707.sam > SRR9925707.bam
samtools view -b SRR9925722.sam > SRR9925722.bam
samtools view -b SRR9925723.sam > SRR9925724.bam
samtools view -b SRR9925724.sam > SRR9925724.bam

# Sort .bam files -------------------------------- #
samtools sort SRR9925703.bam -o SRR9925703_sorted.bam 
samtools sort SRR9925704.bam -o SRR9925704_sorted.bam
samtools sort SRR9925706.bam -o SRR9925706_sorted.bam
samtools sort SRR9925707.bam -o SRR9925707_sorted.bam
samtools sort SRR9925722.bam -o SRR9925722_sorted.bam
samtools sort SRR9925723.bam -o SRR9925723_sorted.bam
samtools sort SRR9925724.bam -o SRR9925724_sorted.bam

# Get stats --------------------------------------- #
samtools flagstat SRR9925703.bam > SRR9925703_stats.txt
samtools flagstat SRR9925704.bam > SRR9925704_stats.txt
samtools flagstat SRR9925706.bam > SRR9925706_stats.txt
samtools flagstat SRR9925707.bam > SRR9925707_stats.txt
samtools flagstat SRR9925722.bam > SRR9925722_stats.txt
samtools flagstat SRR9925723.bam > SRR9925723_stats.txt
samtools flagstat SRR9925724.bam > SRR9925724_stats.txt

# Filter .bam files (low quality and duplicates) ----- #
samtools view -b -q10 SRR9925703_sorted.bam > SRR9925703_sorted_q10.bam
samtools view -b -q10 SRR9925704_sorted.bam > SRR9925704_sorted_q10.bam
samtools view -b -q10 SRR9925706_sorted.bam > SRR9925706_sorted_q10.bam
samtools view -b -q10 SRR9925707_sorted.bam > SRR9925707_sorted_q10.bam
samtools view -b -q10 SRR9925722_sorted.bam > SRR9925722_sorted_q10.bam
samtools view -b -q10 SRR9925723_sorted.bam > SRR9925723_sorted_q10.bam
samtools view -b -q10 SRR9925724_sorted.bam > SRR9925724_sorted_q10.bam

samtools rmdup SRR9925703_sorted_q10.bam SRR9925703_sorted_q10_rmdup.bam
samtools rmdup SRR9925704_sorted_q10.bam SRR9925704_sorted_q10_rmdup.bam
samtools rmdup SRR9925706_sorted_q10.bam SRR9925706_sorted_q10_rmdup.bam
samtools rmdup SRR9925707_sorted_q10.bam SRR9925707_sorted_q10_rmdup.bam
samtools rmdup SRR9925722_sorted_q10.bam SRR9925722_sorted_q10_rmdup.bam
samtools rmdup SRR9925723_sorted_q10.bam SRR9925723_sorted_q10_rmdup.bam
samtools rmdup SRR9925724_sorted_q10.bam SRR9925724_sorted_q10_rmdup.bam

# Index .bam files ----------------------------------- #
samtools index SRR9925703_sorted_q10_rmdup.bam
samtools index SRR9925704_sorted_q10_rmdup.bam
samtools index SRR9925706_sorted_q10_rmdup.bam
samtools index SRR9925707_sorted_q10_rmdup.bam
samtools index SRR9925722_sorted_q10_rmdup.bam
samtools index SRR9925723_sorted_q10_rmdup.bam
samtools index SRR9925724_sorted_q10_rmdup.bam

# Get stats per chromosome --------------------------- #
samtools idxstats SRR9925703_sorted_q10_rmdup.bam > SRR9925703_chr_stats.txt
samtools idxstats SRR9925704_sorted_q10_rmdup.bam > SRR9925704_chr_stats.txt
samtools idxstats SRR9925706_sorted_q10_rmdup.bam > SRR9925706_chr_stats.txt
samtools idxstats SRR9925707_sorted_q10_rmdup.bam > SRR9925707_chr_stats.txt
samtools idxstats SRR9925722_sorted_q10_rmdup.bam > SRR9925722_chr_stats.txt
samtools idxstats SRR9925723_sorted_q10_rmdup.bam > SRR9925723_chr_stats.txt
samtools idxstats SRR9925724_sorted_q10_rmdup.bam > SRR9925724_chr_stats.txt

# Select chromosomes 11, 12 and 13 -------------------- #
samtools view -b SRR9925703_sorted_q10_rmdup.bam NC_027310.1 NC_027311.1 NC_027312.1 > 03_ss11-12-13.bam
samtools view -b SRR9925704_sorted_q10_rmdup.bam NC_027310.1 NC_027311.1 NC_027312.1 > 04_ss11-12-13.bam
samtools view -b SRR9925706_sorted_q10_rmdup.bam NC_027310.1 NC_027311.1 NC_027312.1 > 06_ss11-12-13.bam
samtools view -b SRR9925707_sorted_q10_rmdup.bam NC_027310.1 NC_027311.1 NC_027312.1 > 07_ss11-12-13.bam
samtools view -b SRR9925722_sorted_q10_rmdup.bam NC_027310.1 NC_027311.1 NC_027312.1 > 22_ss11-12-13.bam
samtools view -b SRR9925723_sorted_q10_rmdup.bam NC_027310.1 NC_027311.1 NC_027312.1 > 23_ss11-12-13.bam
samtools view -b SRR9925724_sorted_q10_rmdup.bam NC_027310.1 NC_027311.1 NC_027312.1 > 24_ss11-12-13.bam

# Label samples --------------------------------------- #
java -jar /services/tools/picard-tools/2.20.2/picard.jar  AddOrReplaceReadGroups I=03_ss11-12-13.bam O=03_ss11-12-13_labeled.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU= NONE RGSM=03 CREATE_INDEX=True
java -jar /services/tools/picard-tools/2.20.2/picard.jar  AddOrReplaceReadGroups I=04_ss11-12-13.bam O=04_ss11-12-13_labeled.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU= NONE RGSM=04 CREATE_INDEX=True
java -jar /services/tools/picard-tools/2.20.2/picard.jar  AddOrReplaceReadGroups I=06_ss11-12-13.bam O=06_ss11-12-13_labeled.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU= NONE RGSM=06 CREATE_INDEX=True
java -jar /services/tools/picard-tools/2.20.2/picard.jar  AddOrReplaceReadGroups I=07_ss11-12-13.bam O=07_ss11-12-13_labeled.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU= NONE RGSM=07 CREATE_INDEX=True
java -jar /services/tools/picard-tools/2.20.2/picard.jar  AddOrReplaceReadGroups I=22_ss11-12-13.bam O=22_ss11-12-13_labeled.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU= NONE RGSM=22 CREATE_INDEX=True
java -jar /services/tools/picard-tools/2.20.2/picard.jar  AddOrReplaceReadGroups I=23_ss11-12-13.bam O=23_ss11-12-13_labeled.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU= NONE RGSM=23 CREATE_INDEX=True
java -jar /services/tools/picard-tools/2.20.2/picard.jar  AddOrReplaceReadGroups I=24_ss11-12-13.bam O=24_ss11-12-13_labeled.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=illumina RGPU= NONE RGSM=24 CREATE_INDEX=True

# Index ----------------------------------------------- #
samtools index 03_ss11-12-13_labeled.bam
samtools index 04_ss11-12-13_labeled.bam
samtools index 06_ss11-12-13_labeled.bam
samtools index 07_ss11-12-13_labeled.bam
samtools index 22_ss11-12-13_labeled.bam
samtools index 23_ss11-12-13_labeled.bam
samtools index 24_ss11-12-13_labeled.bam

## ----------- ##
## SNP CALLING ##
## ----------- ##

# Call variants -------------------------------------- #
samtools mpileup -uf Salmo_salar_RFG_11-12-13.fa 03_ss11-12-13_labeled.bam | bcftools call -m -v -O z - > 03_var.raw.vcf.gz
samtools mpileup -uf Salmo_salar_RFG_11-12-13.fa 04_ss11-12-13_labeled.bam | bcftools call -m -v -O z - > 04_var.raw.vcf.gz
samtools mpileup -uf Salmo_salar_RFG_11-12-13.fa 06_ss11-12-13_labeled.bam | bcftools call -m -v -O z - > 06_var.raw.vcf.gz
samtools mpileup -uf Salmo_salar_RFG_11-12-13.fa 07_ss11-12-13_labeled.bam | bcftools call -m -v -O z - > 07_var.raw.vcf.gz
samtools mpileup -uf Salmo_salar_RFG_11-12-13.fa 22_ss11-12-13_labeled.bam | bcftools call -m -v -O z - > 22_var.raw.vcf.gz
samtools mpileup -uf Salmo_salar_RFG_11-12-13.fa 23_ss11-12-13_labeled.bam | bcftools call -m -v -O z - > 23_var.raw.vcf.gz
samtools mpileup -uf Salmo_salar_RFG_11-12-13.fa 24_ss11-12-13_labeled.bam | bcftools call -m -v -O z - > 24_var.raw.vcf.gz

# Index VCF files ------------------------------------- #
tabix -p vcf 03_var.raw.vcf.gz
tabix -p vcf 04_var.raw.vcf.gz
tabix -p vcf 06_var.raw.vcf.gz
tabix -p vcf 07_var.raw.vcf.gz
tabix -p vcf 22_var.raw.vcf.gz
tabix -p vcf 23_var.raw.vcf.gz
tabix -p vcf 24_var.raw.vcf.gz

# Merge samples --------------------------------------- #
vcf-merge 03_var.raw.vcf.gz 04_var.raw.vcf.gz 06_var.raw.vcf.gz 07_var.raw.vcf.gz 22_var.raw.vcf.gz 23_var.raw.vcf.gz 24_var.raw.vcf.gz | bgzip -c > merged_samples.vcf.gz

# Total number of SNPS detected
zcat merged_samples.vcf.gz | grep -v "#" | wc –l # 1992881 total SNPs

# Remove INDELs and low quality SNPs ------------------ #
vcftools --gzvcf merged_samples.vcf.gz --recode --out merged_samples_filtered --minQ 30 --remove-indels 

# Number of SNPs after filtering
zcat merged_samples_filtered.vcf.gz | grep -v "#" | wc –l # 1127216 SNPs

# Run variant effect predictor ------------------------- #
grep -v "##" merged_samples_filtered.recode.vcf | grep -v "#CHROM" | sed 's/NC_027310.1/ssa11/g' | sed 's/NC_027311.1/ssa12/g' | sed 's/NC_027312.1/ssa13/g'| awk '{print $1, $2, $2, $4"/" $5}' > for_vep

# Convert to .ped and .map ----------------------------- # 
vcftools --vcf SRJ.recode.vcf --plink --out SRJ

# How many SNPs have been mapped?
wc -l merged_samples_filtered.recode.map
1119700 SNPs

# How many samples are we working with?
wc -l merged_samples_filtered.recode.ped 
7 samples

# Convert to .tfam and .tped ---------------------------- #
plink --allow-extra-chr --file SRJ.recode --transpose --recode --out SRJ.recode

# Filter out SNPs with low genotyping rates ------------- #
plink --allow-extra-chr --file merged_samples_filtered.recode --geno 0.25 --recode --out 75

# Merge samples by population --------------------------- #
vcf-merge 03_var.raw.vcf.gz 04_var.raw.vcf.gz 06_var.raw.vcf.gz 07_var.raw.vcf.gz | bgzip -c > SRJ.vcf.gz
vcf-merge 22_var.raw.vcf.gz 23_var.raw.vcf.gz 24_var.raw.vcf.gz | bgzip -c > PR.vcf.gz

# Label the .ped file by population --------------------- #
awk '{if ($1 == "SRR9925703" || $1 == "04" || $1 == "06" || $1 == "SRR9925707") printf "SRJ "; for (i=2; i<NF; i++) printf $i " "; print $NF}' 50.ped > relabeled.ped
awk '{if ($1 != "SRJ") {printf "PR "; for (i=1; i<NF; i++) printf $i " "; print $NF} else {for (i=1; i<NF; i++) printf $i " "; print $NF}}' relabeled.ped > rerelabeled.ped

# Label the .ped file by population using plink --------- #
plink --bfile labeled --update-ids renaming.txt --make-bed --out relabeled

# Label the .map file by chromosome --------------------- #
grep "NC_027310.1" 50.map | awk '{print "11", $2, $3, $4 }' > test.map
grep "NC_027311.1" 50.map | awk '{print "11", $2, $3, $4 }' >> test.map
grep "NC_027312.1" 50.map | awk '{print "11", $2, $3, $4 }' >> test.map


## ------------ ##
## FST ANALYSIS ##
## ------------ ##

# Individual SNPs --------------------------------------- #
plink --autosome-num 29 --allow-no-sex --file test --pheno pheno.txt --fst case-control --out FSTtest

# By windows -------------------------------------------- #
sed '1d' test_filtered.fst | awk '{print $1, $3, $2, $4, $5}' > test_filtered_noheader.fst
for el in $(cat /home/LECTURE_09/chromosome_lenght | awk '{print $1}' ); do python /home/LECTURE_09/Fst.py test_filtered_noheader.fst /home/LECTURE_09/chromosome_lenght $el > $el\_dagrep ; done
cat *dagrep |awk '$2!=0 {$4=$3/$2; print}' | sed 's/_/ /g' | sort -nk1 -nk2 > salmon_fst_windows

# Visualization in R script in repository.

## -------------------- ##
## Runs of Homozygosity ##
## -------------------- ##

plink --autosome-num 29 --file SJR --homozyg --homozyg-het 3 --homozyg-window-snp 500 --out roh_500R

# Visualization in R script in repository.

