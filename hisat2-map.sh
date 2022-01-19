#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=genome_indexing
#SBATCH --time=20:00:00
#SBATCH --mail-user=maelle.wannier@unifr.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/home/mwannier/output_%j.o
#SBATCH --error=/home/mwannier/error_%j.e

module add UHTS/Aligner/hisat/2.2.1

mkdir outMap

for reads in Normal HER2; do
	for i in 1 2 3; do
		touch ./outMap/${reads}${i}_hisat2.sam
		hisat2 -q -p 4 -x ./index_files/genome_index -1 /data/courses/rnaseq/breastcancer_de/reads/${reads}${i}_R1.fastq.gz -2 /data/courses/rnaseq/breastcancer_de/reads/${reads}${i}_R2.fastq.gz -S ./outMap/${reads}${i}_hisat2.sam
	done
done

./reference/index_files 

HER21_R1.fastq.gz  HER23_R2.fastq.gz     NonTNBC3_R1.fastq.gz  Normal2_R2.fastq.gz  TNBC2_R1.fastq.gz
HER21_R2.fastq.gz  NonTNBC1_R1.fastq.gz  NonTNBC3_R2.fastq.gz  Normal3_R1.fastq.gz  TNBC2_R2.fastq.gz
HER22_R1.fastq.gz  NonTNBC1_R2.fastq.gz  Normal1_R1.fastq.gz   Normal3_R2.fastq.gz  TNBC3_R1.fastq.gz
HER22_R2.fastq.gz  NonTNBC2_R1.fastq.gz  Normal1_R2.fastq.gz   TNBC1_R1.fastq.gz    TNBC3_R2.fastq.gz
HER23_R1.fastq.gz  NonTNBC2_R2.fastq.gz  Normal2_R1.fastq.gz   TNBC1_R2.fastq.gz


