#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=genome_indexing
#SBATCH --time=20:00:00
#SBATCH --mail-user=james.ackermann@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

module add UHTS/Aligner/hisat/2.2.1

mkdir maps

for reads in Normal NonTNBC TNBC; do
	for i in 1 2 3; do
		touch ./maps/${reads}${i}_hisat2.sam
		hisat2 -q -p 4 -x /data/courses/rnaseq/breastcancer_de/jacke_workspace/reference/index_files/Homo_sapiens.GRCh38.dna.primary_assembly -1 /data/courses/rnaseq/breastcancer_de/reads/${reads}${i}_R1.fastq.gz -2 /data/courses/rnaseq/breastcancer_de/reads/${reads}${i}_R2.fastq.gz -S ./maps/${reads}${i}_hisat2.sam
	done
done
