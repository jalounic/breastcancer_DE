#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=samtobam
#SBATCH --time=4:00:00
#SBATCH --mail-user=james.ackermann@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=./sambam.out
#SBATCH --error=./sambam.err

module add UHTS/Analysis/samtools/1.10

groups=("NonTNBC" "TNBC" "Normal")

for group in ${groups[@]}; do
	for i in 1 2 3; do
		samtools view -hbS /data/courses/rnaseq/breastcancer_de/jacke_workspace/maps/${group}${i}_hisat2.sam > /data/courses/rnaseq/breastcancer_de/jacke_workspace/maps/${group}${i}_hisat2.bam
	done
done
