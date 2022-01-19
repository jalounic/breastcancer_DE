#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=samt_index
#SBATCH --time=1:00:00
#SBATCH --mail-user=james.ackermann@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=./samtools_index.out
#SBATCH --error=./samtools_index.err

module add UHTS/Analysis/samtools/1.10

loc=/data/courses/rnaseq/breastcancer_de/jacke_workspace/maps/
groups=("NonTNBC" "TNBC" "Normal")

for group in ${groups[@]}; do # fixed the loop!!
	for i in 1 2 3; do
		samtools index ${loc}${group}${i}.sorted.bam
	done
done
