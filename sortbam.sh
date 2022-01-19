#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --job-name=sort_bam
#SBATCH --time=6:00:00
#SBATCH --mail-user=james.ackermann@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=./sortbam.out
#SBATCH --error=./sortbam.err

module add UHTS/Analysis/samtools/1.10

loc=/data/courses/rnaseq/breastcancer_de/jacke_workspace/maps/
group=NonTNBC

for i in 1 2 3; do
	samtools sort -m 25G -@ 4 -o ${loc}${group}${i}.sorted.bam -T temp ${loc}${group}${i}_hisat2.bam
done
