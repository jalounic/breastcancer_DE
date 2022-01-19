#!/bin/bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=featC
#SBATCH --time=13:00:00
#SBATCH --mail-user=james.ackermann@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=./featC.out
#SBATCH --error=./featC.err

module add UHTS/Analysis/subread/2.0.1;

loc=/data/courses/rnaseq/breastcancer_de/jacke_workspace/
annot=${loc}/reference/Homo_sapiens.GRCh38.104.gtf
groups=("NonTNBC" "TNBC" "Normal")


featureCounts -a $annot -o  ${loc}fC_out ${loc}/maps/*sorted.bam
