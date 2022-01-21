#!/bin/bash

#SBATCH --mail-user=james.ackermann@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name="hisat2-index"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=12G

module add UHTS/Aligner/hisat/2.2.1

hisat2-build -f ../reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa ../reference/Homo_sapiens.GRCh38.dna.primary_assembly

