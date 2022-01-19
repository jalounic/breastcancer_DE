#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="hisat2-index"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=16G

module add UHTS/Aligner/hisat/2.2.1

module add UHTS/Aligner/hisat/2.2.1
hisat2-build -f ../reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ../reference/Homo_sapiens.GRCh38.dna.primary_assembly


Homo_sapiens.GRCh38.104.gtf.gz
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

