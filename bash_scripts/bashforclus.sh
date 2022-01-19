#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="bc_fastqc"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00
#SBATCH --mem=8G

source /data/courses/rnaseq/breastcancer_de/jacke_workspace/modules.sh

for read in $(ls /data/courses/rnaseq/breastcancer_de/reads/)
do fastqc $read
done

fastqc -t 2 /data/courses/rnaseq/breastcancer_de/reads/*.fastq.gz

mv /data/courses/rnaseq/breastcancer_de/reads/*.html /data/courses/rnaseq/breastcancer_de/jacke_workspace
mv /data/courses/rnaseq/breastcancer_de/reads/*.zip /data/courses/rnaseq/breastcancer_de/jacke_workspace

