#!/bin/bash
#
#SBATCH --job-name=FASTQC
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=47:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

for f in ../FASTQ/*.fastq.gz; do fastqc -t 10 -o . $f; done


