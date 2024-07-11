#!/bin/bash
BASE_FOLDER=~/sc-scratch-cc13-multiplex_ihc-pipeline/HB/Kimberly_WGS
FASTQ_FOLDER=${BASE_FOLDER}/FASTQ
#GENOME_INDEX=${BASE_FOLDER}/../GATK/b37/human_g1k_v37_decoy.fasta
GENOME_INDEX=${BASE_FOLDER}/../GATK/mm10/BWA_index/GRCm38.fa
SAMTOOLS=$(which samtools)
METAFILE_FOLDER=${BASE_FOLDER}/metadata/
function join { local IFS="$1"; shift; echo "$*"; }

cat $METAFILE_FOLDER/files_by_sample_and_lane.txt | while read sample_line
do
  [ -z "$sample_line" ] && continue
  set $sample_line
  sample=$(echo $1)
  reads1=$(echo $2 | sed -e 's/\,/ /g')
  reads2=$(echo $3 | sed -e 's/\,/ /g')

#[ -e ${sample}_sorted.bam ] && continue

  sbatch ./align_BWA_MEM_slurm_runner.sh $sample $GENOME_INDEX ${FASTQ_FOLDER}/$reads1 ${FASTQ_FOLDER}/$reads2 

done

