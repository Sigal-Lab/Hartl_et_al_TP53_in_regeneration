#!/bin/bash
#
#SBATCH --job-name=BWA_MEM
#SBATCH --output=slurm_%x__%j.out
#
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --time=47:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de


#!/bin/bash
SAMTOOLS=$(which samtools)

sample=$1
GENOME_INDEX=$2
reads1=$3
reads2=$4

TMP_FOLDER=./work_tmp/${sample}
mkdir -p $TMP_FOLDER
TMP_READS1_FIFO=${TMP_FOLDER}/tmp_reads1
TMP_READS2_FIFO=${TMP_FOLDER}/tmp_reads2

[ ! -p tmp_reads1 ] && mkfifo ${TMP_READS1_FIFO}
[ ! -p tmp_reads2 ] && mkfifo ${TMP_READS2_FIFO}

#[ -e ${sample}_sorted.bam ] && continue

  zcat $reads1 > ${TMP_READS1_FIFO} &
  zcat $reads2 > ${TMP_READS2_FIFO} &

  file_wo_ext=$(basename $reads1 .fastq.gz) 
  # CHECK IF THIS IS APPROPRIATE FOR THE CURRENT FILE NAMES
  #run_lane=$(echo $file_wo_ext | awk 'BEGIN{FS="_"} {a=NF-2; b=NF; print $a"_"$b;}')
  run_lane=$(echo $file_wo_ext | awk 'BEGIN{FS="_"} {a=NF-2; print $a;}') # this data does not seem to have different runs/FC, so just use the lane
  library=$(echo $file_wo_ext | awk 'BEGIN{FS="_"} {j = $1;  for ( i=1; i<(NF-4); ) {j = j "_" $(++i)};  print j;}')

  # CHECK IF THIS IS APPROPRIATE FOR THE CURRENT FILE NAMES, TOO
  RG=${library}_${run_lane}
  LIBRARY=$library
  PLATFORM="illumina"
  RGPU=$run_lane
  SAMPLE=$library

  #bwa mem -t 16 $GENOME_INDEX ${TMP_READS1_FIFO} ${TMP_READS2_FIFO} | $SAMTOOLS view -u -b - | $SAMTOOLS collate -O -u - | $SAMTOOLS fixmate -m -u - ${sample}_fixmate.bam
  bwa mem -t 16 $GENOME_INDEX ${TMP_READS1_FIFO} ${TMP_READS2_FIFO} \
	  | $SAMTOOLS view -u -b - \
	  | $SAMTOOLS addreplacerg -r "@RG\tID:${RG}\tPL:${PLATFORM}\tPU:${RGPU}\tLB:${LIBRARY}\tSM:${SAMPLE}" -u - \
	  | $SAMTOOLS collate -T ${TMP_FOLDER}/COLLATE -O -u - \
	  | $SAMTOOLS fixmate -m -u - - \
	  | $SAMTOOLS sort -@ 4 -T ${TMP_FOLDER}/CHUNK -O bam -o ${sample}_sorted.bam
  #$SAMTOOLS sort -@ 4 -T ${TMP_FOLDER}/CHUNK -O bam -o ${sample}_sorted.bam ${sample}_fixmate.bam
  $SAMTOOLS index -b ${sample}_sorted.bam
  $SAMTOOLS markdup -f ${sample}_markdup.stats ${sample}_sorted.bam ${sample}_sorted_markdup.bam
  $SAMTOOLS index -b ${sample}_sorted_markdup.bam


  if [ -e ${sample}_sorted_markdup.bam ]; then 
	  rm ${sample}_sorted.bam*
  fi	  


rm ${TMP_READS1_FIFO}
rm ${TMP_READS2_FIFO}

