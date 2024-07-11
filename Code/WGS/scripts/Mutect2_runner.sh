#!/bin/bash
#
#SBATCH --job-name=MUTECT
#
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --time=47:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de


if [ $# -lt 5 ]; then
   echo "Usage:"
   echo "Mutect2_runner.sh BAM_FOLDER TARGET CONTROL ROI RESOURCE_FOLDER"
   exit
fi

BAM_FOLDER=$1
sample=$2
control=$3
REGIONS=$4
RESOURCE_FOLDER=$5

REF_FASTA=$(readlink -f ${RESOURCE_FOLDER}/Mus_musculus.GRCm38.dna.primary_assembly.fa)
DBSNP_FILE=$(readlink -f ${RESOURCE_FOLDER}/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf)
GOLD_INDELS_FILE=$(readlink -f ${RESOURCE_FOLDER}/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf)

export APPTAINER_CACHEDIR=~/sc-scratch-cc13-multiplex_ihc-pipeline/HB/apptainer_cache
     
CWD=$(readlink -f .) # pwd will differ from readlink -f as used later with the full input/output files                                                                                                                                                        
CONTAINER_GATK="apptainer exec --containall --bind ${CWD},$(readlink -f ${RESOURCE_FOLDER}),$(readlink -f $BAM_FOLDER) docker://broadinstitute/gatk:4.2.6.1 gatk "                                                                                                                       
echo $CONTAINER_GATK                                                                                                                                                                                                                                          

READ_MODEL=$(readlink -f ${sample}_f1r2.tar.gz)
OUTFILE_1=$(readlink -f ${sample}.vcf.gz)
LOG_1=$(readlink -f ${sample}.log)

TMP_FOLDER=$(readlink -f .)/tmp  
if ! [ -e $TMP_FOLDER ]; then mkdir $TMP_FOLDER; fi 

# --germline-resource germline_m2.vcf.gz does not work without AF field :(
$CONTAINER_GATK Mutect2 \
--tmp-dir $TMP_FOLDER \
-R $REF_FASTA \
-I $(readlink -f ${BAM_FOLDER}/${sample}_rg_dedup_recal.bam) \
-tumor $sample \
-I $(readlink -f ${BAM_FOLDER}/${control}_rg_dedup_recal.bam) \
-normal $control \
--disable-adaptive-pruning \
--f1r2-tar-gz $READ_MODEL \
-O $OUTFILE_1 &> $LOG_1; 
#-L $REGIONS \
#--germline-resource $DBSNP_FILE \

READ_MODEL_2=$(readlink -f ${sample}_read-orientation-model.tar.gz)
$CONTAINER_GATK LearnReadOrientationModel -I $READ_MODEL -O $READ_MODEL_2 --tmp-dir $TMP_FOLDER

#$GATK GetPileupSummaries \
#-I tumor.bam \
#-V chr17_small_exac_common_3_grch38.vcf.gz \
#-L chr17_small_exac_common_3_grch38.vcf.gz \
#-O getpileupsummaries.table

#$GATK CalculateContamination \
#-I getpileupsummaries.table \
#-tumor-segmentation segments.table \
#-O calculatecontamination.table

OUTFILE_2=$(readlink -f ${sample}_filtered.vcf.gz)
$CONTAINER_GATK FilterMutectCalls -V $OUTFILE_1 -R $REF_FASTA --ob-priors $READ_MODEL_2 -O $OUTFILE_2 --tmp-dir $TMP_FOLDER

bcftools view -f .,PASS -o ${sample}_pass.vcf.gz -O z ${sample}_filtered.vcf.gz
tabix -p vcf ${sample}_pass.vcf.gz

