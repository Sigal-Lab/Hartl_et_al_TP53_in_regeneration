#!/bin/bash
#
#SBATCH --job-name=GATK
#
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --time=47:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de


if [ $# -lt 4 ]; then
   echo "Usage:"
   echo "DNASeq_calling_pipeline_single_sample BAM_FOLDER SAMPLE ROI RESOURCE_FOLDER"
   exit
fi

BAM_FOLDER=$1
sample=$2
ROI_FILE=$3
RESOURCE_FOLDER=$4

BAM_SUFFIX="_sorted_markdup.bam"

#GATK_FOLDER=/data_genome1/SharedSoftware/GATK/GATK_v3.7
#PICARD_FOLDER=/data_genome1/SharedSoftware/Picard
REF_FASTA=$(readlink -f ${RESOURCE_FOLDER}/Mus_musculus.GRCm38.dna.primary_assembly.fa)
#DBSNP_FILE=/data_genome1/References/MusMusculus/Variation/dbSNP/dbSNP_v146_GRCm38/VCF/dbSNP_all.chr.vcf
DBSNP_FILE=$(readlink -f ${RESOURCE_FOLDER}/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf)
GOLD_INDELS_FILE=$(readlink -f ${RESOURCE_FOLDER}/C57BL_6NJ.mgp.v5.indels.dbSNP142.normed.vcf)

# put the temp folder where we know there will be enough space available 
TMP_FOLDER=$(readlink -f .)/tmp
if ! [ -e $TMP_FOLDER ]; then mkdir $TMP_FOLDER; fi

#BASH specific
#http://stackoverflow.com/questions/1527049/bash-join-elements-of-an-array
function join { local IFS="$1"; shift; echo "$*"; }

# This is the final BAM
MERGED_RECAL_BAM=${sample}_rg_dedup_recal.bam
MERGED_RECAL_BAI=${sample}_rg_dedup_recal.bai

CWD=$(readlink -f .) # pwd will differ from readlink -f as used later with the full input/output files
CONTAINER_GATK="apptainer exec --containall --bind ${CWD},$(readlink -f ${RESOURCE_FOLDER}) docker://broadinstitute/gatk:4.2.6.1 gatk "
echo $CONTAINER_GATK

MERGED_BAM=${sample}_merged.bam
if ! [ -e $MERGED_RECAL_BAM ]; then

  if ! [ -e $MERGED_BAM ]; then 
    all_sample_lanes=""
    all_sample_lane_cnt=0

#    echo Adding Read groups...
    ls ${BAM_FOLDER}/${sample}*${BAM_SUFFIX} 
    for sample_lane in ${BAM_FOLDER}/${sample}*${BAM_SUFFIX}; do

	INPUT_BAM_FILE=$sample_lane
	i=$(basename $sample_lane $BAM_SUFFIX)

	# skip per lane files that already have been generated
	if [ -e ${i}_rg.bam -o -e $MERGED_BAM ]; then 
            all_sample_lanes+="INPUT="
            all_sample_lanes+=$INPUT_BAM_FILE
            all_sample_lanes+=" "
            all_sample_lane_cnt+=1
	    continue
	fi

#	run_lane=$(echo $i | awk 'BEGIN{FS="_"} {a=NF-1; b=NF; print $a"_"$b;}')
#	RG=${run_lane}_${sample}
#	LIBRARY=$sample
#	PLATFORM="illumina"
#	RGPU=$run_lane
#	SAMPLE=$sample
#
#	RG_ADDED_SORTED_BAM=${i}_rg.bam
#	if ! [ -e $RG_ADDED_SORTED_BAM ]; then
#	   samtools view -F 8 -u -O BAM $INPUT_BAM_FILE | samtools addreplacerg -r "@RG\tID:${RG}\tPL:${PLATFORM}\tPU:${RGPU}\tLB:${LIBRARY}\tSM:${SAMPLE}" -@ 4 -O BAM -o $RG_ADDED_SORTED_BAM - 
#	fi

	all_sample_lanes+="INPUT="
	all_sample_lanes+=$INPUT_BAM_FILE
	all_sample_lanes+=" "
    all_sample_lane_cnt+=1

   done 

   # now merge per-lane BAMs to per-sample BAMs
   #inputfiles=$(join ' INPUT=' $all_sample_lanes)
   echo $all_sample_lanes
   inputfiles=$all_sample_lanes
   
   # If we have more than one read group/lane/run file per sample
   if [ $all_sample_lane_cnt -gt 1 ]; then
       # If merged bam does not exist yet, merge individual files
       if ! [ -e $MERGED_BAM ]; then
	   #java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MergeSamFiles.jar $inputfiles OUTPUT=$MERGED_BAM TMP_DIR=$TMP_FOLDER 
	   picard -Djava.io.tmpdir=$TMP_FOLDER MergeSamFiles $inputfiles OUTPUT=$MERGED_BAM TMP_DIR=$TMP_FOLDER 

       fi
    
   # single RG/lane/run per sample: change name only
   else
	if ! [ -e $MERGED_BAM ]; then        
	  mv $(echo $all_sample_lanes | sed -e 's/INPUT=//' | sed -e 's/"//g' | sed -e 's/ //g') $MERGED_BAM
	  samtools index $MERGED_BAM
	fi
   fi
  fi
fi
  



#DEDUPPED_BAM=${sample}_rg_dedup.bam
#DEDUP_METRICS_FILE=${sample}.dedup.metrics
#if ! [ -e $DEDUPPED_BAM ]; then
#    #java -Djava.io.tmpdir=$TMP_FOLDER -jar ${PICARD_FOLDER}/MarkDuplicates.jar I=$RG_ADDED_SORTED_BAM O=$DEDUPPED_BAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$DEDUP_METRICS_FILE TMP_DIR=$TMP_FOLDER
#    $CONTAINER_GATK MarkDuplicatesSpark \
#	    -I $(readlink -f $MERGED_BAM) \
#	    -O $(readlink -f $DEDUPPED_BAM) \
#	    --tmp-dir $TMP_FOLDER \
#	    --create-output-bam-index \
#            --conf 'spark.executor.cores=8'
#fi

DEDUPPED_BAM=$MERGED_BAM

RECALIBRATED_BAM=${sample}_rg_dedup_recal.bam
RECAL_TAB_FILE=${sample}_rg_dedup.recal
POST_RECAL_TAB_FILE=${sample}_rg_dedup_post.recal
RECAL_PLOTS=${sample}_rg_dedup_recal.pdf

if ! [ -e $RECALIBRATED_BAM ]; then
   $CONTAINER_GATK BaseRecalibrator \
	    -R $REF_FASTA \
	    -I $(readlink -f $DEDUPPED_BAM) \
	    -O $(readlink -f $RECAL_TAB_FILE) \
	    --known-sites $DBSNP_FILE \
	    --known-sites $GOLD_INDELS_FILE 
	    
   $CONTAINER_GATK ApplyBQSR \
	    -R $REF_FASTA \
	    -I $(readlink -f $DEDUPPED_BAM) \
	    -O $(readlink -f $RECALIBRATED_BAM) \
	    --bqsr-recal-file $(readlink -f $RECAL_TAB_FILE)

#   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $REF_FASTA  -I $REALIGNED_BAM  -knownSites $DBSNP_FILE -knownSites $GOLD_INDELS_FILE -o $RECAL_TAB_FILE
#   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R $REF_FASTA  -I $REALIGNED_BAM  -knownSites $DBSNP_FILE -knownSites $GOLD_INDELS_FILE -BQSR $RECAL_TAB_FILE -o $POST_RECAL_TAB_FILE
#   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF_FASTA   -before $RECAL_TAB_FILE -after $POST_RECAL_TAB_FILE -plots $RECAL_PLOTS
#   java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T PrintReads -nct 8 -I $REALIGNED_BAM -R $REF_FASTA  -BQSR $RECAL_TAB_FILE -o $RECALIBRATED_BAM
fi 

if [ -e $RECALIBRATED_BAM ]; then
    #only remove intermediate files if we reached the last step
    rm $MERGED_BAM 
fi

exit


##########################################################################################################################################
	
   ############################################################################
   # Call variants
   ############################################################################

   OUTPUT_VCF=${sample}.vcf

   if ! [ -e $OUTPUT_VCF ]; then
	   java -Xmx16g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF_FASTA  -nct 1 -nt 1 -glm BOTH -XL hs37d5 --dbsnp $DBSNP_FILE -I $REALIGNED_MERGED_BAM -stand_call_conf 30.0 -o $OUTPUT_VCF 
   fi

   RECAL_FILE=${sample}.recal
   TRANCHES_FILE=${sample}_var_recal.tranches
   RSCRIPT_FILE=${sample}_var_recal.plots.R
   RECAL_SNP_FILE=${sample}_recal_snp_raw_indel.vcf

    if ! [ -e $RECAL_SNP_FILE ]; then
    java -Xmx4g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_FASTA -input $OUTPUT_VCF \
	   -XL hs37d5 \
	   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP_FILE \
	   -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -mode SNP \
	   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
	   --recal_file $RECAL_FILE \
	   --tranches_file $TRANCHES_FILE \
	   --rscript_file $RSCRIPT_FILE

    java -Xmx6g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_FASTA -input $OUTPUT_VCF \
	   -XL hs37d5 -mode SNP --ts_filter_level 99.0 \
	   --recal_file $RECAL_FILE \
	   --tranches_file $TRANCHES_FILE \
	   -o $RECAL_SNP_FILE

    fi

    RECAL_INDEL_FILE=${sample}_indel.recal
    TRANCHES_INDEL_FILE=${sample}_indel_var_recal.tranches
    RSCRIPT_INDEL_FILE=${sample}_indel_var_recal.plots.R
    RECAL_SNP_INDEL_FILE=${sample}_recal_snp_indel.vcf

    if ! [ -e $RECAL_SNP_INDEL_FILE ]; then
	   java -Xmx4g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantRecalibrator -R $REF_FASTA -input $RECAL_SNP_FILE \
		   -XL hs37d5 \
		   -resource:mills,known=true,training=true,truth=true,prior=12.0 $GOLD_INDELS_FILE \
		   -an DP -an FS -an MQRankSum -an ReadPosRankSum -mode INDEL \
		   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 \
		   --recal_file $RECAL_INDEL_FILE \
		   --tranches_file $TRANCHES_INDEL_FILE \
		   --rscript_file $RSCRIPT_INDEL_FILE

	   java -Xmx6g -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T ApplyRecalibration -R $REF_FASTA -input $RECAL_SNP_FILE \
		   -XL hs37d5 -mode INDEL --ts_filter_level 99.0 \
		   --recal_file $RECAL_INDEL_FILE \
		   --tranches_file $TRANCHES_INDEL_FILE \
		   -o $RECAL_SNP_INDEL_FILE
    fi

#   FILTERED_VCF=${sample}_filtered.vcf
   #java -jar ${GATK_FOLDER}/GenomeAnalysisTK.jar -T VariantFiltration -R $REF_FASTA -V $OUTPUT_VCF -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $FILTERED_VCF


