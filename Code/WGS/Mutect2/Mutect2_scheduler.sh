BAM_FOLDER=../GATK

BASE_FOLDER=~/sc-scratch-cc13-multiplex_ihc-pipeline/HB/Kimberly_WGS
DATA_FOLDER=${BASE_FOLDER}/GATK
RESOURCES_FOLDER=${BASE_FOLDER}/../GATK/mm10

REGIONS=${BASE_FOLDER}/../GATK/mm10/call_regions_mm10_without_MaSat.bed


control="P3397_DNA_01_5001WT"
samples="P3397_DNA_03_5001KOWIK1 P3397_DNA_04_5001KOWIK3"

for sample in $samples; do
    sbatch -o GATK_${sample}.out -e GATK_${sample}.err ${BASE_FOLDER}/scripts/Mutect2_runner.sh $DATA_FOLDER $sample $control $REGIONS $RESOURCES_FOLDER
done

