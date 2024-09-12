#!/bin/bash

# Script to annotate VCFs with annovar for Hartl et al. (2024)

# Dependencies:
# * bcftools (can be installed through conda: https://anaconda.org/bioconda/bcftools)
# * annovar (read documentation and install here: https://annovar.openbioinformatics.org/ )
# * then install annovar refseq like this: `./annovar/ annovar_script_share.shannotate_variation.pl -buildver mm9 -downdb -webfrom annovar refGene mousedb/`

# INPUT: $1 should be vcf.gz file to annotate

input_file=$1
echo "annotating this file: $input_file"
filename=$(echo $input_file | rev | cut -d "/" -f 1 | rev)
prefix=$(echo $filename | sed "s/.vcf.gz//g")
avinput_prefix=${prefix}_avinput
avoutput_prefix=${prefix}_avoutput

if [ ! -f "$input_file" ] ; then
        echo "ERROR: This input file does not exist: ${input_file}"
	exit
fi

# set work directory
workdir=variant_annotation/${prefix}
mkdir -p $workdir

# make input file
Mutect_usage=$(bcftools view -h $input_file | grep 'MutectVersion')

if [ -n "$Mutect_usage" ]
then
        echo "Mutect usage detected"
	normal_sample=$(bcftools view -h $input_file | grep '^##normal_sample=' | sed 's/##normal_sample=//')
        echo "running annovar on $tumour_sample only and not the normal sample"
        bcftools annotate --remove "INFO/PON" $input_file \
		| bcftools view -s^$normal_sample -Oz \
		> $workdir/input_filt.vcf.gz
        ./annovar/convert2annovar.pl \
                -format vcf4 <(bgzip -cd -@ $(nproc) $workdir/input_filt.vcf.gz) \
                -outfile ${workdir}/${avinput_prefix}.avinput -includeinfo -comment
else
        echo "Mutect usage *NOT* detected"
        ./annovar/convert2annovar.pl \
                -format vcf4 <(bgzip -cd -@ $(nproc) $input_file) \
                -outfile ${workdir}/${avinput_prefix}.avinput -includeinfo -comment \
		-allsample -withfreq
fi
echo "annovar input file written to: ${workdir}/${avinput_prefix}"


# correct input file (convert2annovar.pl adds duplicated columns) and run annovar
echo "running correction of input"
# correct format
grep '^#' ${workdir}/${avinput_prefix}.avinput \
        > ${workdir}/vcf_header.txt

wrongnum=$(cat ${workdir}/${avinput_prefix}.avinput | grep -v '^#' | head -n1 | wc -w)
rightnum=$(cat ${workdir}/${avinput_prefix}.avinput | grep '^#CHROM' | wc -w)
mydiff=$(expr $wrongnum - $rightnum )
mystr=$(seq $(expr $mydiff + 1) $wrongnum | tr '\n' ',' )
mystredit=$(echo ${mystr%?})

cat ${workdir}/${avinput_prefix}.avinput \
        | grep -v '^#' \
        | cut -f $mystredit \
        > ${workdir}/vcf_corrected.txt

cat ${workdir}/vcf_header.txt > ${workdir}/${avinput_prefix}.corrected.avinput
cat ${workdir}/vcf_corrected.txt >> ${workdir}/${avinput_prefix}.corrected.avinput


echo "input file corrected, written to: ${workdir}/${avinput_prefix}.corrected.avinput"

# run Annovar annotation
./annovar/table_annovar.pl \
        ${workdir}/${avinput_prefix}.corrected.avinput \
        /fast/work/groups/ag_sanders/tools/annovar/mousedb \
        --buildver mm10 \
        --out ${workdir}/${avoutput_prefix} --remove \
        --protocol refGene \
        --operation g --otherinfo --nastring . --vcfinput --thread $(nproc) \
	--maxgenethread $(nproc)

# add dbnsfp42a in --protocol for all variant effect predictions

# process output and remove unneeded columsn that arent present for all vars and f up the format
bcftools view ${workdir}/${avoutput_prefix}*.vcf \
        | bcftools sort \
	| bcftools annotate --remove "INFO/INDEL,INFO/VDB,INFO/RPB,INFO/MQB,INFO/MQSB,INFO/BQB,INFO/PV4,INFO/HOB,INFO/ICB,INFO/IDV,INFO/IMF,INFO/MQ0F,INFO/SGB,INFO/AC,INFO/AN" \
        | bcftools view \
		-O z \
		--threads \
		$(nproc) \
		> ${workdir}/${prefix}_annotated.vcf.gz

echo "final output written to ${workdir}/${prefix}_annotated.vcf.gz"
