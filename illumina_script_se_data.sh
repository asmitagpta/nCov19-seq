#!/usr/bin/bash

start="date +%D"
usage()
{
  echo "Usage: ./illumina_script_se_data.sh -d <Directory of fastq files> -h <help>"
  exit 2
}

#Provide path to directory of input fastq.gz files
while getopts d:h: option
do
 case "${option}"
 in
 d) DIRECTORY=${OPTARG};;
 h|?) usage ;; esac
done


FileList="$(ls $DIRECTORY/Xi* | awk -F'.fastq.gz' '{print $1}' | awk -F'/' '{print $NF}')"

for sampleName in $FileList
do

echo -e "Now processing file $sampleName...\n"

ref="NC_045512.fasta" # reference genome (for which index have been generated)
#bedfile="ref_bed/nCoV-2019.bed" # V3

# some default parameters used in artic illumina pipeline
illuminaKeepLen=20 #Length of illumina reads to keep after primer trimming
illuminaQualThreshold=20 #Sliding window quality threshold for keeping reads after primer trimming (illumina)
mpileupDepth=0 #Mpileup depth for ivar (although undocumented in mpileup, setting to zero removes limit)
ivarFreqThreshold="0.75" # iVar frequency threshold for consensus variant (ivar consensus: -t)
ivarMinDepth=10 # // Minimum coverage depth to call variant (ivar consensus: -m; ivar variants -m)
ivarMinFreqThreshold="0.25" # // iVar frequency threshold to call variant (ivar variants: -t )
ivarMinVariantQuality=20 # // iVar minimum mapQ to call variant (ivar variants: -q)

# make output directory in output
outdirname=`echo $sampleName`
mkdir $outdirname 

# trim adapter
#echo "trim galore!.."
#trim_galore --cores 4 --fastqc $sampleName"_R1_001.fastq.gz" -o $outdir/
echo -e "Trimmomatic- trimming adapters and bad quality reads...\n\n"
TRIM_READS=`echo -e "java -jar /apps/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 16 \
                $DIRECTORY/$sampleName".fastq.gz" $outdirname/$sampleName"_trimmed.fq.gz" \
                ILLUMINACLIP:/apps/Trimmomatic-0.38/adapters/merged_adapters.fa:2:30:30 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:30\n\n"`
eval $TRIM_READS

# run bwa mem, mapping, filter and sorting, generate stats
echo -e "bwa mem mapping..\n\n"

BWA_MAP=`echo -e "bwa mem -t 32 $ref $outdirname/$sampleName"_trimmed.fq.gz" | samtools sort -@16 -O BAM -o $outdirname/$sampleName".sort.bam" -\n\n"`

eval $BWA_MAP


samtools view -@16 -bh -F4 -o $outdirname/$sampleName".sort.filt.bam" $outdirname/$sampleName".sort.bam" 
samtools index $outdirname/$sampleName".sort.filt.bam"  
samtools stat -@16 $outdirname/$sampleName".sort.filt.bam" > $outdirname/$sampleName".stat"
	 
# mask primers using ivar
#echo "ivar primer masking.."
#ivar trim -e -i $outdir/$outdirname".sort.filt.bam" -b ${bedfile} -m ${illuminaKeepLen} -q ${illuminaQualThreshold} -p $outdir/ivar.out
#samtools sort -o $outdir/$outdirname".primertrimmed.sorted.bam" $outdir/ivar.out.bam

# remove intermediate files
rm -rf $outdirname/$sampleName".sort.bam"

# make consensus sequence
samtools mpileup -aa -A -B -d ${mpileupDepth} -Q0 $outdirname/$sampleName".sort.filt.bam" |  \
ivar consensus -t ${ivarFreqThreshold} -m ${ivarMinDepth} -n N -p $outdirname/$sampleName".consensus"


# variant calling using iVar
echo -e "ivar variant calling from samtools mpileup..\n\n"
VAR_CALL=`echo "samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 $outdirname/$sampleName".sort.filt.bam" | \
	ivar variants -r ${ref} -m ${ivarMinDepth} -p $outdirname/$sampleName".variants" -q
${ivarMinVariantQuality} -t ${ivarMinFreqThreshold}"`

eval $VAR_CALL

#convert ivar .tsv to .vcf 
echo -e "convert .tsv to .vcf (iVar)\n\n"
TSVTOVCF=`echo -e python3 ivar_variants_to_vcf.py $outdirname/$sampleName".variants.tsv" $outdirname/$sampleName".variants.vcf"` 

eval $TSVTOVCF

#annotate VCF file using snpEff
echo -e "snpEff annotation...\n\n"
java -jar ../../snpEff/snpEff.jar ncov $outdirname/$sampleName".variants.vcf"  \
        -formatEff -classic -no-downstream -no-intergenic -no-intron -no-upstream -no-utr \
        -stats $outdirname/$sampleName".ivar.dat" > $outdirname/$sampleName".ivar_ann.vcf"

done
