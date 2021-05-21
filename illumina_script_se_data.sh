#!/usr/bin/bash

#Get the file list in the directory
FileList=`ls *.fastq.gz | sed 's/.fastq.gz//'`

for file in $FileList
do

# Input parameters to be provided
sampleName=$file

ref="NC_045512.fasta" # reference genome (for which index have been generated)
bedfile="ref_bed/nCoV-2019.bed" # V3

# some default parameters used in artic illumina pipeline
illuminaKeepLen=20 #Length of illumina reads to keep after primer trimming
illuminaQualThreshold=20 #Sliding window quality threshold for keeping reads after primer trimming (illumina)
mpileupDepth=0 #Mpileup depth for ivar (although undocumented in mpileup, setting to zero removes limit)
ivarFreqThreshold="0.75" # iVar frequency threshold for consensus variant (ivar consensus: -t)
ivarMinDepth=10 # // Minimum coverage depth to call variant (ivar consensus: -m; ivar variants -m)
ivarMinFreqThreshold="0.25" # // iVar frequency threshold to call variant (ivar variants: -t )
ivarMinVariantQuality=20 # // iVar minimum mapQ to call variant (ivar variants: -q)

# make output directory in output
outdirname=`echo $file`
outdir="output_nomask/$outdirname"
mkdir $outdir 

# trim adapter
echo -e "Trimmomatic- trimming adapters and bad quality reads.......\n\n"
TRIM_READS=`echo -e "java -jar /apps/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 16 \
                $sampleName".fastq.gz" $outdir/$sampleName"_trimmed.fastq.gz" \
                ILLUMINACLIP:/apps/Trimmomatic-0.38/adapters/merged_adapters.fa:2:30:30 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:30\n"`
eval $TRIM_READS

# run bwa mem, mapping, filter and sorting, generate stats
echo -e "bwa mem mapping..\n\n"

BWA_MAP=`echo -e "bwa mem -t 32 $ref $outdir/$sampleName"_trimmed.fastq.gz" |\
	samtools sort -@16 -O BAM -o $outdir/$sampleName".sort.bam" -\n\n"`

eval $BWA_MAP


samtools view -@16 -bh -F4 -o $outdir/$sampleName".sort.filt.bam" $outdir/$sampleName".sort.bam" 
samtools index $outdir/$sampleName".sort.filt.bam"  
samtools stat -@16 $outdir/$file".sort.bam" > $file".stat"
	 
# mask primers using ivar
#echo "ivar primer masking.."

####################################################################################################
#sabari included the previous samtools command here##						   #		
##samtools view -F4 -o $outdir/${sampleName}.mapped.bam $outdir/${sampleName}.sorted.bam	   #
##samtools index $outdir/${sampleName}.mapped.bam						   #
####################################################################################################

#ivar trim -e -i $outdir/$outdirname".sort.filt.bam" -b ${bedfile} -m ${illuminaKeepLen} -q ${illuminaQualThreshold} -p $outdir/ivar.out
#samtools sort -o $outdir/$outdirname".primertrimmed.sorted.bam" $outdir/ivar.out.bam
# remove intermediate files
rm -rf $outdir/$outdirname".sort.bam" 
#rm $outdir/ivar.out.bam

# variant calling using ${sampleName}.mapped.primertrimmed.sorted.bam
echo -e "ivar variant calling from samtools mpileup..\n\n"
samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 $outdir/$outdirname".sort.filt.bam" | ivar variants -r ${ref} -m ${ivarMinDepth} -p $outdir/$outdirname".variants" -q ${ivarMinVariantQuality} -t ${ivarMinFreqThreshold}

# use this ${sampleName}.mapped.primertrimmed.sorted.bam file to mark duplicated, and run lofreq analysis
# mark duplicates using picard
#java -jar /apps/picard.jar MarkDuplicates \
#	INPUT=$outdir/$outdirname".primertrimmed.sorted.bam" \
#	OUTPUT=$outdir/$outdirname".primertrimmed.sorted.markdup.bam" \
#	METRICS_FILE=$outdir/$outdirname".dup_metrics.dat" \
#	REMOVE_DUPLICATES=true \
#	ASSUME_SORTED=true \
#	DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
#	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
#	VERBOSITY=ERROR \
#	QUIET=true \
#	VALIDATION_STRINGENCY=LENIENT \
#	MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \
#	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \
#	SORTING_COLLECTION_SIZE_RATIO=0.25 \
#	TAG_DUPLICATE_SET_MEMBERS=false \
#	REMOVE_SEQUENCING_DUPLICATES=false \
#	TAGGING_POLICY=DontTag \
#	CLEAR_DT=true \
#	ADD_PG_TAG_TO_READS=true \
#	PROGRAM_RECORD_ID=MarkDuplicates \
#	PROGRAM_GROUP_NAME=MarkDuplicates \
#	MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 \
#	COMPRESSION_LEVEL=5 \
#	MAX_RECORDS_IN_RAM=500000 \
#	CREATE_INDEX=false \
#	CREATE_MD5_FILE=false \
#	USE_JDK_DEFLATER=false \
#	USE_JDK_INFLATER=false

# lofreq viterbi realignment
#lofreq viterbi -f ${ref} -q 2 $outdir/$outdirname".primertrimmed.sorted.markdup.bam" | samtools sort -@8 -o $outdir/$outdirname".lofreq.realign.bam" -
# lofreq variant calls
#lofreq call --verbose --ref ${ref} --call-indels --min-cov 50 --max-depth 1000000 --min-bq 30 \
#	    --min-alt-bq 30 --def-alt-bq 0 --min-mq 20 --max-mq 255 --min-jq 0 --min-alt-jq 0	\
#	    --def-alt-jq 0 --sig 0.01 --bonf dynamic --no-default-filter \
#	    -o $outdir/$outdirname".lofreq_calls.vcf" $outdir/$outdirname".primertrimmed.sorted.bam" 

# make consensus sequence
samtools mpileup -aa -A -B -d ${mpileupDepth} -Q0 $outdir/$outdirname".sort.filt.bam" | ivar consensus -t ${ivarFreqThreshold} -m ${ivarMinDepth} -n N -p $outdir/$outdirname".consensus"

#snpeff annotations
#snpEff SnpEff4.3_database_for_ncov.snpeffdb -i $outdir/$outdirname".lofreq_calls.vcf" -o $outdir/$outdirname".lofreq_snpeff.vcf" -formatEff -classic -no-downstream -no-intergenic -no-intron -no-upstream \
#	-no-utr -t 8 -stats $outdir/$outdirname".snpeff.dat" 
done


