# nCov19-seq

Custom scripts and pipelines used to analyze Illumina COVIDSeq data generated using NextSeq2000 platform

Dependencies

1. Trimmomatic
   - Trimming the adapters 
2. bwa-mem
   - Mapping 
3. samtools-1.10
   - Processing and filtering the aligned BAM files
4. iVar
   - calling the variants from samtools mpileup
5. snpEff
   - Annotating the variants called by iVar

The pipeline expects above dependencies to be in the user's working path. The pipeline also utilizes an external python script to convert the .tsv files generated from iVar to VCF format. The external pythn scripts has been taken from another github repository-

https://github.com/jts/ncov-tools/blob/master/workflow/scripts/ivar_variants_to_vcf.py
