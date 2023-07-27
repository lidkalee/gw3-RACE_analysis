#! /bin/bash
#######################################################################
# this script filters fastq and makes genome alingment for R1 and R1 de
# rived from gwRACE.
# Usage: bash test_script -i R1 -I R2 -o output_folder_to_be_created
#######################################################################
# this bit sets up user input variables
while getopts i:I:o: flag
do
    case "${flag}" in
        i) R1=${OPTARG};;
        I) R2=${OPTARG};;
	o) OUT=${OPTARG};;
    esac
done
mkdir $OUT
#######################################################################
# create report file
touch $OUT'/report_f.txt'
#######################################################################
# report read number in input files
echo "$R1 read number";
grep "^@" $R1 | wc -l 
echo "$R2 read number";
grep "^@" $R2 | wc -l
printf "input reads\t%s\n" $(grep -c "^@" $R2 | bc) >> $OUT'/report_f.txt'
#######################################################################
# filtering of R2 using adapter sequence
grep -E -A 2 -B 1 --no-group-separator '^[[:alpha:]]{6}GTCAG' $R2 > $R2'_FA.fq'
#######################################################################
# report read number after adapter filtering
echo "$R2 read number after adapter filtering"
grep "^@" $R2'_FA.fq' | wc -l
printf "with adapter\t%s\n" $(grep -c "^@" $R2'_FA.fq' | bc) >> $OUT'/report_f.txt'
#######################################################################
# filter R1 based on R2 filtering results
fastq_pair $R2'_FA.fq' $R1
#######################################################################
# remove input and not paired read files
rm *single*
rm *_FA.fq
#######################################################################
# fastp with quality option disabled !!
# (R2 have low qual because of the tails)
fastp -i $R1'.paired.fq' -I $R2'_FA.fq.paired.fq' -o $OUT'/'$R1'.fp.fq' -O $OUT'/'$R2'.fp.fq' -Q  -U --umi_loc=read2 --umi_len=6 --trim_front2=5 -l  70 -h $OUT'/' -j $OUT'/'
#######################################################################
# remove input files
rm *paired*
#######################################################################
# report read number in output files
echo "$R1 after fastp:"
grep "^@" $OUT'/'$R1'.fp.fq' | wc -l
echo "$R2 after fastp:"
grep "^@" $OUT'/'$R2'.fp.fq' | wc -l
printf "%s after fastp\t%s\n" $R1 $(grep -c "^@" $OUT'/'$R1'.fp.fq' | bc) >> $OUT'/report_f.txt'
########################################################################
# align R1
STAR --genomeDir genome/ --readFilesIn $OUT'/'$R1'.fp.fq' --alignIntronMax 1000 --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1  --outFileNamePrefix $OUT'/R1_'
# align R2
STAR --genomeDir genome/ --readFilesIn $OUT'/'$R2'.fp.fq' --alignIntronMax 1000 --outSAMtype BAM Unsorted --outSAMmultNmax 1 --outSAMunmapped Within --outFileNamePrefix $OUT'/R2_'
grep "Uniquely mapped reads number" $OUT'/R1_Log.final.out' >> $OUT'/report_f.txt'
grep "Number of reads mapped to multiple loci" $OUT'/R1_Log.final.out' >> $OUT'/report_f.txt'
grep "Number of reads unmapped: too short" $OUT'/R1_Log.final.out' >> $OUT'/report_f.txt'
grep "Number of reads unmapped: other" $OUT'/R1_Log.final.out' >> $OUT'/report_f.txt'
printf "%s after fastp\t%s\n" $R2 $(grep -c "^@" $OUT'/'$R2'.fp.fq' | bc) >> $OUT'/report_f.txt'
grep "Uniquely mapped reads number" $OUT'/R2_Log.final.out' >> $OUT'/report_f.txt'
grep "Number of reads mapped to multiple loci" $OUT'/R2_Log.final.out' >> $OUT'/report_f.txt'
grep "Number of reads unmapped: too short" $OUT'/R2_Log.final.out' >> $OUT'/report_f.txt'
grep "Number of reads unmapped: other" $OUT'/R2_Log.final.out' >> $OUT'/report_f.txt'
# as a result I have bam sorted for R1 and unsorted for R2...
