#! /bin/bash
#######################################################################
# this script joins mapped R1 reads with part of sam for R2 read
#######################################################################
# this bit sets up user input variables
while getopts i:I:a: flag
do
    case "${flag}" in
        i) R1=${OPTARG};;
        I) R2=${OPTARG};;
	a) AD=${OPTARG};;
    esac
done
#######################################################################
# create report file
touch report_join.txt
#######################################################################
# working on R1
#######################################################################
# this filters uniquely mapped reads (based on MAPQ value - -q > 5) and
# pipes result to bedtools which reformats bam to bed file.
samtools view -b -q 5 $R1 | bedtools bamtobed > $R1'_unique.bed'
echo "BED file from R1.bam created succesfully"
#######################################################################
# sorting bed and annotation file by coordinates and chromosomes
sort -k1,1 -k2,2n $R1'_unique.bed' > $R1'_sorted_unique.bed'
echo "BED file sorted bycoordinate"
sort -k1,1 -k2,2n $AD > sorted_annotation.bed
echo "annotation file sorted by coordinate"
# removing unsorted file to save space
rm $R1'_unique.bed'
echo "unsorted R1.bed file romoved"
#######################################################################
# intersecting R1 reads with annotation forcing same strendedness to 
# link reads to genes
bedtools intersect -sorted -loj -wb -s -a $R1'_sorted_unique.bed' -b sorted_annotation.bed > merged_R1_annotation.bed
echo "bed file intersected with annotation succesfully"
#  sorting this by read name
sort -k 4 merged_R1_annotation.bed > sorted_merged_R1_annotation.bed
echo "merged file sorted by read name"
# remove unsorted file
rm merged_R1_annotation.bed
echo "unsorted file removed"
#######################################################################
# report number of lines in bed and in merged file in bed should be the 
# same number as number of uniquely aligned reads from STAR
# after merging there will be more lines because some genes overlap thus 
# reads will be repeated
printf "R1_bed file reads\t%s\n" $(wc -l < $R1'_sorted_unique.bed' | bc) >> report_join.txt
printf "R1_bed merged with annotation\t%s\n" $(wc -l < sorted_merged_R1_annotation.bed | bc) >> report_join.txt
#######################################################################
# working our R2
#######################################################################
# this bit takes bam file for R2 (with unmapped reads reported!) and ta
# kes chosen columns (awk bit) to a newly created ahort version of sam
# with read names CIGAR flags and reads sequence...
samtools view $R2 | awk -F"\t" '{ print $1,$2,$3,$4,$5,$6,$10}' OFS='\t' > R2_short.sam
echo "short .sam file created succesfully from R2.bam file"
# sorting by read name
sort -k 1 R2_short.sam > sorted_R2_short.sam
echo "short .sam sorted by read name"
# removing unsorted file
rm R2_short.sam
echo "unsorted file removed"
# report read number in R2 (should be the same as all filtered reads STAR)
printf "R2_short_sam file reads\t%s\n" $(wc -l < sorted_R2_short.sam | bc) >> report_join.txt
########################################################################
# final step - joining two files
########################################################################
join -1 4 -2 1 sorted_merged_R1_annotation.bed sorted_R2_short.sam > merged_output.tab
echo "R1 and R2 files joined succesfully"
rm sorted_merged_R1_annotation.bed
rm sorted_R2_short.sam
printf "merged_output file reads\t%s\n" $(wc -l < merged_output.tab | bc) >> report_join.txt
# remove columns you do not need in the output
awk -F" " '{ print $1,$2,$3,$4,$6,$8,$9,$10,$15,$17,$18}' OFS='\t' merged_output.tab > output_short.tab
echo "unused columns removed succesfully"
rm merged_output.tab
echo "original output removed"

