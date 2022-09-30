#! /bin/bash
#######################################################################
# this script takes bam file R2 and created bw files for + and - strand
# and sorted bed with R2 reads (uniquely aligned)
# -i input bed;; -o genome sizes (fai file)
#######################################################################
# this bit sets up user input variables
while getopts i:g: flag
do
    case "${flag}" in
        i) F1=${OPTARG};;
        g) F2=${OPTARG};;
    esac
done
# plus strand
samtools view -b -q 5 $F1 | bedtools bamtobed > $F1'.bed'
sort -k1,2 -k2,2n $F1'.bed' > $F1'_sorted.bed'
grep "+" $F1'_sorted.bed' > $F1'_p.bed'
bedtools genomecov -bg -i $F1'_p.bed' -g $F2 > $F1'_p.bg'
rm $F1'_p.bed'
sort -k1,1 -k2,2n $F1'_p.bg' > $F1'_sorted_p.bg'
rm $F1'_p.bg'
bedGraphToBigWig $F1'_sorted_p.bg' $F2 $F1'_p.bw'
rm $F1'_sorted_p.bg'
echo "R2 plus strand bigWig done"
# minus strand
grep "-" $F1'_sorted.bed' > $F1'_m.bed'
bedtools genomecov -bg -i $F1'_m.bed' -g $F2 > $F1'_m.bg'
rm $F1'_m.bed'
sort -k1,1 -k2,2n $F1'_m.bg' > $F1'_sorted_m.bg'
rm $F1'_m.bg'
bedGraphToBigWig $F1'_sorted_m.bg' $F2 $F1'_m.bw'
rm $F1'_sorted_m.bg'
echo "R2 minus strand bigWig done"
