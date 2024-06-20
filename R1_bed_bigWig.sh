#! /bin/bash
#######################################################################
# this script takes bed file R1 and created bw files for + and - strand 
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
grep "+" $F1 > $F1'_p.bed'
bedtools genomecov -bg -i $F1'_p.bed' -g $F2 > $F1'_p.bg'
rm $F1'_p.bed'
sort -k1,1 -k2,2n $F1'_p.bg' > $F1'_sorted_p.bg'
rm $F1'_p.bg'
bedGraphToBigWig $F1'_sorted_p.bg' $F2 $F1'_p.bw'
rm $F1'_sorted_p.bg'
echo "R2 plus strand bigWig done"
# minus strand
grep "-" $F1 > $F1'_m.bed'
bedtools genomecov -bg -i $F1'_m.bed' -g $F2 > $F1'_m.bg'
rm $F1'_p.bed'
sort -k1,1 -k2,2n $F1'_m.bg' > $F1'_sorted_m.bg'
rm $F1'_m.bg'
bedGraphToBigWig $F1'_sorted_m.bg' $F2 $F1'_m.bw'
rm $F1'_sorted_m.bg'
echo "R2 minus strand bigWig done"
