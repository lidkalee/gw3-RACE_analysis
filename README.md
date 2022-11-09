# Pipeline to analyze Tail-seq data with an emphasis on 3' mRNA tails identification
Script to analysis of TAIL-seq data from fastq files.


### 0. PrepareSTAR index
STAR --runMode genomeGenerate --genomeDir genome/ --genomeFastaFiles genome/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa


### 1. test.sh
bash test_script.sh -i Wild_type_clone2_R1_001.fastq -I Wild_type_clone2_R2_001.fastq -o output_Wild_type_clone2_20220830


bash joining_R1R2.sh -i Wild_type_clone2_R1_001.fastq -I Wild_type_clone2_R2_001.fastq -a genome/annotation_6k_clean.bed 


### 2. Joining 
  *run the script in the directory with the output)
bash ../joining_R1R2.sh -i R1_Aligned.sortedByCoord.out.bam -I R2_Aligned.out.bam -a ../genome/annotation_6k_clean.bed


### 3.bigwigs 
  * run the script in the directory with the output
  * R1 input: bed, R2 input bam!
bash ../R1_bed_bigWig.sh -i R1_Aligned.sortedByCoord.out.bam_sorted_unique.bed -g ../genome/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.fai 

bash ../R2_bed_bigWig.sh -i R2_Aligned.out.bam  -g ../genome/Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.fai 

### 4. Analysis of tails using Python script
