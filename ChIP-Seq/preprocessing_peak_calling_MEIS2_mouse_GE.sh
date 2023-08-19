#!/bin/bash


path="path/to/file"
sample="GE_meis_IP"
input="GE_meis_input"

# cutadapt 

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o $path/raw_data/"$sample"_cutadapt.fastq.gz  \
$path/raw_data/meis_GE_IP.fastq.gz   \
--minimum-length 30 --cores=1 &> $path/out/"$sample"_cutadapt_report.txt

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -o $path/raw_data/"$input"_cutadapt.fastq.gz  \
$path/raw_data/meis_GE_input.fastq.gz   \
--minimum-length 30 --cores=1 &> $path/out/"$input"_cutadapt_report.txt

# mapping
bowtie2 -p 6 -x /UCSC/mm10/Sequence/Bowtie2Index/genome \
 -U $path/raw_data/"$sample"_cutadapt.fastq.gz \
 --no-mixed 2> $path/mapped/"$sample"_mapping.log | samtools view -bq 30 - > $path/mapped/"$sample"_mm10.bam

 bowtie2 -p 6 -x /UCSC/mm10/Sequence/Bowtie2Index/genome \
  -U $path/raw_data/"$input"_cutadapt.fastq.gz \
  --no-mixed 2> $path/mapped/"$input"_mapping.log | samtools view -bq 30 - > $path/mapped/"$input"_mm10.bam

# sorting_indexing
samtools sort -O bam -T tmp -@ 4 $path/mapped/"$sample"_mm10.bam  -o $path/mapped/"$sample"_mm10_sorted.bam &&
samtools index $path/mapped/"$sample"_mm10_sorted.bam

samtools sort -O bam -T tmp -@ 4 $path/mapped/"$input"_mm10.bam  -o $path/mapped/"$input"_mm10_sorted.bam &&
samtools index $path/mapped/"$input"_mm10_sorted.bam

# remove duplicates with picardtools
java -jar /jar/picard.jar MarkDuplicates I= $path/mapped/"$sample"_mm10_sorted.bam \
O= $path/mapped/"$sample"_mm10_sorted_rmdup.bam M= $path/out/"$sample"_mm10_sorted_metrics.txt REMOVE_DUPLICATES=true

java -jar /jar/picard.jar MarkDuplicates I= $path/mapped/"$input"_mm10_sorted.bam \
O= $path/mapped/"$input"_mm10_sorted_rmdup.bam M= $path/out/"$input"_mm10_sorted_metrics.txt REMOVE_DUPLICATES=true

# peak calling
nohup macs2 callpeak -t $path/mapped/"$sample"_mm10_sorted_rmdup.bam -c $path/mapped/"$input"_mm10_sorted_rmdup.bam \
-n "$sample"_q0.01  -B -f BAM -g mm -q 0.01 --outdir $path/called_peaks/

