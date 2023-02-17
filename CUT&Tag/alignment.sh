#!/bin/bash

mkdir -p alignment/sam/bowtie2_summary
mkdir -p alignment/bam
mkdir -p alignment/bed
mkdir -p rmdup/sortbam
mkdir -p rmdup/FixandPos
mkdir -p rmdup/rmdupbam
mkdir -p rmdup/rmdupbam/rmdup_summary

ref='/data4/cut_TAG/ref/Hg_index'

GATA3='control-3 Treat-1 Treat-2'
HIF2="B-2 B-3"

for name in $GATA3
do
	bowtie2 --very-sensitive-local --no-mixed --no-discordant -X 1000 -x ${ref} \
	-1 /data4/cut_TAG/HIF2_Gata3/cleandata/GATA3/${name}_R1_clean_data.fastq.gz \
 	-2 /data4/cut_TAG/HIF2_Gata3/cleandata/GATA3/${name}_R2_clean_data.fastq.gz -p 8 -S alignment/sam/${name}.sam &> alignment/sam/bowtie2_summary/${name}_bowtie2.txt
	samtools view -o alignment/bam/${name}.bam alignment/sam/${name}.sam	
done

for name in $HIF2
do
	bowtie2 --very-sensitive-local --no-mixed --no-discordant -X 1000 -x ${ref} \
	-1 /data4/cut_TAG/HIF2_Gata3/cleandata/HIF2/${name}_R1_clean_data.fastq.gz \
 	-2 /data4/cut_TAG/HIF2_Gata3/cleandata/HIF2/${name}_R2_clean_data.fastq.gz -p 8 -S alignment/sam/${name}.sam &> alignment/sam/bowtie2_summary/${name}_bowtie2.txt
	samtools view -o alignment/bam/${name}.bam alignment/sam/${name}.sam
	
done


for name in $GATA3
do
	samtools sort -n alignment/bam/${name}.bam -o rmdup/sortbam/${name}.sorted.bam
	samtools fixmate -m rmdup/sortbam/${name}.sorted.bam rmdup/FixandPos/${name}.fixmate.bam
	samtools sort  rmdup/FixandPos/${name}.fixmate.bam -o rmdup/FixandPos/${name}.positionsort.bam
	samtools markdup -r rmdup/FixandPos/${name}.positionsort.bam rmdup/rmdupbam/${name}.sorted.markdup.bam -s -f rmdup/rmdupbam/rmdup_summary/${name}.rmdup.txt
done

for name in $HIF2
do
	samtools sort -n alignment/bam/${name}.bam -o rmdup/sortbam/${name}.sorted.bam
	samtools fixmate -m rmdup/sortbam/${name}.sorted.bam rmdup/FixandPos/${name}.fixmate.bam
	samtools sort  rmdup/FixandPos/${name}.fixmate.bam -o rmdup/FixandPos/${name}.positionsort.bam
	samtools markdup -r rmdup/FixandPos/${name}.positionsort.bam rmdup/rmdupbam/${name}.sorted.markdup.bam -s -f rmdup/rmdupbam/rmdup_summary/${name}.rmdup.txt
done
