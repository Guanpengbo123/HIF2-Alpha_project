#!/bin/bash

mkdir -p mergebam/rmbam

samtools merge mergebam/HIF2.bam \
rmdup/rmdupbam/B-2.sorted.markdup.bam \
rmdup/rmdupbam/B-3.sorted.markdup.bam 


samtools merge mergebam/GATA3.bam \
rmdup/rmdupbam/Treat-1.sorted.markdup.bam \
rmdup/rmdupbam/Treat-2.sorted.markdup.bam 

SAMPLE_NAME='HIF2 GATA3'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done


