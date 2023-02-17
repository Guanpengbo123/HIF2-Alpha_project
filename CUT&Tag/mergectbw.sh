#!/bin/bash
SAMPLE_NAME='HIF2 GATA3'
mkdir mergebwdata
for name in $SAMPLE_NAME
do
	samtools index mergebam/rmbam/${name}.sorted.markdup.bam
	bamCoverage -b mergebam/rmbam/${name}.sorted.markdup.bam -o mergebwdata/${name}.sorted.markdup.bw --extendReads --normalizeUsing RPGC --skipNAs --effectiveGenomeSize 2913022398 
done
