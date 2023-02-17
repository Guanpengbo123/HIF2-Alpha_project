#!/bin/bash


mkdir -p SEACR_callpeak/Top_005_new
SAMPLE_NAME='HIF2 GATA3'


SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"
for name in $SAMPLE_NAME
do
        bash $SEACR bgfile/${name}.bedGraph 0.05 non stringent SEACR_callpeak/Top_005_new/${name}_005
done


blacklist="/data1/00.software/01.database/peak_blacklist/hg38-blacklist.v2.2.bed"

for name in $SAMPLE_NAME
do
	bedtools intersect -a /data4/cut_TAG/HIF2_Gata3/SEACR_callpeak/Top_005_new/${name}_005.stringent.bed -b ${blacklist} -v > /data4/cut_TAG/HIF2_Gata3/SEACR_callpeak/Top_005_new/${name}_005.stringent.rmbl.bed
done

