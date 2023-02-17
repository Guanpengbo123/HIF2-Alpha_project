#！/bin/bash
blacklist="/data1/00.software/01.database/peak_blacklist/hg38-blacklist.v2.2.bed"
SAMPLE="B-2 B-3 Treat-1 Treat-2 control-3"

mkdir /data4/cut_TAG/HIF2_Gata3/rmblbed

for name in $SAMPLE
do
	bedtools intersect -a /data4/cut_TAG/HIF2_Gata3/bedgraph/callpeak/Top_005/${name}_005.stringent.bed -b ${blacklist} -v > /data4/cut_TAG/HIF2_Gata3/rmblbed/${name}_005.stringent.rmbl.bed
done

