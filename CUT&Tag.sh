################CUT&Tag  analysis code
RAWDATA_FILE="/data/02.private/guanpb/12、13、15/ABFC20221500-12/230716_X516_0883_BHJT7TCCX2"
SAMPLE="KO-30064-C0-GATA3 KO-30064-C1-GATA3 KO-30064-C1-IgG KO-C1-IgG WT-30062-C0-GATA3 WT-30062-C1-GATA3"
for name in $SAMPLE
do
	fileR1="${RAWDATA_FILE}/${name}*.R1.fastq.gz" 
	fileR2="${RAWDATA_FILE}/${name}*.R2.fastq.gz" 
	fastp -i ${fileR1} -o cleandata/${name}_R1_clean_data.fastq.gz \
	-I ${fileR2} -O cleandata/${name}_R2_clean_data.fastq.gz \
	--thread=5 --html fastp_report/${name}_fastp.html --json fastp_report/${name}_fastp.json
done
RAWDATA_FILE="/data/02.private/guanpb/12、13、15/ABFC20221500-15/230719_X499_0815_BHJMYFCCX2"
SAMPLE="KO30036C0HIF KO30039C0HIF KOC0IgG WT30037C0HIF WT30038C0HIF WTC0IgG"
for name in $SAMPLE
do
	fileR1="${RAWDATA_FILE}/${name}*.R1.fastq.gz" 
	fileR2="${RAWDATA_FILE}/${name}*.R2.fastq.gz" 
	fastp -i ${fileR1} -o cleandata/${name}_R1_clean_data.fastq.gz \
	-I ${fileR2} -O cleandata/${name}_R2_clean_data.fastq.gz \
	--thread=5 --html fastp_report/${name}_fastp.html --json fastp_report/${name}_fastp.json
done

RAWDATA_FILE="/data/02.private/guanpb/HIF2a/230703_X517_0894_AHJM3JCCX2"
SAMPLE="KO30049C0HIF KO30049C0IgG KO30049C1HIF KO30049C1IgG KO71746C0GATA3 KO71746C0IgG KO71746C1GATA3 KO71845C0GATA3 KO71845C0IgG KO71845C1GATA3 WT30048C0HIF WT30048C0IgG WT30048C1HIF WT30048C1IgG WT71851C0GATA3 WT71851C0IgG WT71851C1GATA3 WT71859C0GATA3 WT71859C0IgG WT71859C1GATA3 WT71859C1IgG"
for name in $SAMPLE
do
	fileR1="${RAWDATA_FILE}/${name}*.R1.fastq.gz" 
	fileR2="${RAWDATA_FILE}/${name}*.R2.fastq.gz" 
	fastp -i ${fileR1} -o cleandata/${name}_R1_clean_data.fastq.gz \
	-I ${fileR2} -O cleandata/${name}_R2_clean_data.fastq.gz \
	--thread=5 --html fastp_report/${name}_fastp.html --json fastp_report/${name}_fastp.json
done


SAMPLE="KO-30064-C0-GATA3 KO-30064-C1-GATA3 KO-30064-C1-IgG KO-C1-IgG WT-30062-C0-GATA3 WT-30062-C1-GATA3 KO30036C0HIF KO30039C0HIF KOC0IgG WT30037C0HIF WT30038C0HIF WTC0IgG KO30049C0HIF KO30049C0IgG KO30049C1HIF KO30049C1IgG KO71746C0GATA3 KO71746C0IgG KO71746C1GATA3 KO71845C0GATA3 KO71845C0IgG KO71845C1GATA3 WT30048C0HIF WT30048C0IgG WT30048C1HIF WT30048C1IgG WT71851C0GATA3 WT71851C0IgG WT71851C1GATA3 WT71859C0GATA3 WT71859C0IgG WT71859C1GATA3 WT71859C1IgG"

for name in $SAMPLE
do	
	bowtie2 --very-sensitive-local --no-mixed --no-discordant -X 1000 -x /data4/cut_TAG/H4K77ac/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome \
	-1 cleandata/${name}_R1_clean_data.fastq.gz \
	-2 cleandata/${name}_R2_clean_data.fastq.gz -p 8 -S alignment/sam/${name}.sam &> alignment/sam/bowtie2_summary/${name}_bowtie2.txt
	samtools view -o alignment/bam/${name}.bam alignment/sam/${name}.sam
	rm alignment/sam/${name}.sam
	samtools sort -n alignment/bam/${name}.bam -o rmdup/sortbam/${name}.sorted.bam
	samtools fixmate -m rmdup/sortbam/${name}.sorted.bam rmdup/FixandPos/${name}.fixmate.bam
	samtools sort  rmdup/FixandPos/${name}.fixmate.bam -o rmdup/FixandPos/${name}.positionsort.bam
	samtools markdup -r rmdup/FixandPos/${name}.positionsort.bam rmdup/rmdupbam/${name}.sorted.markdup.bam -s -f rmdup/rmdupbam/rmdup_summary/${name}.rmdup.txt	
done


cd rmdup/rmdupbam/
SAMPLE="KO-30064-C0-GATA3 KO-30064-C1-GATA3 KO-30064-C1-IgG KO-C1-IgG WT-30062-C0-GATA3 WT-30062-C1-GATA3 KO30036C0HIF KO30039C0HIF KOC0IgG WT30037C0HIF WT30038C0HIF WTC0IgG KO30049C0HIF KO30049C0IgG KO30049C1HIF KO30049C1IgG KO71746C0GATA3 KO71746C0IgG KO71746C1GATA3 KO71845C0GATA3 KO71845C0IgG KO71845C1GATA3 WT30048C0HIF WT30048C0IgG WT30048C1HIF WT30048C1IgG WT71851C0GATA3 WT71851C0IgG WT71851C1GATA3 WT71859C0GATA3 WT71859C0IgG WT71859C1GATA3 WT71859C1IgG"

for name in $SAMPLE
do
	samtools index ${name}.sorted.markdup.bam
done

#####先全部都callpeaks一下

SAMPLE_NAME="KO-30064-C0-GATA3 KO-30064-C1-GATA3 KO-30064-C1-IgG KO-C1-IgG WT-30062-C0-GATA3 WT-30062-C1-GATA3 KO30036C0HIF KO30039C0HIF KOC0IgG WT30037C0HIF WT30038C0HIF WTC0IgG KO30049C0HIF KO30049C0IgG KO30049C1HIF KO30049C1IgG KO71746C0GATA3 KO71746C0IgG KO71746C1GATA3 KO71845C0GATA3 KO71845C0IgG KO71845C1GATA3 WT30048C0HIF WT30048C0IgG WT30048C1HIF WT30048C1IgG WT71851C0GATA3 WT71851C0IgG WT71851C1GATA3 WT71859C0GATA3 WT71859C0IgG WT71859C1GATA3 WT71859C1IgG"

for name in $SAMPLE_NAME
do
        bedtools genomecov -bg -ibam /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done

mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"
for name in $SAMPLE_NAME
do
        bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}.bedGraph 0.05 non stringent SEACR_005/${name}_005
        bedtools intersect -a SEACR_005/${name}_005.stringent.bed -b ${mmbt} -v >SEACR_005/${name}_005.rmbl.bed
done
@@@@@@@@@@@@@@@@@@@@@@@@@@@
#######带control callpeak
@@@@@@@@@@@@@@@@@@@@@@@@@@@

mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"
###HIF2 C1
mkdir -p SEACR_control/HIF2_C1
sample="WT71839 WT71867 WT71869 WT30048C1 KO71840 KO71865 KO71870 KO30049C1"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}HIF.bedGraph /data/02.private/guanpb/HIF2a/bgfile/${name}IgG.bedGraph non stringent SEACR_control/HIF2_C1/${name}HIF_SEACR_control
        bedtools intersect -a SEACR_control/HIF2_C1/${name}HIF_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C1/${name}HIF_control.rmbl.bed
done

### HIF2 C0 WT KO
mkdir -p SEACR_control/HIF2_C0
sample="WT30037C0 WT30038C0 WT30048C0"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}HIF.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT30048C0IgG.bedGraph non stringent SEACR_control/HIF2_C0/${name}HIF_SEACR_control
        bedtools intersect -a SEACR_control/HIF2_C0/${name}HIF_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C0/${name}HIF_control.rmbl.bed
done
sample="KO30036C0 KO30039C0 KO30049C0"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}HIF.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO30049C0IgG.bedGraph non stringent SEACR_control/HIF2_C0/${name}HIF_SEACR_control
        bedtools intersect -a SEACR_control/HIF2_C0/${name}HIF_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C0/${name}HIF_control.rmbl.bed
done

#####GATA3 C1 WT KO

mkdir -p SEACR_control/GATA3_C1
sample="WT71859C1 WT71851C1 WT-30062-C1-"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}GATA3.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT71859C1IgG.bedGraph non stringent SEACR_control/GATA3_C1/${name}GATA3_SEACR_control
        bedtools intersect -a SEACR_control/GATA3_C1/${name}GATA3_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C1/${name}GATA3_control.rmbl.bed
done

sample="KO71845C1 KO71746C1 KO-30064-C1-"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}GATA3.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO-30064-C1-IgG.bedGraph non stringent SEACR_control/GATA3_C1/${name}GATA3_SEACR_control
        bedtools intersect -a SEACR_control/GATA3_C1/${name}GATA3_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C1/${name}GATA3_control.rmbl.bed
done

#####GATA3 C0 
mkdir -p SEACR_control/GATA3_C0

sample="WT71859C0 WT71851C0 WT-30062-C0- KO71845C0 KO71746C0 KO-30064-C0-"

for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}GATA3.bedGraph /data/02.private/guanpb/HIF2a/bgfile/${name}IgG.bedGraph non stringent SEACR_control/GATA3_C0/${name}GATA3_SEACR_control
        bedtools intersect -a SEACR_control/GATA3_C0/${name}GATA3_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C0/${name}GATA3_control.rmbl.bed
done




###############################################
######筛选一个control 来进行 callpeak########
#############################################
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"
###HIF2 C1  ###KO用 KO30049C1IgG ###WT用 WT30048C1IgG
mkdir -p SEACR_control/HIF2_C1
sample="WT71839 WT71867 WT71869 WT30048C1"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}HIF.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT30048C1IgG.bedGraph norm stringent SEACR_control/HIF2_C1/${name}HIF_SEACR_control
        bedtools intersect -a SEACR_control/HIF2_C1/${name}HIF_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C1/${name}HIF_control.rmbl.bed
done
sample="KO71840 KO71865 KO71870 KO30049C1"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}HIF.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO30049C1IgG.bedGraph norm stringent SEACR_control/HIF2_C1/${name}HIF_SEACR_control
        bedtools intersect -a SEACR_control/HIF2_C1/${name}HIF_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C1/${name}HIF_control.rmbl.bed
done

### HIF2 C0 WT KO #####WT用WT30048C0IgG KO用KO30049C0IgG
mkdir -p SEACR_control/HIF2_C0
sample="WT30037C0 WT30038C0 WT30048C0"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}HIF.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT30048C0IgG.bedGraph norm stringent SEACR_control/HIF2_C0/${name}HIF_SEACR_control
        bedtools intersect -a SEACR_control/HIF2_C0/${name}HIF_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C0/${name}HIF_control.rmbl.bed
done
sample="KO30036C0 KO30039C0 KO30049C0"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}HIF.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO30049C0IgG.bedGraph norm stringent SEACR_control/HIF2_C0/${name}HIF_SEACR_control
        bedtools intersect -a SEACR_control/HIF2_C0/${name}HIF_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C0/${name}HIF_control.rmbl.bed
done

########################

#####GATA3 C1 WT KO  #####WT用 WT71859C1IgG KO用 KO-C1-IgG（GATA3C1）

mkdir -p SEACR_control/GATA3_C1
sample="WT71859C1 WT71851C1 WT-30062-C1-"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}GATA3.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT71859C1IgG.bedGraph norm stringent SEACR_control/GATA3_C1/${name}GATA3_SEACR_control
        bedtools intersect -a SEACR_control/GATA3_C1/${name}GATA3_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C1/${name}GATA3_control.rmbl.bed
done
sample="KO71845C1 KO71746C1 KO-30064-C1-"
for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}GATA3.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO-C1-IgG.bedGraph norm stringent SEACR_control/GATA3_C1/${name}GATA3_SEACR_control
        bedtools intersect -a SEACR_control/GATA3_C1/${name}GATA3_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C1/${name}GATA3_control.rmbl.bed
done

#####GATA3 C0  ###WT用 WT71859C0IgG KO用 KO71845C0IgG

mkdir -p SEACR_control/GATA3_C0
sample="WT71859C0 WT71851C0 WT-30062-C0-"

for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}GATA3.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT71859C0IgG.bedGraph non stringent SEACR_control/GATA3_C0/${name}GATA3_SEACR_control
        bedtools intersect -a SEACR_control/GATA3_C0/${name}GATA3_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C0/${name}GATA3_control.rmbl.bed
done

sample="KO71845C0 KO71746C0 KO-30064-C0-"

for name in $sample
do
		bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/${name}GATA3.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO71845C0IgG.bedGraph non stringent SEACR_control/GATA3_C0/${name}GATA3_SEACR_control
        bedtools intersect -a SEACR_control/GATA3_C0/${name}GATA3_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C0/${name}GATA3_control.rmbl.bed
done


######
samtools merge mergebam/ \
#####
mkdir -p mergebam/rmbam

samtools merge mergebam/HIF2a_C1_WT.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71839HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71867HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71869HIF.sorted.markdup.bam


samtools merge mergebam/HIF2a_C1_KO.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71840HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71865HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C1HIF.sorted.markdup.bam
###############################
samtools merge mergebam/HIF2a_C0_WT.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30037C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30038C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30048C0HIF.sorted.markdup.bam

samtools merge mergebam/HIF2a_C0_KO.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30036C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30039C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C0HIF.sorted.markdup.bam
#################################
samtools merge mergebam/GATA3_C1_WT.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71859C1GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71851C1GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT-30062-C1-GATA3.sorted.markdup.bam

samtools merge mergebam/GATA3_C1_KO.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71845C1GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71746C1GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO-30064-C1-GATA3.sorted.markdup.bam 

##############################
samtools merge mergebam/GATA3_C0_WT.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71859C0GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71851C0GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT-30062-C0-GATA3.sorted.markdup.bam

samtools merge mergebam/GATA3_C0_KO.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71845C0GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71746C0GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO-30064-C0-GATA3.sorted.markdup.bam


SAMPLE_NAME='HIF2a_C1_WT HIF2a_C1_KO HIF2a_C0_WT HIF2a_C0_KO GATA3_C1_WT GATA3_C1_KO GATA3_C0_WT GATA3_C0_KO'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done

cd mergebam/rmbam/
SAMPLE_NAME='HIF2a_C1_WT HIF2a_C1_KO HIF2a_C0_WT HIF2a_C0_KO GATA3_C1_WT GATA3_C1_KO GATA3_C0_WT GATA3_C0_KO'

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

cd /data/02.private/guanpb/HIF2a
SAMPLE_NAME='HIF2a_C1_WT HIF2a_C1_KO HIF2a_C0_WT HIF2a_C0_KO GATA3_C1_WT GATA3_C1_KO GATA3_C0_WT GATA3_C0_KO'
mkdir bgfile
for name in $SAMPLE_NAME
do
	bedtools genomecov -bg -ibam /data/02.private/guanpb/HIF2a/mergebam/rmbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done


SAMPLE_NAME='HIF2a_C1_WT HIF2a_C1_KO HIF2a_C0_WT HIF2a_C0_KO GATA3_C1_WT GATA3_C1_KO GATA3_C0_WT GATA3_C0_KO'
mkdir SEACR_005
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"

for name in $SAMPLE_NAME
do
	bash $SEACR bgfile/${name}.bedGraph 0.05 non stringent SEACR_005/${name}_005
	bedtools intersect -a SEACR_005/${name}_005.stringent.bed -b ${mmbt} -v >SEACR_005/${name}_005.rmbl.bed
done
##########################
#HIF2a C1 contral :KO用 KO30049C1IgG ###WT用 WT30048C1IgG
#HIF2a C0 contral :WT用WT30048C0IgG KO用KO30049C0IgG
#GATA3 C1 contral :WT用 WT71859C1IgG KO用 KO-C1-IgG（GATA3C1）
#GATA3 C0 contral :WT用 WT71859C0IgG KO用 KO71845C0IgG
###MASC2 进行差异peaks分析#####
mkdir MACS2_Diff
##HIF2a C1
macs2 bdgdiff --t1 ../bgfile/HIF2a_C1_WT.bedGraph --c1 ../bgfile/WT30048C1IgG.bedGraph --t2 ../bgfile/HIF2a_C1_KO.bedGraph --c2 ../bgfile/KO30049C1IgG.bedGraph --outdir HIF2a_C1_diff -o HIF2a_C1_WT.bed HIF2a_C1_KO.bed HIF2a_C1_common.bed
##HIF2a C0
macs2 bdgdiff --t1 ../bgfile/HIF2a_C0_WT.bedGraph --c1 ../bgfile/WT30048C0IgG.bedGraph --t2 ../bgfile/HIF2a_C0_KO.bedGraph --c2 ../bgfile/KO30049C0IgG.bedGraph --outdir HIF2a_C0_diff -o HIF2a_C0_WT.bed HIF2a_C0_KO.bed HIF2a_C0_common.bed
###GATA3 C1
macs2 bdgdiff --t1 ../bgfile/GATA3_C1_WT.bedGraph --c1 ../bgfile/WT71859C1IgG.bedGraph --t2 ../bgfile/GATA3_C1_KO.bedGraph --c2 ../bgfile/KO-C1-IgG.bedGraph --outdir GATA3_C1_diff -o GATA3_C1_WT.bed GATA3_C1_KO.bed GATA3_C1_common.bed
###GATA3 C0
macs2 bdgdiff --t1 ../bgfile/GATA3_C0_WT.bedGraph --c1 ../bgfile/WT71859C0IgG.bedGraph --t2 ../bgfile/GATA3_C0_KO.bedGraph --c2 ../bgfile/KO71845C0IgG.bedGraph --outdir GATA3_C0_diff -o GATA3_C0_WT.bed GATA3_C0_KO.bed GATA3_C0_common.bed
#####
################################
###common peaks 分析
####先用control callpeaks

#############################################
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"
###HIF2 C1  ###WT用 WT30048C1IgG ###KO用 KO30049C1IgG 
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C1_WT.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT30048C1IgG.bedGraph norm stringent SEACR_control/HIF2_C1/HIF2a_C1_WT_SEACR_control
bedtools intersect -a SEACR_control/HIF2_C1/HIF2a_C1_WT_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C1/HIF2a_C1_WT_SEACR_control.rmbl.bed
#KO
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C1_KO.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO30049C1IgG.bedGraph norm stringent SEACR_control/HIF2_C1/HIF2a_C1_KO_SEACR_control
bedtools intersect -a SEACR_control/HIF2_C1/HIF2a_C1_KO_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C1/HIF2a_C1_KO_SEACR_control.rmbl.bed


### HIF2 C0 WT KO #####WT用WT30048C0IgG KO用KO30049C0IgG
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C0_WT.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT30048C0IgG.bedGraph norm stringent SEACR_control/HIF2_C0/HIF2a_C0_WT_SEACR_control
bedtools intersect -a SEACR_control/HIF2_C0/HIF2a_C0_WT_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C0/HIF2a_C0_WT_SEACR_control.rmbl.bed
#KO
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C0_KO.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO30049C0IgG.bedGraph norm stringent SEACR_control/HIF2_C0/HIF2a_C0_KO_SEACR_control
bedtools intersect -a SEACR_control/HIF2_C0/HIF2a_C0_KO_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/HIF2_C0/HIF2a_C0_KO_SEACR_control.rmbl.bed


########################

#####GATA3 C1 WT KO  #####WT用 WT71859C1IgG KO用 KO-C1-IgG（GATA3C1）
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/GATA3_C1_WT.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT71859C1IgG.bedGraph norm stringent SEACR_control/GATA3_C1/GATA3_C1_WT_SEACR_control
bedtools intersect -a SEACR_control/GATA3_C1/GATA3_C1_WT_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C1/GATA3_C1_WT_SEACR_control.rmbl.bed
#KO
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/GATA3_C1_KO.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO-C1-IgG.bedGraph norm stringent SEACR_control/GATA3_C1/GATA3_C1_KO_SEACR_control
bedtools intersect -a SEACR_control/GATA3_C1/GATA3_C1_KO_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C1/GATA3_C1_KO_SEACR_control.rmbl.bed


#####GATA3 C0  ###WT用 WT71859C0IgG KO用 KO71845C0IgG
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/GATA3_C0_WT.bedGraph /data/02.private/guanpb/HIF2a/bgfile/WT71859C0IgG.bedGraph norm stringent SEACR_control/GATA3_C0/GATA3_C0_WT_SEACR_control
bedtools intersect -a SEACR_control/GATA3_C0/GATA3_C0_WT_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C0/GATA3_C0_WT_SEACR_control.rmbl.bed
#KO
bash $SEACR /data/02.private/guanpb/HIF2a/bgfile/GATA3_C0_KO.bedGraph /data/02.private/guanpb/HIF2a/bgfile/KO71845C0IgG.bedGraph norm stringent SEACR_control/GATA3_C0/GATA3_C0_KO_SEACR_control
bedtools intersect -a SEACR_control/GATA3_C0/GATA3_C0_KO_SEACR_control.stringent.bed -b ${mmbt} -v >SEACR_control/GATA3_C0/GATA3_C0_KO_SEACR_control.rmbl.bed

###################

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C1/HIF2a_C1_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C1/HIF2a_C1_KO_SEACR_control.rmbl.bed -wa -wb >c_HIF2a_C1.bed

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C1/HIF2a_C1_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/HIF2a_C1_common.bed -v > HIF2a_C1_only_WT.bed
bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C1/HIF2a_C1_KO_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/HIF2a_C1_common.bed -v > HIF2a_C1_only_KO.bed

############

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C0/HIF2a_C0_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C0/HIF2a_C0_KO_SEACR_control.rmbl.bed -wa -wb >c_HIF2a_C0.bed

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C0/HIF2a_C0_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/HIF2a_C0_common.bed -v > HIF2a_C0_only_WT.bed
bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/HIF2_C0/HIF2a_C0_KO_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/HIF2a_C0_common.bed -v > HIF2a_C0_only_KO.bed
#####

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C1/GATA3_C1_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C1/GATA3_C1_KO_SEACR_control.rmbl.bed -wa -wb >c_GATA3_C1.bed

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C1/GATA3_C1_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/GATA3_C1_common.bed -v > GATA3_C1_only_WT.bed
bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C1/GATA3_C1_KO_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/GATA3_C1_common.bed -v > GATA3_C1_only_KO.bed
######

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C0/GATA3_C0_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C0/GATA3_C0_KO_SEACR_control.rmbl.bed -wa -wb >c_GATA3_C0.bed

bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C0/GATA3_C0_WT_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/GATA3_C0_common.bed -v > GATA3_C0_only_WT.bed
bedtools intersect -a /data/02.private/guanpb/HIF2a/SEACR_control/GATA3_C0/GATA3_C0_KO_SEACR_control.rmbl.bed -b /data/02.private/guanpb/HIF2a/SEACR_control/common_peaks/GATA3_C0_common.bed -v > GATA3_C0_only_KO.bed

###############
/data/02.private/guanpb/HIF2a/MACS2_Diff/common_regulation

###HIF2a  GATA3 C1 共调控
#/data/02.private/guanpb/HIF2a/MACS2_Diff/HIF2a_C1_diff/HIF2a_C1_WT_nohead.bed
#/data/02.private/guanpb/HIF2a/MACS2_Diff/GATA3_C1_diff/GATA3_C1_WT_nohead.bed

bedtools intersect -a /data/02.private/guanpb/HIF2a/MACS2_Diff/HIF2a_C1_diff/HIF2a_C1_WT_nohead.bed -b /data/02.private/guanpb/HIF2a/MACS2_Diff/GATA3_C1_diff/GATA3_C1_WT_nohead.bed -wa -wb >c_regulation_C1.bed


###HIF2a  GATA3 C0 共调控
#/data/02.private/guanpb/HIF2a/MACS2_Diff/HIF2a_C0_diff/HIF2a_C0_WT_nohead.bed
#/data/02.private/guanpb/HIF2a/MACS2_Diff/GATA3_C0_diff/GATA3_C0_WT_nohead.bed
bedtools intersect -a /data/02.private/guanpb/HIF2a/MACS2_Diff/HIF2a_C0_diff/HIF2a_C0_WT_nohead.bed -b /data/02.private/guanpb/HIF2a/MACS2_Diff/GATA3_C0_diff/GATA3_C0_WT_nohead.bed -wa -wb >c_regulation_C0.bed


WT30048C1HIF WT71839HIF WT30048C0HIF WT30038C0HIF WT71859C1GATA3 WT71851C1GATA3 WT71859C0GATA3 WT71851C0GATA3
######分开一对一对计算
##########################
#HIF2a C1 contral :WT用 WT30048C1IgG KO用 KO30049C1IgG
#HIF2a C0 contral :WT用 WT30048C0IgG KO用KO30049C0IgG
#GATA3 C1 contral :WT用 WT71859C1IgG KO用 KO-C1-IgG（GATA3C1）
#GATA3 C0 contral :WT用 WT71859C0IgG KO用 KO71845C0IgG
###HIF2a C1
WT71839HIF KO71840HIF  #1
WT71867HIF KO71865HIF  #2
WT71869HIF KO71870HIF  #3
WT30048C1HIF KO30049C1HIF #4 
##HIF2a C0
WT30037C0HIF KO30036C0HIF  #1
WT30038C0HIF KO30039C0HIF  #2
WT30048C0HIF KO30049C0HIF  #3  
###GATA3 C1
WT71859C1GATA3 KO71845C1GATA3  #1
WT71851C1GATA3 KO71746C1GATA3  #2
WT-30062-C1-GATA3 KO-30064-C1-GATA3  #3
###GATA3 C0
WT71859C0GATA3 KO71845C0GATA3  #1
WT71851C0GATA3 KO71746C0GATA3  #2
WT-30062-C0-GATA3 KO-30064-C0-GATA3  #3

##HIF2a C1
bgfile='/data/02.private/guanpb/HIF2a/bgfile'

macs2 bdgdiff \
--t1 ${bgfile}/WT71839HIF.bedGraph \
--c1 ${bgfile}/WT30048C1IgG.bedGraph \
--t2 ${bgfile}/KO71840HIF.bedGraph \
--c2 ${bgfile}/KO30049C1IgG.bedGraph \
--outdir HIF2a_C1 -o  HIF2a_C1_WT_1.bed HIF2a_C1_KO_1.bed HIF2a_C1_common_1.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT71867HIF.bedGraph \
--c1 ${bgfile}/WT30048C1IgG.bedGraph \
--t2 ${bgfile}/KO71865HIF.bedGraph \
--c2 ${bgfile}/KO30049C1IgG.bedGraph \
--outdir HIF2a_C1 -o  HIF2a_C1_WT_2.bed HIF2a_C1_KO_2.bed HIF2a_C1_common_2.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT71869HIF.bedGraph \
--c1 ${bgfile}/WT30048C1IgG.bedGraph \
--t2 ${bgfile}/KO71870HIF.bedGraph \
--c2 ${bgfile}/KO30049C1IgG.bedGraph \
--outdir HIF2a_C1 -o  HIF2a_C1_WT_3.bed HIF2a_C1_KO_3.bed HIF2a_C1_common_3.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT30048C1HIF.bedGraph \
--c1 ${bgfile}/WT30048C1IgG.bedGraph \
--t2 ${bgfile}/KO30049C1HIF.bedGraph \
--c2 ${bgfile}/KO30049C1IgG.bedGraph \
--outdir HIF2a_C1 -o  HIF2a_C1_WT_4.bed HIF2a_C1_KO_4.bed HIF2a_C1_common_4.bed

#####HIF2a C0
macs2 bdgdiff \
--t1 ${bgfile}/WT30037C0HIF.bedGraph \
--c1 ${bgfile}/WT30048C0IgG.bedGraph \
--t2 ${bgfile}/KO30036C0HIF.bedGraph \
--c2 ${bgfile}/KO30049C0IgG.bedGraph \
--outdir HIF2a_C0 -o  HIF2a_C0_WT_1.bed HIF2a_C0_KO_1.bed HIF2a_C0_common_1.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT30038C0HIF.bedGraph \
--c1 ${bgfile}/WT30048C0IgG.bedGraph \
--t2 ${bgfile}/KO30039C0HIF.bedGraph \
--c2 ${bgfile}/KO30049C0IgG.bedGraph \
--outdir HIF2a_C0 -o  HIF2a_C0_WT_2.bed HIF2a_C0_KO_2.bed HIF2a_C0_common_2.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT30048C0HIF.bedGraph \
--c1 ${bgfile}/WT30048C0IgG.bedGraph \
--t2 ${bgfile}/KO30049C0HIF.bedGraph \
--c2 ${bgfile}/KO30049C0IgG.bedGraph \
--outdir HIF2a_C0 -o  HIF2a_C0_WT_3.bed HIF2a_C0_KO_3.bed HIF2a_C0_common_3.bed

####GATA3 C1
macs2 bdgdiff \
--t1 ${bgfile}/WT71859C1GATA3.bedGraph \
--c1 ${bgfile}/WT71859C1IgG.bedGraph \
--t2 ${bgfile}/KO71845C1GATA3.bedGraph \
--c2 ${bgfile}/KO-C1-IgG.bedGraph \
--outdir GATA3_C1 -o  GATA3_C1_WT_1.bed GATA3_C1_KO_1.bed GATA3_C1_common_1.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT71851C1GATA3.bedGraph \
--c1 ${bgfile}/WT71859C1IgG.bedGraph \
--t2 ${bgfile}/KO71746C1GATA3.bedGraph \
--c2 ${bgfile}/KO-C1-IgG.bedGraph \
--outdir GATA3_C1 -o  GATA3_C1_WT_2.bed GATA3_C1_KO_2.bed GATA3_C1_common_2.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT-30062-C1-GATA3.bedGraph \
--c1 ${bgfile}/WT71859C1IgG.bedGraph \
--t2 ${bgfile}/KO-30064-C1-GATA3.bedGraph \
--c2 ${bgfile}/KO-C1-IgG.bedGraph \
--outdir GATA3_C1 -o  GATA3_C1_WT_3.bed GATA3_C1_KO_3.bed GATA3_C1_common_3.bed

#####GATA3 C0 
macs2 bdgdiff \
--t1 ${bgfile}/WT71859C0GATA3.bedGraph \
--c1 ${bgfile}/WT71859C0IgG.bedGraph \
--t2 ${bgfile}/KO71845C0GATA3.bedGraph \
--c2 ${bgfile}/KO71845C0IgG.bedGraph \
--outdir GATA3_C0 -o  GATA3_C0_WT_1.bed GATA3_C0_KO_1.bed GATA3_C0_common_1.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT71851C0GATA3.bedGraph \
--c1 ${bgfile}/WT71859C0IgG.bedGraph \
--t2 ${bgfile}/KO71746C0GATA3.bedGraph \
--c2 ${bgfile}/KO71845C0IgG.bedGraph \
--outdir GATA3_C0 -o  GATA3_C0_WT_2.bed GATA3_C0_KO_2.bed GATA3_C0_common_2.bed

macs2 bdgdiff \
--t1 ${bgfile}/WT-30062-C0-GATA3.bedGraph \
--c1 ${bgfile}/WT71859C0IgG.bedGraph \
--t2 ${bgfile}/KO-30064-C0-GATA3.bedGraph \
--c2 ${bgfile}/KO71845C0IgG.bedGraph \
--outdir GATA3_C0 -o  GATA3_C0_WT_3.bed GATA3_C0_KO_3.bed GATA3_C0_common_3.bed

############################

SAMPLE="WT71839HIF WT71867HIF WT71869HIF WT30048C1HIF KO71840HIF KO71865HIF KO71870HIF KO30049C1HIF WT30037C0HIF WT30038C0HIF WT30048C0HIF KO30036C0HIF KO30039C0HIF KO30049C0HIF WT71859C1GATA3 WT71851C1GATA3 WT-30062-C1-GATA3 KO71845C1GATA3 KO71746C1GATA3 KO-30064-C1-GATA3 WT71859C0GATA3 WT71851C0GATA3 WT-30062-C0-GATA3 KO71845C0GATA3 KO71746C0GATA3 KO-30064-C0-GATA3"
SAMPLE_b="WT30048C1HIF KO71840HIF KO71865HIF KO71870HIF KO30049C1HIF WT30037C0HIF WT30038C0HIF WT30048C0HIF KO30036C0HIF KO30039C0HIF KO30049C0HIF WT71859C1GATA3 WT71851C1GATA3 WT-30062-C1-GATA3 KO71845C1GATA3 KO71746C1GATA3 KO-30064-C1-GATA3 WT71859C0GATA3 WT71851C0GATA3 WT-30062-C0-GATA3 KO71845C0GATA3 KO71746C0GATA3 KO-30064-C0-GATA3"

for name in $SAMPLE_b
do
	samtools view -bF 4 /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/${name}.sorted.markdup.bam > ${name}.bam 
	bedtools bamtobed -i ${name}.bam -bedpe > ${name}.bed
	awk '$1==$4 && $6-$2 < 1000 {print $0}' ${name}.bed > ${name}.clean.bed
	cut -f 1,2,6 ${name}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${name}.fragments.bed
binLen=300
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${name}.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > ${name}.fragmentsCount.bin$binLen.bed
done

######HIF2a C1 2960 885 3386 5159
bedtools intersect -a HIF2a_C1_WT_1.bed -b HIF2a_C1_WT_2.bed -wa -wb >c12_HIF2a_C1.bed
bedtools intersect -a HIF2a_C1_WT_1.bed -b HIF2a_C1_WT_3.bed -wa -wb >c13_HIF2a_C1.bed
bedtools intersect -a HIF2a_C1_WT_1.bed -b HIF2a_C1_WT_4.bed -wa -wb >c14_HIF2a_C1.bed
bedtools intersect -a HIF2a_C1_WT_2.bed -b HIF2a_C1_WT_3.bed -wa -wb >c23_HIF2a_C1.bed
bedtools intersect -a HIF2a_C1_WT_2.bed -b HIF2a_C1_WT_4.bed -wa -wb >c24_HIF2a_C1.bed
bedtools intersect -a HIF2a_C1_WT_3.bed -b HIF2a_C1_WT_4.bed -wa -wb >c34_HIF2a_C1.bed

84 c12_HIF2a_C1.bed
304 c13_HIF2a_C1.bed
267 c14_HIF2a_C1.bed
104 c23_HIF2a_C1.bed
168 c24_HIF2a_C1.bed
437 c34_HIF2a_C1.bed
###############################
###HIF2a C0 4253  6006  5528 
bedtools intersect -a HIF2a_C0_WT_1.bed -b HIF2a_C0_WT_2.bed -wa -wb >c12_HIF2a_C0.bed
bedtools intersect -a HIF2a_C0_WT_1.bed -b HIF2a_C0_WT_3.bed -wa -wb >c13_HIF2a_C0.bed
bedtools intersect -a HIF2a_C0_WT_2.bed -b HIF2a_C0_WT_3.bed -wa -wb >c23_HIF2a_C0.bed

192 c12_HIF2a_C0.bed
212 c13_HIF2a_C0.bed
457 c23_HIF2a_C0.bed

###############################
###GATA3 C1 25797 25544 3473
bedtools intersect -a GATA3_C1_WT_1.bed -b GATA3_C1_WT_2.bed -wa -wb >c12_GATA3_C1.bed
bedtools intersect -a GATA3_C1_WT_1.bed -b GATA3_C1_WT_3.bed -wa -wb >c13_GATA3_C1.bed
bedtools intersect -a GATA3_C1_WT_2.bed -b GATA3_C1_WT_3.bed -wa -wb >c23_GATA3_C1.bed

13165 c12_GATA3_C1.bed
1061 c13_GATA3_C1.bed
658 c23_GATA3_C1.bed
#############################
###GATA3 C0 27895 14036 6748 
bedtools intersect -a GATA3_C0_WT_1.bed -b GATA3_C0_WT_2.bed -wa -wb >c12_GATA3_C0.bed
bedtools intersect -a GATA3_C0_WT_1.bed -b GATA3_C0_WT_3.bed -wa -wb >c13_GATA3_C0.bed
bedtools intersect -a GATA3_C0_WT_2.bed -b GATA3_C0_WT_3.bed -wa -wb >c23_GATA3_C0.bed

13068 c12_GATA3_C0.bed
3821 c13_GATA3_C0.bed
3121 c23_GATA3_C0.bed

########


####################

sampleList<-c("WT30048C1HIF","WT71839HIF",
             "WT30048C0HIF","WT30038C0HIF",
             "WT71859C1GATA3","WT71851C1GATA3",
             "WT71859C0GATA3","WT71851C0GATA3")


sampleList<-c("WT30048C1HIF","WT71839HIF","KO30049C1HIF","KO71840HIF",
             "WT30048C0HIF","WT30038C0HIF","KO30049C0HIF","KO30039C0HIF",
             "WT71859C1GATA3","WT71851C1GATA3","KO71845C1GATA3","KO71746C1GATA3",
             "WT71859C0GATA3","WT71851C0GATA3","KO71845C0GATA3","KO71746C0GATA3")


#####

#####

samtools merge mergebam/HIF2a_C1_WT_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71839HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30048C1HIF.sorted.markdup.bam 


samtools merge mergebam/HIF2a_C1_KO_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71840HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C1HIF.sorted.markdup.bam
###############################
samtools merge mergebam/HIF2a_C0_WT_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30038C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30048C0HIF.sorted.markdup.bam

samtools merge mergebam/HIF2a_C0_KO_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30039C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C0HIF.sorted.markdup.bam
#################################
samtools merge mergebam/GATA3_C1_WT_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71859C1GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71851C1GATA3.sorted.markdup.bam 

samtools merge mergebam/GATA3_C1_KO_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71845C1GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71746C1GATA3.sorted.markdup.bam 
##############################
samtools merge mergebam/GATA3_C0_WT_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71859C0GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71851C0GATA3.sorted.markdup.bam 

samtools merge mergebam/GATA3_C0_KO_2.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71845C0GATA3.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71746C0GATA3.sorted.markdup.bam 


SAMPLE_NAME='HIF2a_C1_WT_2 HIF2a_C1_KO_2 HIF2a_C0_WT_2 HIF2a_C0_KO_2 GATA3_C1_WT_2 GATA3_C1_KO_2 GATA3_C0_WT_2 GATA3_C0_KO_2'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done

cd mergebam/rmbam/
SAMPLE_NAME='HIF2a_C1_WT_2 HIF2a_C1_KO_2 HIF2a_C0_WT_2 HIF2a_C0_KO_2 GATA3_C1_WT_2 GATA3_C1_KO_2 GATA3_C0_WT_2 GATA3_C0_KO_2'

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

cd /data/02.private/guanpb/HIF2a
SAMPLE_NAME='HIF2a_C1_WT_2 HIF2a_C1_KO_2 HIF2a_C0_WT_2 HIF2a_C0_KO_2 GATA3_C1_WT_2 GATA3_C1_KO_2 GATA3_C0_WT_2 GATA3_C0_KO_2'

#mkdir bgfile
for name in $SAMPLE_NAME
do
	bedtools genomecov -bg -ibam /data/02.private/guanpb/HIF2a/mergebam/rmbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done


SAMPLE_NAME='HIF2a_C1_WT_2 HIF2a_C1_KO_2 HIF2a_C0_WT_2 HIF2a_C0_KO_2 GATA3_C1_WT_2 GATA3_C1_KO_2 GATA3_C0_WT_2 GATA3_C0_KO_2'

#mkdir SEACR_005
mmbt='/data1/00.software/01.database/peak_blacklist/mm10-blacklist.v2.bed'
SEACR="/data1/00.software/00.common/SEACR/SEACR-1.3/SEACR_1.3.sh"

for name in $SAMPLE_NAME
do
	bash $SEACR bgfile/${name}.bedGraph 0.05 non stringent SEACR_005/${name}_005
	bedtools intersect -a SEACR_005/${name}_005.stringent.bed -b ${mmbt} -v >SEACR_005/${name}_005.rmbl.bed
done



gmtlist='c2.cp.reactome.v7.5.1.symbols c5.go.v7.5.1.symbols Total_kegg'

mkdir common_peaks_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.cls'#T_naive_versus_Th2_C0 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out common_peaks_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA -rpt_label ${gmt}
done

###############

mkdir common_peaks_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.cls'#T_naive_versus_Th2_C0 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out common_peaks_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA -rpt_label ${gmt}
done

###################

mkdir only_GATA3_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.cls'#T_naive_versus_Th2_C0 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_GATA3_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA -rpt_label ${gmt}
done

###################

mkdir only_GATA3_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.cls'#T_naive_versus_Th2_C0 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_GATA3_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA -rpt_label ${gmt}
done

###################

mkdir only_HIF2a_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.cls'#T_naive_versus_Th2_C0 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_HIF2a_C0_gene_MASC2_diff_T_naive_Th2_C0_GSEA -rpt_label ${gmt}
done

###################

mkdir only_HIF2a_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA_matrix.cls'#T_naive_versus_Th2_C0 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_HIF2a_C0_gene_TSS_MASC2_diff_T_naive_Th2_C0_GSEA -rpt_label ${gmt}
done


################

mkdir common_peaks_C1_gene_MASC2_diff_Th2_C01_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C1_gene_MASC2_diff_Th2_C01_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C1_gene_MASC2_diff_Th2_C01_GSEA_matrix.cls'#0_versus_1 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out common_peaks_C1_gene_MASC2_diff_Th2_C01_GSEA -rpt_label ${gmt}
done

################

mkdir common_peaks_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/common_peaks_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA_matrix.cls'#0_versus_1 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out common_peaks_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA -rpt_label ${gmt}
done

################

mkdir only_GATA3_C1_gene_MASC2_diff_Th2_C01_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C1_gene_MASC2_diff_Th2_C01_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C1_gene_MASC2_diff_Th2_C01_GSEA_matrix.cls'#0_versus_1 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_GATA3_C1_gene_MASC2_diff_Th2_C01_GSEA -rpt_label ${gmt}
done

################

mkdir only_GATA3_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_GATA3_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA_matrix.cls'#0_versus_1 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_GATA3_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA -rpt_label ${gmt}
done

################

mkdir only_HIF2a_C1_gene_MASC2_diff_Th2_C01_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C1_gene_MASC2_diff_Th2_C01_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C1_gene_MASC2_diff_Th2_C01_GSEA_matrix.cls'#0_versus_1 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_HIF2a_C1_gene_MASC2_diff_Th2_C01_GSEA -rpt_label ${gmt}
done

################

mkdir only_HIF2a_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA
for gmt in ${gmtlist}
do
        gsea-cli.sh GSEA -res '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA_matrix.txt' \
 -cls '/data1/02.private/guanpb/Jupyter_R_4.2.0/HIF2a_CUTTag/first_finally_analysis/GSEA/input/only_HIF2a_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA_matrix.cls'#0_versus_1 \
 -gmx /data1/02.private/guanpb/GSEA_GMT/new/${gmt}.gmt \
 -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted  -metric Diff_of_Classes -sort real -order descending -chip /data1/02.private/guanpb/GSEA_GMT/Mouse_Gene_Symbol_Remapping_Human_Orthologs_MSigDB.v7.4.chip \
 -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 5 \
-zip_report false -out only_HIF2a_C1_gene_TSS_MASC2_diff_Th2_C01_GSEA -rpt_label ${gmt}
done

#################

mkdir bwfile
SAMPLE_NAME='HIF2a_C1_WT_2 HIF2a_C1_KO_2 HIF2a_C0_WT_2 HIF2a_C0_KO_2 GATA3_C1_WT_2 GATA3_C1_KO_2 GATA3_C0_WT_2 GATA3_C0_KO_2 '
for name in $SAMPLE_NAME
do
	bamCoverage -b /data/02.private/guanpb/HIF2a/mergebam/rmbam/${name}.sorted.markdup.bam -o /data/02.private/guanpb/HIF2a/bwfile/${name}.sorted.markdup.bw --extendReads --normalizeUsing RPGC --skipNAs --effectiveGenomeSize 2652783500 
done

SAMPLE_NAME_IgG='WT30048C1IgG WT30048C0IgG WT71859C1IgG WT71859C0IgG'
for name in $SAMPLE_NAME_IgG
do
	bamCoverage -b /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/${name}.sorted.markdup.bam -o /data/02.private/guanpb/HIF2a/bwfile/${name}.sorted.markdup.bw --extendReads --normalizeUsing RPGC --skipNAs --effectiveGenomeSize 2652783500 
done

############绘制peaks分布热图
mkdir -p heatmap/mat
mkdir -p heatmap/plot

bwfile="/data/02.private/guanpb/HIF2a/bwfile"
computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 5 \
	-R '/data/02.private/guanpb/HIF2a/MACS2_Diff/new_result/HIF2a_C1/HIF2a_C1_WT_nohead.bed' \
	-S ${bwfile}/HIF2a_C1_WT_2.sorted.markdup.bw ${bwfile}/HIF2a_C1_KO_2.sorted.markdup.bw \
	-o heatmap/mat/HIF2a_C1_WT_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_HIF2a_C1_WT_TSS_3Kb.bed



plotHeatmap -m heatmap/mat/HIF2a_C1_WT_TSS_3Kb.mat.gz  \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out heatmap/plot/HIF2a_C1_WT_TSS_3Kb.svg \
--colorList "white,#009A9A" --missingDataColor 1 \
--dpi 600 --regionsLabel 'HIF2a peaks in C1_WT'

################


computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 5 \
	-R '/data/02.private/guanpb/HIF2a/MACS2_Diff/new_result/HIF2a_C0/HIF2a_C0_WT_nohead.bed' \
	-S ${bwfile}/HIF2a_C0_WT_2.sorted.markdup.bw ${bwfile}/HIF2a_C0_KO_2.sorted.markdup.bw \
	-o heatmap/mat/HIF2a_C0_WT_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_HIF2a_C0_WT_TSS_3Kb.bed



plotHeatmap -m heatmap/mat/HIF2a_C0_WT_TSS_3Kb.mat.gz  \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out heatmap/plot/HIF2a_C0_WT_TSS_3Kb.svg \
--colorList "white,#5B9090" --missingDataColor 1 \
--dpi 600 --regionsLabel 'HIF2a peaks in C0_WT'
########################

computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 5 \
	-R '/data/02.private/guanpb/HIF2a/MACS2_Diff/new_result/GATA3_C1/GATA3_C1_WT_nohead.bed' \
	-S ${bwfile}/GATA3_C1_WT_2.sorted.markdup.bw ${bwfile}/GATA3_C1_KO_2.sorted.markdup.bw \
	-o heatmap/mat/GATA3_C1_WT_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_GATA3_C1_WT_TSS_3Kb.bed




plotHeatmap -m heatmap/mat/GATA3_C1_WT_TSS_3Kb.mat.gz  \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out heatmap/plot/GATA3_C1_WT_TSS_3Kb.svg \
--colorList "white,#8F75AF" --missingDataColor 1 \
--dpi 600 --regionsLabel 'GATA3 peaks in C1_WT'

#############################
computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 5 \
	-R '/data/02.private/guanpb/HIF2a/MACS2_Diff/new_result/GATA3_C0/GATA3_C0_WT_nohead.bed' \
	-S ${bwfile}/GATA3_C0_WT_2.sorted.markdup.bw ${bwfile}/GATA3_C0_KO_2.sorted.markdup.bw \
	-o heatmap/mat/GATA3_C0_WT_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_GATA3_C0_WT_TSS_3Kb.bed




plotHeatmap -m heatmap/mat/GATA3_C0_WT_TSS_3Kb.mat.gz  \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out heatmap/plot/GATA3_C0_WT_TSS_3Kb.svg \
--colorList "white,#AE69A5" --missingDataColor 1 \
--dpi 600 --regionsLabel 'GATA3 peaks in C0_WT'



##################





samtools merge mergebam/HIF2a_C1_KO_3.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71870HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C1HIF.sorted.markdup.bam
###############################
samtools merge mergebam/HIF2a_C0_WT_3.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30037C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/.sorted.markdup.bam
WT30048C0HIF
samtools merge mergebam/HIF2a_C0_KO_3.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30036C0HIF.sorted.markdup.bam \
/data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C0HIF.sorted.markdup.bam

#################
SAMPLE_NAME='HIF2a_C1_KO_3 HIF2a_C0_WT_3 HIF2a_C0_KO_3'

for name in $SAMPLE_NAME
do
	samtools sort -n mergebam/${name}.bam -o mergebam/${name}.sorted.bam
	samtools fixmate -m mergebam/${name}.sorted.bam mergebam/${name}.fixmate.bam
	samtools sort  mergebam/${name}.fixmate.bam -o mergebam/${name}.positionsort.bam
	samtools markdup -r mergebam/${name}.positionsort.bam mergebam/rmbam/${name}.sorted.markdup.bam -s -f mergebam/rmbam/${name}.rmdup.txt
done

cd mergebam/rmbam/
SAMPLE_NAME='HIF2a_C1_KO_3 HIF2a_C0_WT_3 HIF2a_C0_KO_3'

for name in $SAMPLE_NAME
do
        samtools index ${name}.sorted.markdup.bam
done

cd /data/02.private/guanpb/HIF2a
SAMPLE_NAME='HIF2a_C1_KO_3 HIF2a_C0_WT_3 HIF2a_C0_KO_3'

#mkdir bgfile
for name in $SAMPLE_NAME
do
	bedtools genomecov -bg -ibam /data/02.private/guanpb/HIF2a/mergebam/rmbam/${name}.sorted.markdup.bam > bgfile/${name}.bedGraph
done

##############用MACS2生成的文件进行diff
macs2 callpeak -t /data/02.private/guanpb/HIF2a/mergebam/rmbam/HIF2a_C1_WT_2.sorted.markdup.bam \
               -c /data/02.private/guanpb/HIF2a/bgfile/WT30048C1IgG.bedGraph  \
               -f BAM \
               -g mm \
               -n HIF2a_C1_WT 
               --outdir /data/02.private/guanpb/HIF2a/MACS2_Diff/MACS2_analysis \
               --min-length 100 -B \
               --nomodel --extsize 200 2>macs.log


macs2 callpeak -t /data/02.private/guanpb/HIF2a/mergebam/rmbam/HIF2a_C1_WT_2.sorted.markdup.bam \
      -c /data/02.private/guanpb/HIF2a/bgfile/WT30048C1IgG.bedGraph \
      -g mm -f BAMPE -n HIF2a_C1_WT --outdir /data/02.private/guanpb/HIF2a/MACS2_Diff/MACS2_analysis --keep-dup all 2>/data/02.private/guanpb/HIF2a/MACS2_Diff/MACS2_analysis/macs2Peak_summary.txt -B

macs2 callpeak -t /data/02.private/guanpb/HIF2a/mergebam/rmbam/HIF2a_C1_WT_2.sorted.markdup.bam -c /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30048C1IgG.sorted.markdup.bam -f BAMPE -g mm -B -q 0.01 -n HIF2a_C1_WT 



macs2 callpeak -t /data/02.private/guanpb/HIF2a/mergebam/rmbam/HIF2a_C1_KO_2.sorted.markdup.bam -c /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C1IgG.sorted.markdup.bam -f BAMPE -g mm -B -q 0.01 -n HIF2a_C1_KO 

macs2 bdgdiff \
--t1 /data/02.private/guanpb/HIF2a/bgfile/WT71839HIF.bedGraph \
--c1 /data/02.private/guanpb/HIF2a/bgfile/WT30048C1IgG.bedGraph \
--t2 /data/02.private/guanpb/HIF2a/bgfile/KO71870HIF.bedGraph \
--c2 /data/02.private/guanpb/HIF2a/bgfile/KO30049C1IgG.bedGraph \
--outdir test -o  HIF2a_C1_WT_test.bed HIF2a_C1_KO_test.bed HIF2a_C1_common_test.bed




macs2 bdgdiff \
--t1 /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C1_WT_2.bedGraph \
--c1 /data/02.private/guanpb/HIF2a/bgfile/WT30048C1IgG.bedGraph \
--t2 /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C1_KO_3.bedGraph \
--c2 /data/02.private/guanpb/HIF2a/bgfile/KO30049C1IgG.bedGraph \
--outdir HIF2a_C1 -o  HIF2a_C1_WT.bed HIF2a_C1_KO.bed HIF2a_C1_common.bed

macs2 bdgdiff \
--t1 /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C0_WT_3.bedGraph \
--c1 /data/02.private/guanpb/HIF2a/bgfile/WT30048C0IgG.bedGraph \
--t2 /data/02.private/guanpb/HIF2a/bgfile/HIF2a_C0_KO_3.bedGraph \
--c2 /data/02.private/guanpb/HIF2a/bgfile/KO30049C0IgG.bedGraph \
--outdir HIF2a_C0 -o  HIF2a_C0_WT.bed HIF2a_C0_KO.bed HIF2a_C0_common.bed

SAMPLE_NAME='HIF2a_C1_KO_3 HIF2a_C0_KO_3 HIF2a_C0_WT_3'
for name in $SAMPLE_NAME
do
	bamCoverage -b /data/02.private/guanpb/HIF2a/mergebam/rmbam/${name}.sorted.markdup.bam -o /data/02.private/guanpb/HIF2a/bwfile/${name}.sorted.markdup.bw --extendReads --normalizeUsing RPGC --skipNAs --effectiveGenomeSize 2652783500 
done


###########################################
bwfile="/data/02.private/guanpb/HIF2a/bwfile"
computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 5 \
	-R '/data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C1/HIF2a_C1_WT_nohead.bed' \
	-S ${bwfile}/HIF2a_C1_WT_2.sorted.markdup.bw ${bwfile}/HIF2a_C1_KO_3.sorted.markdup.bw \
	-o heatmap/mat/HIF2a_C1_WT_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_HIF2a_C1_WT_TSS_3Kb.bed



plotHeatmap -m heatmap/mat/HIF2a_C1_WT_TSS_3Kb.mat.gz  \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out heatmap/plot/HIF2a_C1_WT_TSS_3Kb.svg \
--colorList "white,#009A9A" --missingDataColor 1 \
--dpi 600 --regionsLabel 'HIF2a peaks in C1_WT'

################


computeMatrix reference-point \
	--referencePoint TSS \
	-b 3000 -a 3000  -bs 5 \
	-R '/data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C0/HIF2a_C0_WT_nohead.bed' \
	-S ${bwfile}/HIF2a_C0_WT_3.sorted.markdup.bw ${bwfile}/HIF2a_C0_KO_3.sorted.markdup.bw \
	-o heatmap/mat/HIF2a_C0_WT_TSS_3Kb.mat.gz -p 8 --missingDataAsZero \
	--outFileSortedRegions heatmap/mat/regions_HIF2a_C0_WT_TSS_3Kb.bed



plotHeatmap -m heatmap/mat/HIF2a_C0_WT_TSS_3Kb.mat.gz  \
--whatToShow 'heatmap and colorbar' \
--heatmapWidth 8 \
-out heatmap/plot/HIF2a_C0_WT_TSS_3Kb.svg \
--colorList "white,#5B9090" --missingDataColor 1 \
--dpi 600 --regionsLabel 'HIF2a peaks in C0_WT'

#######根据peaks 计算相关性
#首先保留callpeak后区域内的bam文件
#之后计算这些bam文件的相似性
## HIF2a_WT_C1
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71839HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C1/HIF2a_C1_WT_nohead.bed > HIF2a_C1—WT_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30048C1HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C1/HIF2a_C1_WT_nohead.bed > HIF2a_C1—WT_sample2.bam
samtools index HIF2a_WT_C1_sample1.bam
samtools index HIF2a_WT_C1_sample2.bam

#####HIF2a_KO_C1
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71870HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C1/HIF2a_C1_KO_nohead.bed > HIF2a_KO_C1_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C1HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C1/HIF2a_C1_KO_nohead.bed > HIF2a_KO_C1_sample2.bam
samtools index HIF2a_KO_C1_sample1.bam
samtools index HIF2a_KO_C1_sample2.bam

##HIF2a_WT_C0
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30037C0HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C0/HIF2a_C0_WT_nohead.bed > HIF2a_WT_C0_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT30048C0HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C0/HIF2a_C0_WT_nohead.bed > HIF2a_WT_C0_sample2.bam
samtools index HIF2a_WT_C0_sample1.bam
samtools index HIF2a_WT_C0_sample2.bam

####HIF2a_KO_C0
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30036C0HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C0/HIF2a_C0_KO_nohead.bed > HIF2a_KO_C0_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO30049C0HIF.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/HIF2a_C0/HIF2a_C0_KO_nohead.bed > HIF2a_KO_C0_sample2.bam
samtools index HIF2a_KO_C0_sample1.bam
samtools index HIF2a_KO_C0_sample2.bam

##GATA3_WT_C1
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71859C1GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C1/GATA3_C1_WT_nohead.bed > GATA3_WT_C1_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71851C1GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C1/GATA3_C1_WT_nohead.bed > GATA3_WT_C1_sample2.bam
samtools index GATA3_WT_C1_sample1.bam
samtools index GATA3_WT_C1_sample2.bam

###GATA3_KO_C1
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71845C1GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C1/GATA3_C1_KO_nohead.bed > GATA3_KO_C1_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71746C1GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C1/GATA3_C1_KO_nohead.bed > GATA3_KO_C1_sample2.bam
samtools index GATA3_KO_C1_sample1.bam
samtools index GATA3_KO_C1_sample2.bam

##GATA3_WT_C0
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71859C0GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C0/GATA3_C0_WT_nohead.bed > GATA3_WT_C0_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/WT71851C0GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C0/GATA3_C0_WT_nohead.bed > GATA3_WT_C0_sample2.bam
samtools index GATA3_WT_C0_sample1.bam
samtools index GATA3_WT_C0_sample2.bam

###GATA3_KO_C0
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71845C0GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C0/GATA3_C0_KO_nohead.bed > GATA3_KO_C0_sample1.bam
bedtools intersect -a /data/02.private/guanpb/HIF2a/rmdup/rmdupbam/KO71746C0GATA3.sorted.markdup.bam -b /data/02.private/guanpb/HIF2a/MACS2_Diff/test_diff/GATA3_C0/GATA3_C0_KO_nohead.bed > GATA3_KO_C0_sample2.bam
samtools index GATA3_KO_C0_sample1.bam
samtools index GATA3_KO_C0_sample2.bam

####

multiBamSummary bins \
--bamfiles /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C1_WT/HIF2a_C1_WT_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C1_WT/HIF2a_C1_WT_sample2.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C0_WT/HIF2a_WT_C0_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C0_WT/HIF2a_WT_C0_sample2.bam \
--binSize 300 \
--numberOfProcessors 10 \
--outRawCounts HIF2a_WT_results.txt \
-o HIF2a_WT_results.npz 

plotCorrelation \
-in HIF2a_WT_results.npz \
--corMethod spearman \
--skipZeros \
--plotTitle "Sperman Correlation of Read Counts" \
--whatToPlot heatmap \
--colorMap RdYlBu \
--plotNumbers \
-o HIF2a_WT_heatmap_SpearmanCorr.pdf \
--outFileCorMatrix HIF2a_WT_SpearmanCorr_readCounts.tab

###
multiBamSummary bins \
--bamfiles /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C1_KO/HIF2a_KO_C1_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C1_KO/HIF2a_KO_C1_sample2.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C0_KO/HIF2a_KO_C0_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C0_KO/HIF2a_KO_C0_sample2.bam \
--binSize 300 \
--numberOfProcessors 10 \
--outRawCounts HIF2a_KO_results.txt \
-o HIF2a_KO_results.npz 

plotCorrelation \
-in HIF2a_KO_results.npz \
--corMethod spearman \
--skipZeros \
--plotTitle "Sperman Correlation of Read Counts" \
--whatToPlot heatmap \
--colorMap RdYlBu \
--plotNumbers \
-o HIF2a_KO_heatmap_SpearmanCorr.pdf \
--outFileCorMatrix HIF2a_KO_SpearmanCorr_readCounts.tab

####

multiBamSummary bins \
--bamfiles /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C1_KO/GATA3_KO_C1_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C1_KO/GATA3_KO_C1_sample2.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C0_KO/GATA3_KO_C0_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C0_KO/GATA3_KO_C0_sample2.bam \
--binSize 300 \
--numberOfProcessors 10 \
--outRawCounts GATA3_KO_results.txt \
-o GATA3_KO_results.npz 

plotCorrelation \
-in GATA3_KO_results.npz \
--corMethod spearman \
--skipZeros \
--plotTitle "Sperman Correlation of Read Counts" \
--whatToPlot heatmap \
--colorMap RdYlBu \
--plotNumbers \
-o GATA3_KO_heatmap_SpearmanCorr.pdf \
--outFileCorMatrix GATA3_KO_SpearmanCorr_readCounts.tab

####
multiBamSummary bins \
--bamfiles /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C1_WT/GATA3_WT_C1_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C1_WT/GATA3_WT_C1_sample2.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C0_WT/GATA3_WT_C0_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C0_WT/GATA3_WT_C0_sample2.bam \
--binSize 300 \
--numberOfProcessors 10 \
--outRawCounts GATA3_WT_results.txt \
-o GATA3_WT_results.npz 

plotCorrelation \
-in GATA3_WT_results.npz \
--corMethod spearman \
--skipZeros \
--plotTitle "Sperman Correlation of Read Counts" \
--whatToPlot heatmap \
--colorMap RdYlBu \
--plotNumbers \
-o GATA3_WT_heatmap_SpearmanCorr.pdf \
--outFileCorMatrix GATA3_WT_SpearmanCorr_readCounts.tab








重新绘制
#####

multiBamSummary bins \
--bamfiles /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C1_WT/HIF2a_C1_WT_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C1_WT/HIF2a_C1_WT_sample2.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C0_WT/HIF2a_WT_C0_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/HIF2a_C0_WT/HIF2a_WT_C0_sample2.bam \
--binSize 300 \
--numberOfProcessors 10 \
--outRawCounts HIF2a_WT_results.txt \
-o HIF2a_WT_results.npz 

plotCorrelation \
-in HIF2a_WT_results.npz \
--corMethod spearman \
--skipZeros \
--whatToPlot heatmap \
--colorMap seismic \
--plotNumbers \
-o HIF2a_WT_heatmap_SpearmanCorr_seismic.svg \
--outFileCorMatrix HIF2a_WT_SpearmanCorr_readCounts.tab \
--zMin -1 --zMax 1 --labels HIF2a_C1_rep1 HIF2a_C1_rep2 HIF2a_C0_rep1 HIF2a_C0_rep2

plotCorrelation \
-in HIF2a_WT_results.npz \
--corMethod spearman \
--skipZeros \
--whatToPlot heatmap \
--colorMap bwr \
--plotNumbers \
-o HIF2a_WT_heatmap_SpearmanCorr_bwr.svg \
--outFileCorMatrix HIF2a_WT_SpearmanCorr_readCounts.tab \
--zMin -1 --zMax 1 --labels HIF2a_C1_rep1 HIF2a_C1_rep2 HIF2a_C0_rep1 HIF2a_C0_rep2

######
multiBamSummary bins \
--bamfiles /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C1_WT/GATA3_WT_C1_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C1_WT/GATA3_WT_C1_sample2.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C0_WT/GATA3_WT_C0_sample1.bam /data/02.private/guanpb/HIF2a/peaks_cor_analysis/GATA3_C0_WT/GATA3_WT_C0_sample2.bam \
--binSize 300 \
--numberOfProcessors 10 \
--outRawCounts GATA3_WT_results.txt \
-o GATA3_WT_results.npz 

plotCorrelation \
-in GATA3_WT_results.npz \
--corMethod spearman \
--skipZeros \
--whatToPlot heatmap \
--colorMap seismic \
--plotNumbers \
-o GATA3_WT_heatmap_SpearmanCorr_seismic.svg \
--outFileCorMatrix GATA3_WT_SpearmanCorr_readCounts.tab \
--zMin -1 --zMax 1 --labels GATA3_C1_rep1 GATA3_C1_rep2 GATA3_C0_rep1 GATA3_C0_rep2

plotCorrelation \
-in GATA3_WT_results.npz \
--corMethod spearman \
--skipZeros \
--whatToPlot heatmap \
--colorMap bwr \
--plotNumbers \
-o GATA3_WT_heatmap_SpearmanCorr_bwr.svg \
--outFileCorMatrix GATA3_WT_SpearmanCorr_readCounts.tab \
--zMin -1 --zMax 1 --labels GATA3_C1_rep1 GATA3_C1_rep2 GATA3_C0_rep1 GATA3_C0_rep2





