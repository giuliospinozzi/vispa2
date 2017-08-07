#!/bin/bash

CHANGEDIR="/home/andrea/test/LR40"
VECTORLTR="/opt/genome/vector/ltr/LTR.32bp.fa"
PREFIX="LR40"
R1_FASTQ="/storage/dx/backup/nas/LabDevelopment/HIV_patients/data/raw/20141003_OSR_MiSeq_LR40/LR40_S1_L001_R1_001.fastq.gz"
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;

### plot raw data stats
cd ${CHANGEDIR}
fastx_quality_stats -i <(zcat ${R1_FASTQ}) -o ${BNAME_R1}.statslog -Q 33
fastq_quality_boxplot_graph.sh -i ${BNAME_R1}.statslog -o ${BNAME_R1}.perbasequality.png -t "${BNAME_R1}"
fastx_nucleotide_distribution_graph.sh -i ${BNAME_R1}.statslog -o ${BNAME_R1}.nucleocomposition.png -t "${BNAME_R1}"


### plot LTR stats after alignment
# exmaple from my home, testing LR40
cd ${CHANGEDIR}

# align sequences to LTR genome VECTORLTR
bwa-7.5 mem -r 1 -T 15 -R "@RG\tID:${PREFIX}.RawLTR\tSM:LTR32bp\tCN:${PREFIX}.RawDataR1" -t 16 ${VECTORLTR} <(zcat ${R1_FASTQ} ) > ${PREFIX}.LTRgenome.rawR1.alm.sam

# remove reads that are:
#- not primary
#- not mapped
#- in reverse orientation
#- under quality alignment of 3 (phred)
samtools view -F 276 -q 3 -uS ${PREFIX}.LTRgenome.rawR1.alm.sam | samtools sort - ${PREFIX}.LTRgenome.rawR1.alm.sorted;
# index bam and do the MD filling
samtools index ${PREFIX}.LTRgenome.rawR1.alm.sorted.bam ;
samtools fillmd -b ${PREFIX}.LTRgenome.rawR1.alm.sorted.bam ${VECTORLTR} > ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.bam

# count per base (nucleotide) chenge and frequency (from a piled-up alignment)
# sembra funzionare ma NON benissimo....
samtools mpileup -d 1000000000 -f /opt/genome/vector/ltr/LTR.32bp.fa ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.bam > ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.pileup
/home/andrea/test/LR40/./sequenza-utils.py pileup2acgt ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.pileup









# java -jar /opt/applications/bin/picard/picard-tools-1.79/SortSam.jar I=LR40.LTRgenome.rawR1.alm.sorted.md.bam O=LR40.LTRgenome.rawR1.alm.sorted.md.coord.bam SORT_ORDER=coordinate
# samtools index LR40.LTRgenome.rawR1.alm.sorted.md.coord.bam ;

# # samstat LR40.LTRgenome.rawR1.alm.sorted.md.coord.bam

# # now marking duplicates
# MarkDuplicates INPUT=LR40.LTRgenome.rawR1.alm.sorted.md.coord.bam OUTPUT=LR40.LTRgenome.rawR1.alm.sorted.md.coord.dupflag.bam METRICS_FILE=LR40.LTRgenome.rawR1.alm.sorted.md.coord.dupflag.metrics REMOVE_DUPLICATES=false
# samtools index LR40.LTRgenome.rawR1.alm.sorted.md.coord.dupflag.bam

# # samstat LR40.LTRgenome.rawR1.alm.sorted.md.dupflag.bam

# # given the LTR bed file LTR32.bed: LTR	1	32	LTR32bp	100	+
# # coverage per base
# coverageBed -d -hist -b LTR32.bed -abam LR40.LTRgenome.rawR1.alm.sorted.md.bam


# # java jvm-args -jar /opt/applications/bin/picard/picard-tools-1.79/CollectWgsMetrics.jar I=LR40.LTRgenome.rawR1.alm.sorted.md.coord.bam O=LR40.LTRgenome.rawR1.alm.sorted.md.coord.wgsmetrics REFERENCE_SEQUENCE=/opt/genome/vector/ltr/LTR.32bp.fa COVERAGE_CAP='null'


# samtools mpileup -uf ${VECTORLTR} LR40.LTRgenome.rawR1.alm.sorted.md.coord.bam > LR40.LTRgenome.rawR1.alm.sorted.md.coord.vcf



