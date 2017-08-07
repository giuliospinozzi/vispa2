#!/bin/bash -x

#ONLY BWA-MEM HIV

POOL="LR22"
SOURCEDIR="/opt/NGS/data/HIV/454/LR22"


GENOME="/opt/genome/human/hg19/index/bwa_7/hg19.fa"
OUTDIR="/storage/d3/ngs/pipelinetmpdir/HIV_purealigment_${POOL}"
mkdir ${OUTDIR}
cd ${OUTDIR}

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Align"
for k in $( ls ${SOURCEDIR}/*.fa ); do
TAG=`basename ${k} | cut -d'.' -f1 `
bwa-7.5 mem -r 1 -M -T 15 -R "@RG\tID:${POOL}.${TAG}\tSM:${POOL}.${TAG}\tCN:${POOL}" -t 10 ${GENOME} ${k} > ${TAG}.sam
samtools view -F 260 -q 5 -uS ${TAG}.sam > ${TAG}.bam;
samtools sort ${TAG}.bam ${TAG}.sorted;
samtools index ${TAG}.sorted.bam ;
samtools fillmd -b ${TAG}.sorted.bam ${GENOME} > ${TAG}.sorted.md.bam
samtools index ${TAG}.sorted.bam
gzip ${TAG}.sam
rm ${TAG}.bam
done

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Merge BAM with PICARD"
INPUT="INPUT";
INSTR="";
for k in $( ls *sorted.bam ); do
	INSTR+=" ${INPUT}=${k}";
done
MergeSamFiles ${INSTR} OUTPUT="${POOL}.sorted.wholepool.bam" SORT_ORDER=coordinate ;
# index grouped bam
samtools index ${POOL}.sorted.wholepool.bam
