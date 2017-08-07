#!/bin/bash
### quantify plasmids from RAW files (gemini server)
### Usage example:
### ~/test/shearsites$ ./quantify_plasmid_bypool.sh test100r.r1.fastq.gz test100r.r2.fastq.gz /home/andrea/test/shearsites PLASMID.log
### /home/andrea/test/shearsites/./quantify_plasmid_bypool.sh /storage/dx/backup/nas/LabDevelopment/Monica/data/20141015_OSR_MiSeq_MV3/MV3_S1_L001_R1_001.fastq.gz /storage/dx/backup/nas/LabDevelopment/Monica/data/20141015_OSR_MiSeq_MV3/MV3_S1_L001_R2_001.fastq.gz /home/andrea/test/shearsites PLASMID.log

R1_FASTQ="${1}"
R2_FASTQ="${2}"
MAXTHREADS=10
DDIR="${3}" ; #"/home/andrea/test/shearsites"
LOGFILE="${4}"

##### ================ PLASMIDS quantification POOL based - start ======================== #####
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
for PLASMID in kana amp; do 
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference ${PLASMID} genome"
	## --- paired-end reads ---
	# low seed value due to small genome
	bwa-stable mem -k 14 -r 1 -T 15 -c 1 -R "@RG\tID:plasmid\tSM:${PLASMID}\tCN:Andrea" -t ${MAXTHREADS} /opt/genome/vector/plasmid/plasmids.${PLASMID}.fa <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) > ${DDIR}/${BNAME_R1}.${PLASMID}.sam
	# quality filter of seed-1 value
	samtools view -F 2308 -q 20 -uS ${DDIR}/${BNAME_R1}.${PLASMID}.sam | samtools sort - ${DDIR}/${BNAME_R1}.${PLASMID}.sorted
	rm ${DDIR}/${BNAME_R1}.${PLASMID}.sam
	# fastqc -o ${OUTDIR_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${BASEDIR}/${DISEASE}_${PATIENT}_${POOL}.LTR${b}.LC${k}.noLTRLC.${PLASMID}.merge.bam
	samtools view ${DDIR}/${BNAME_R1}.${PLASMID}.sorted.bam | cut -f1 > ${DDIR}/${BNAME_R1}.${PLASMID}.list
done
cat ${DDIR}/${BNAME_R1}.*.list | sort | uniq | gzip > ${DDIR}/${BNAME_R1}.plasmids.list.gz
rm ${DDIR}/${BNAME_R1}.*.sorted.bam
rm ${DDIR}/${BNAME_R1}.*.list

# get the number of remaining sequences
RAW_READS=$((`zcat ${R1_FASTQ} | wc -l | cut -d' ' -f1 `/4)) ;
PLASMID_MAPPING_READS=`zcat ${DDIR}/${BNAME_R1}.plasmids.list.gz | wc -l  | cut -d' ' -f1 `
##### ================ PLASMIDS quantification POOL based - end ======================== #####
echo -e "PoolFile\tPoolName\tRawReads\tPlasmidReads" >> ${LOGFILE}
echo -e "${R1_FASTQ}\t${BNAME_R1}\t${RAW_READS}\t${PLASMID_MAPPING_READS}" >> ${LOGFILE}
