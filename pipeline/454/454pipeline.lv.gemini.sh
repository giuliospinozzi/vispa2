#!/bin/bash -x
source /etc/environment
source /etc/profile

PIPEVERSION="3.0.5"
STARTTIME=`date +'%y-%m-%d %H:%M:%S'`

echo "
=====================================================================================
		NGS pipeline for 454 Roche -> GAMMA AND LENTI RV
=====================================================================================
Date:	03rd December 2012 v0.1
		07th May 2013 v0.2
		27th August 2013 v0.2.1 -> reduced to accept 1 barcode per time and FastA
		04th April 2014 v1 -> BWA-MEM (ref to Illumina 3.0.4)
Author:	Andrea Calabria (andrea.calabria@hsr.it)
Results: available on Alien web server at: http://172.25.39.2/Workbench/...
=====================================================================================
	
	REQUIRED VARS WITH ORDER POSITION -> REMEMBER NO SPACES!!!!!!!!!
	1.	disease id (MLD)
	2.	patient id (MLD02)
	3.	local server working path [/opt/NGS/results] -> it will be added DISEASE and PATIENT ID [/opt/NGS/results/MLD/MLD02/]
	4.	FASTA file (ex SFF file) -> from the file name I will derive the TAG/barcode ID/sequence
	5.	pool id
	6.	barcode list
	7.	reference genome (NB: indexed by BWA)
	8.	tmp dir
	9.	association file
	10.	db host id (write local as default)
	11.	target db schema for redundant iss
	12.	target db table for redundant iss
	13.	LTR file (gamma or lenti?)
	14.	LC file
	15.	Genome ID for CIGAR-MD program. NB: this option collects only whole genome so far, not single regions (even the program allows it!). Example: hg19, mm9, mfa5, ce, lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada. [CIGARGENOMEID]
	16.	Genome ID for CIGAR-MD program for VECTOR genome. Example: lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada [VECTORCIGARGENOMEID]
	17.	SUBOPTIMALTHRESHOLD to compare first vs second hit (AS vs XS tags), this is a percentage value (int).

"

echo "
Preconditions
- input file format: FASTA
- input folder must contain all FASTA files named with the barcode: givean a barcode AATTCC the fasta file will be AACCTT.fa

"

DISEASE="${1}";
PATIENT="${2}";
POOL="${5}";
DBHOSTID="${10}";
DBSCHEMA="${11}";
DBTABLE="${12}";
SERVERWORKINGPATH="${3}"; # without final /

FASTA="${4}";
#ASSOCIATIONFILE="/home/andrea/Dropbox/TIGET/Workspace/PROJECTS/ADA/FileLegge/1.association_files/CB05/asso.cb05.gb.csv";
ASSOCIATIONFILE="${9}";
GENOME="${7}"; # the right one
TMPDIR="${8}";
LTR="${13}";
LC="${14}";
BASEPATH="${SERVERWORKINGPATH}/${DISEASE}/${PATIENT}/${POOL}"; #/opt/NGS/results/ADA/poolCB05";
OUTDIR="mux";
BASENAME="${DISEASE}.${PATIENT}.${POOL}";
BARCODELIST="${6}";
##BCF=`cat "$BARCODELIST" | cut -f1`; 
##BARCODEARRAY=( ${BCF} ); # array of barcodes only
TAG=`basename ${FASTA} | cut -d'.' -f 1`; # NB: now the TAG id (that is mainly the tag sequence) is derived from the FASTA name!!!!!
CIGARGENOMEID="${15}" ;
VECTORCIGARGENOMEID="${16}";
SUBOPTIMALTHRESHOLD="${17}"
MAXTHREADS="10"
VECTORGENOME="/opt/genome/vector/lv/bwa_7/lv.backbone.fa" # gemini set up

# create tmp dir
mkdir ${SERVERWORKINGPATH}
mkdir ${SERVERWORKINGPATH}/${DISEASE}
mkdir ${SERVERWORKINGPATH}/${DISEASE}/${PATIENT}
mkdir ${SERVERWORKINGPATH}/${DISEASE}/${PATIENT}/${POOL}
#mkdir ${SERVERWORKINGPATH}/${DISEASE}/${PATIENT}/${POOL}/bam
mkdir ${BASEPATH}/${OUTDIR}
mkdir ${TMPDIR}
mkdir ${TMPDIR}/reads
mkdir ${TMPDIR}/bam
mkdir ${TMPDIR}/bed
mkdir ${TMPDIR}/sam
mkdir ${TMPDIR}/lv
mkdir ${TMPDIR}/${OUTDIR}
# checking vars
RUN_ID=`date +"%Y%m%d%H%M%S"`
RUN_NAME="${DISEASE}|${PATIENT}|${POOL}|${TAG}"
PHIX_MAPPING="0"
BARCODE_MUX="0" # to add in case of barcoding demux


echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Preprocessing input variables (delimiters:<>)"
## print input variables (check for log utils)
INPUTVARNUM=0
for INPUTVAR in "$@"; do
	let INPUTVARNUM++; 
	printf -v INPUTNUM '%02d' $INPUTVARNUM;
    echo "  => Input Variable: Order number = <${INPUTNUM}> ; Var Content = <${INPUTVAR}>";
done


### ------------------ DB SCHEMA AND STATS SETUP ------------- ###
mysql -u andrea --password=andrea -e "create database if not exists ${DBSCHEMA} ; " ;

mysql -u andrea --password=andrea -e "
	CREATE TABLE IF NOT EXISTS stats_summary_454 (
		RUN_ID varchar(255) default null,
		RUN_NAME varchar(1000) default null,
		DISEASE varchar(255) default null,
		PATIENT varchar(255) default null,
		POOL varchar(255) default null,
		TAG varchar(255) default null,
		LTR_ID varchar(255) default null,
		LC_ID varchar(255) default null,
		PHIX_MAPPING int(30) DEFAULT NULL,
		BARCODE_MUX int(30) DEFAULT NULL,
		LTR_IDENTIFIED int(30) DEFAULT NULL,
		LTR_SEQAFTERTRIMMING int(30) DEFAULT NULL,
		LC_SEQAFTERTRIMMING int(30) DEFAULT NULL,
		BWA_RESULTS varchar(1000) default null,
		BWA_INPUT int(30) DEFAULT NULL,
		BWA_MAPPED int(30) DEFAULT NULL,
		BWA_MAPPED_PP int(30) DEFAULT NULL,
		BWA_MAPPED_ST int(30) DEFAULT NULL,
		BWA_MAPPED_OVERALL int(30) DEFAULT NULL,
		FILTER_QUAL_RESULTS varchar(1000) default null,
		FILTER_QUAL_PP int(30) DEFAULT NULL,
		FILTER_QUAL_ST int(30) DEFAULT NULL,
		FILTER_QUAL_OVERALL int(30) DEFAULT NULL,
		FILTER_ISSR_RESULTS varchar(1000) default null,
		FILTER_ISSR_PP int(30) DEFAULT NULL,
		FILTER_ISSR_ST int(30) DEFAULT NULL,
		FILTER_ISSR_OVERALL int(30) DEFAULT NULL,
		CIGARREADS_REMOVED int(30) DEFAULT NULL,
		ISS_IDENTIFIED int(30) DEFAULT NULL,
		LV_MAPPING_RESULTS varchar(1000) default null,
		LV_MAPPING_PP int(30) DEFAULT NULL,
		LV_MAPPING_ST int(30) DEFAULT NULL,
		LV_MAPPING_OVERALL int(30) DEFAULT NULL,
		LV_MAPPING_MAPPED int(30) DEFAULT NULL
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;
	" ${DBSCHEMA} ;


### ------------------ CONVERTING FASTA TO FASTQ ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Converting data FastA to FastQ"
# step 1 convert fasta into fastq, NB to fastq encoding!
fasta_to_fastq -qv 9 ${FASTA} > ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq

### ------------------ MARKING LTR READS ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Marking LTR"
# mark LTR reads and create list of only LTR reads (SAM accepted)
## ORIGINAL: # flexbar --source <(zcat ${TMPDIR}/${BASENAME}.fastq.gz) --target ${TMPDIR}/flexbar -f fastq-i1.8 -a ${LTR} --threads 8 -ae ANY -at 1 -ao 58 -aa -m 20 -q 10 --removal-tag ## this is for LENTI
flexbar --source ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq --target ${TMPDIR}/reads/${BASENAME}.${TAG}.flexbar -f fastq-i1.8 -a ${LTR} --threads 8 -ae ANY -at 3 -ao 25 -aa -m 20 -q 10 --removal-tag ## this is for LENTI
###flexbar --source <(zcat ${TMPDIR}/${BASENAME}.fastq.gz) --target ${TMPDIR}/flexbar -f fastq-i1.8 -a ${LTR} --threads 8 -ae ANY -at 1 -ao 53 -aa -m 20 -q 10 --removal-tag ## this is for GAMMA retro
grep "Flexbar" ${TMPDIR}/reads/${BASENAME}.${TAG}.flexbar.fastq > ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.only
cut -d' ' -f1 ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.only | sed -e 's/@//g' > ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.sam.list
# change input data -> leave only recognized LTR reads
mv ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq ${TMPDIR}/reads/${BASENAME}.${TAG}.allreads.fastq
cat ${TMPDIR}/reads/${BASENAME}.${TAG}.allreads.fastq | fqextract_pureheader ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.sam.list | gzip -f > ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq.gz
rm ${TMPDIR}/reads/${BASENAME}.${TAG}.allreads.fastq

# #### =========== COMMENTED - START ===== ANDREA AUGUST 2013 ========= #####
# # mux barcodes
# ## this is our FINAL version of barcode MUX ...but here is not working fine with 0 mismatches
# #fastq-multx -m 0 -B ${BARCODELIST} ${TMPDIR}/${BASENAME}.fastq.gz -o ${TMPDIR}/${OUTDIR}/${BASENAME}.%.fastq.gz 
# fastq-multx -B ${BARCODELIST} ${TMPDIR}/${BASENAME}.fastq.gz -o ${TMPDIR}/${OUTDIR}/${BASENAME}.%.fastq.gz 
# #### =========== COMMENTED - END ===== ANDREA AUGUST 2013 ========= #####

### ------------------ start TRIMMING LTR : GAMMA and LENTI ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Trimming LTR"
#fastq-mcf ${LTR} mld02.454.pool6.${TAG}.fastq.gz -o mld02.454.pool6.${TAG}.noLTR.fastq.gz -m 61 -p 2 -l 20 -q 10 -P 33; 
#flexbar --source <(zcat mld02.454.pool6.${TAG}.fastq.gz) --target mld02.454.pool6.${TAG}.noLTR -f fastq-i1.8 -as "TGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA" --threads 8 -ae LEFT -at 2 -ao 61 -aa -m 20 -q 10
### con questi parametri ottimizzi il riconoscimento LTR. attenzione ad un comportamento indesiderato di flexbar: i gap e mismatches fanno shiftare la sequenza di output!! per questo motivo ho impostato at a 0.4 perchè significa 0.4 errori ogni 10 basi e quindi possiamo arrivare a tollerare LTR di lunghezza massima 74 bp con al piu 2 errori
##### --->>>> 63 bp LTR:: flexbar --source <(zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq.gz ) --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR -f fastq-i1.8 -a ${LTR} --threads 8 -ae LEFT -at 0.5 -ai -10 -ao 28 -m 20 -q 10 ## this is for LENTI
flexbar --source <(zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq.gz ) --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR -f fastq-i1.8 -a ${LTR} --threads 8 -ae LEFT -at 0.72 -ai -10 -ao 28 -m 20 -q 10 ## this is for LENTI

##now:: use eautils
#fastq-mcf ${LTR} ${TMPDIR}/${TAG}.fastq.gz -o ${TMPDIR}/${OUTDIR}/${BASENAME}.${TAG}.noLTR.fastq.gz -m 26 -p 80 -l 20 -q 10 -P 33 -R -k 0 -x 0 -S

## se vuoi 1 solo errore -> cosi eviti derive eccessive anche se ne butti di piu
#flexbar --source <(zcat mld02.454.pool6.${TAG}.fastq.gz) --target mld02.454.pool6.${TAG}.noLTR -f fastq-i1.8 -a ${LTR} --threads 8 -ae LEFT -at 0.3 -ai -10 -ao 50 -aa -m 20 -q 10

### per gamma molto stringente:
##flexbar --source <(zcat ${TMPDIR}/${OUTDIR}/${BASENAME}.${TAG}.fastq.gz) --target ${TMPDIR}/${OUTDIR}/${BASENAME}.${TAG}.noLTR -f fastq-i1.8 -a ${LTR} --threads 8 -ae LEFT -at 0.2 -ai -10 -ao 55 -aa -m 20 -q 10 ## this is for GAMMA

### ------------------ end TRIMMING LTR ------------- ###

# flexbar --> you must then gzip files
gzip -f ${TMPDIR}/reads/${BASENAME}.*.fastq

### ------------------ TRIMMING LC ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Trimming LC"
#fastq-mcf ${LC} mld02.454.pool6.${TAG}.noLTR.fastq.gz -o mld02.454.pool6.${TAG}.noLTRLC.fastq.gz -m 12 -p 20 -l 30 -q 10 -P 33
## la precedente istruzione fallisce nel caso di più LC concatenate, esempio la read HHAUOBH02JFYTO
flexbar --source <(zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR.fastq.gz ) --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC -f fastq-i1.8 -a ${LC} --threads 8 -ae RIGHT -at 4 -ao 8 -aa -m 10 -q 10 ;

# flexbar --> you must then gzip files
gzip -f ${TMPDIR}/reads/${BASENAME}.*.noLTRLC.fastq


    
##### ================ LV removal TAG based - start ======================== #####
OUTDIR_VECTOR_POOL="${TMPDIR}/lv";
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference VECTOR genome"
echo "bwa-7.5 mem -r 1 -T 15 -R \"@RG\tID:${TAG}\tSM:Vector\tCN:${DISEASE}.${PATIENT}.${POOL}.Andrea\" -t ${MAXTHREADS} ${VECTORGENOME} <(zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq.gz) > ${TMPDIR}/${DISEASE}_${PATIENT}_${POOL}.${TAG}.noLTRLC.sam"
bwa-7.5 mem -r 1 -T 15 -R "@RG\tID:${TAG}\tSM:Vector\tCN:${DISEASE}.${PATIENT}.${POOL}.Andrea" -t ${MAXTHREADS} ${VECTORGENOME} <(zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq.gz) > ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sam

# create BAM and sort them
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter LV from here using only valid reads: mapped and primary)"
samtools view -F 260 -q 30 -uS ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sam > ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.bam ;
samtools sort ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.bam ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted;
samtools index ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bam;
samtools fillmd -b ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bam ${VECTORGENOME} > ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam
samtools index ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam ;

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filter by CIGAR and MD tag (only for starting position and CIGAR value) and then select only LV reads from BAM."
filter_by_cigar_bam --bam ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam -o ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list --minStartingMatches 15 -p ${MAXTHREADS} -g ${VECTORCIGARGENOMEID} --minStartingBasesNoIndels 20 --endClipThreshold 5 --singleEnd

CHIMERANROWS=`wc -l ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list | cut -d' ' -f1`
if [ $CHIMERANROWS -gt 0 ] 
	then
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- found reads not LV to remove, now removing them all from BAM file"
	FilterSamReads INPUT=${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam FILTER=excludeReadList RLF=${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list SO=coordinate O=${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam
else
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- reads to remove not found, LV BAM is clean, go ahead (are you sure to go ahead...?! all reads are LV!!!)"
	cp ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam
fi
samtools index ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Statistics for LV BAM (both fastqc and samstat)"
samstat ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam 

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import data into DB from BAM (new table different from the standard one)" 
import_iss --bam ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam --singleEnd --locusRead alnstart -p ${MAXTHREADS} -g ${VECTORCIGARGENOMEID} --dbhost 172.25.39.57 --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_lv --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG};

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no LV reads from raw data"
zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq.gz | fqreverseextract_pureheader ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list | gzip -f > ${TMPDIR}/reads/${BASENAME}.no12.noLTRLC.noVector.fastq.gz 

# clean lv folder
rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sam;
rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.bam;
#rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam
rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bam ;
gzip -f ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list

##### ================ LV removal TAG based - end ======================== #####


### ------------------ ALIGNMENT ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference genome"
#---
bwa-7.5 mem -r 1 -M -T 15 -R "@RG\tID:${DISEASE}.${PATIENT}.${POOL}.${TAG}\tSM:${TAG}\tCN:Andrea.${DISEASE}.${PATIENT}.${POOL}" -t ${MAXTHREADS} ${GENOME} <(zcat ${TMPDIR}/reads/${BASENAME}.no12.noLTRLC.noVector.fastq.gz ) > ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam

# create BAM and sort them
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter from here the dataset using only valid reads: mapped and primary)"
echo "samtools view -F 260 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam > ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam;"
samtools view -F 260 -q 5 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam > ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam;
samtools sort ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted;
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.bam ;
chmod -R 755 ${TMPDIR}/bam/*
#rm ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam
rm ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam

### ------------------ RECALIBRATION ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Recreate/Recalibrate/Fill the MD tag"
samtools fillmd -b ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.bam ${GENOME} > ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam 
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam
rm ${TMPDIR}/bam/${BASENAME}.${TAG}.sorted.bam

### ------------------ FILTERING ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filter by CIGAR and MD tag"
filter_by_cigar_bam --bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam -o ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list --minStartingMatches 3 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --minStartingBasesNoIndels 2 --compareSubOptimal --suboptimalThreshold ${SUBOPTIMALTHRESHOLD} --SAalnAction ignore --XSlikeTag XS --ASlikeTag AS --endClipThreshold 18 --singleEnd 

CHIMERANROWS=`wc -l ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list | cut -d' ' -f1`
if [ $CHIMERANROWS -gt 0 ] 
	then
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- found reads to remove, now removing them all from BAM file"
	FilterSamReads INPUT=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam FILTER=excludeReadList RLF=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list SO=coordinate O=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam
else
	echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] --- reads to remove not found, go ahead"
	cp ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam
fi
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam
	
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filter BAM with only LTR reads -> out is a BAM file with only valid integration sites (redundant)"
FilterSamReads INPUT=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam FILTER=includeReadList READ_LIST_FILE=${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.sam.list OUTPUT=${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam
samtools index {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data (Bamtools)"
bamtools filter -in {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam -mapQuality ">=12" -out {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.hq.bam
#samtools index {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.hq.bam

# rename files and index them
mv {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.lq.bam
mv {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.hq.bam {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam
samtools index {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.lq.bam
samtools index {TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam

### ------------------ DATA FORMATTING ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Convert BAM to BED (returns score as CIGAR)"
bamToBed -cigar -i ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam > ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed

### ------------------ IMPORT DATA INTO DB ------------- ###
NOWIS=`date +'%y-%m-%d %H:%M:%S'`
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import BED data into DB (my script - running only for MY USER! so far)" 
isa_importrediss_frombed -b ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed -a ${ASSOCIATIONFILE} --patient ${PATIENT} --pool ${POOL} --tag ${TAG} -d ${DBHOSTID} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}
# add comment to DB
mysql -u andrea --password=andrea -e "ALTER TABLE ${DBSCHEMA}.${DBTABLE} COMMENT = 'Pipe v.${PIPEVERSION}, last import at ${NOWIS}' ; " ;

### import data using the new DB structure
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import data into DB from BAM (new table different from the standard one)" 
import_iss --bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam --singleEnd --locusRead alnstart -p ${MAXTHREADS} -g ${CIGARGENOMEID} --dbhost 172.25.39.57 --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_refactored --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG};
# add comment to DB
mysql -u andrea --password=andrea -e "ALTER TABLE ${DBSCHEMA}.${DBTABLE}_refactored COMMENT = 'Pipe v.${PIPEVERSION}, last import at ${NOWIS}' ; " ;

### ------------------ ALIGNMENT STATS ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Statistics for reliable merged BAM (both fastqc and samstat)"
for k in $( ls ${TMPDIR}/bam/${BASENAME}.${TAG}*.bam ); do
	samstat $k;
done

### ------------------ PIPELINE (to DB) STATS ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Pipeline Statistics"

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";
LV_MAPPING_RESULTS=`samtools flagstat ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam ` ;
LV_MAPPING_PP=$((`echo ${LV_MAPPING_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
LV_MAPPING_ST=$((`echo ${LV_MAPPING_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
LV_MAPPING_OVERALL=$((${LV_MAPPING_ST}+${LV_MAPPING_PP})) ;
LV_MAPPING_MAPPED=$((`echo ${LV_MAPPING_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;

LTR_IDENTIFIED=`wc -l ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.sam.list | cut -d' ' -f1 ` ;
LTR_SEQAFTERTRIMMING=$((`zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
LC_SEQAFTERTRIMMING=$((`zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq.gz | wc -l `/4)) ;

BWA_RESULTS=`samtools flagstat ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam `
BWA_INPUT=$((`echo ${BWA_RESULTS} | grep ' in total ' | cut -d' ' -f1`/2)) ;
BWA_MAPPED=$((`echo ${BWA_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
BWA_MAPPED_PP=$((`echo ${BWA_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
BWA_MAPPED_ST=$((`echo ${BWA_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
BWA_MAPPED_OVERALL=$((${BWA_MAPPED_PP}+${BWA_MAPPED_ST})) ;
	
FILTER_QUAL_RESULTS=`samtools flagstat ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.bam ` ;
FILTER_QUAL_PP=$((`echo ${FILTER_QUAL_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
FILTER_QUAL_ST=$((`echo ${FILTER_QUAL_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
FILTER_QUAL_OVERALL=$((${FILTER_QUAL_ST}+${FILTER_QUAL_PP})) ;

FILTER_ISSR_RESULTS=`samtools flagstat ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam ` ;
FILTER_ISSR_PP=$((`echo ${FILTER_ISSR_RESULTS} | grep ' properly paired ' | cut -d' ' -f34`/2)) ;
FILTER_ISSR_ST=$((`echo ${FILTER_ISSR_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
FILTER_ISSR_OVERALL=$((${FILTER_ISSR_ST}+${FILTER_ISSR_PP})) ;

CIGARREADS_REMOVED=`zcat ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list.gz | wc -l | cut -d' ' -f1 ` ;
ISS_IDENTIFIED=`wc -l ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed | cut -d' ' -f1 ` ;

echo "
::VARIABLE SUMMARY:: dlimiters<>
	RUN_ID=<${RUN_ID}>
	RUN_NAME=<${RUN_NAME}>
	PHIX_MAPPING=<${PHIX_MAPPING}>
	BARCODE_MUX=<${BARCODE_MUX}>
	LTR_IDENTIFIED=<${LTR_IDENTIFIED}>
	LTR_SEQAFTERTRIMMING=<${LTR_SEQAFTERTRIMMING}>
	LC_SEQAFTERTRIMMING=<${LC_SEQAFTERTRIMMING}>
	BWA_RESULTS=<${BWA_RESULTS}>
	BWA_INPUT=<${BWA_INPUT}>
	BWA_MAPPED=<${BWA_MAPPED}>
	BWA_MAPPED_PP=<${BWA_MAPPED_PP}>
	BWA_MAPPED_ST=<${BWA_MAPPED_ST}>
	BWA_MAPPED_OVERALL=<${BWA_MAPPED_OVERALL}>
	FILTER_QUAL_RESULTS=<${FILTER_QUAL_RESULTS}>
	FILTER_QUAL_PP=<${FILTER_QUAL_PP}>
	FILTER_QUAL_ST=<${FILTER_QUAL_ST}>
	FILTER_QUAL_OVERALL=<${FILTER_QUAL_OVERALL}>
	FILTER_ISSR_RESULTS=<${FILTER_ISSR_RESULTS}>
	FILTER_ISSR_PP=<${FILTER_ISSR_PP}>
	FILTER_ISSR_ST=<${FILTER_ISSR_ST}>
	FILTER_ISSR_OVERALL=<${FILTER_ISSR_OVERALL}>
	CIGARREADS_REMOVED=<${CIGARREADS_REMOVED}>
	ISS_IDENTIFIED=<${ISS_IDENTIFIED}>
	LV_MAPPING_RESULTS=<${LV_MAPPING_RESULTS}>
	LV_MAPPING_PP=<${LV_MAPPING_PP}>
	LV_MAPPING_ST=<${LV_MAPPING_ST}>
	LV_MAPPING_OVERALL=<${LV_MAPPING_OVERALL}>
	LV_MAPPING_MAPPED=<${LV_MAPPING_MAPPED}>
" ;
	
	echo "[TIGET] Import STATS into SUMMARY table"
mysql -u andrea --password=andrea -e "INSERT INTO ${DBSCHEMA}.stats_summary_454 
(RUN_ID, RUN_NAME, DISEASE, PATIENT, POOL, TAG, LTR_ID, LC_ID, PHIX_MAPPING, BARCODE_MUX, LTR_IDENTIFIED, LTR_SEQAFTERTRIMMING, LC_SEQAFTERTRIMMING, BWA_RESULTS, BWA_INPUT, BWA_MAPPED, BWA_MAPPED_PP, BWA_MAPPED_ST, BWA_MAPPED_OVERALL, FILTER_QUAL_RESULTS, FILTER_QUAL_PP, FILTER_QUAL_ST, FILTER_QUAL_OVERALL, FILTER_ISSR_RESULTS, FILTER_ISSR_PP, FILTER_ISSR_ST, FILTER_ISSR_OVERALL, CIGARREADS_REMOVED, ISS_IDENTIFIED, LV_MAPPING_RESULTS, LV_MAPPING_PP, LV_MAPPING_ST, LV_MAPPING_OVERALL, LV_MAPPING_MAPPED)
VALUES
('${RUN_ID}', '${RUN_NAME}', '${DISEASE}', '${PATIENT}', '${POOL}', '${TAG}', '${TAG}', '${TAG}', '${PHIX_MAPPING}', '${BARCODE_MUX}', '${LTR_IDENTIFIED}', '${LTR_SEQAFTERTRIMMING}', '${LC_SEQAFTERTRIMMING}', '${BWA_RESULTS}', '${BWA_INPUT}', '${BWA_MAPPED}', '${BWA_MAPPED_PP}', '${BWA_MAPPED_ST}', '${BWA_MAPPED_OVERALL}', '${FILTER_QUAL_RESULTS}', '${FILTER_QUAL_PP}', '${FILTER_QUAL_ST}', '${FILTER_QUAL_OVERALL}', '${FILTER_ISSR_RESULTS}', '${FILTER_ISSR_PP}', '${FILTER_ISSR_ST}', '${FILTER_ISSR_OVERALL}', '${CIGARREADS_REMOVED}', '${ISS_IDENTIFIED}', '${LV_MAPPING_RESULTS}', '${LV_MAPPING_PP}', '${LV_MAPPING_ST}', '${LV_MAPPING_OVERALL}', '${LV_MAPPING_MAPPED}')
" ${DBSCHEMA} ;

# add comment to DB
NOWIS=`date +'%y-%m-%d %H:%M:%S'`
mysql -u andrea --password=andrea -e "ALTER TABLE ${DBSCHEMA}.${DBTABLE} COMMENT = 'Pipe v.${PIPEVERSION}, from ${STARTTIME} to ${NOWIS}' ; " ;


### ------------------ COPY RESULTS AND CLEANING ------------- ###
cp -r ${TMPDIR}/* ${SERVERWORKINGPATH}/${DISEASE}/${PATIENT}/${POOL}
# cp ${TMPDIR}/bam/${TAG}.sorted.md.*.ba* ${BASEPATH}/bam/
# cp ${TMPDIR}/${TAG}*.htm* ${BASEPATH}/bam/
# cp ${TMPDIR}/${OUTDIR}/${TAG}*.bed ${BASEPATH}/
# cp ${TMPDIR}/bam/${TAG}*.list ${BASEPATH}/


# a little final cleaning process
rm -fr ${TMPDIR}
