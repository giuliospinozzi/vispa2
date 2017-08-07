#!/bin/bash -x
source /etc/environment
source /etc/profile

PIPEVERSION="3.0.8d"
STARTTIME=`date +'%y-%m-%d %H:%M:%S'`
STATSDBNAME='stats_summary_454v38d'
DBHOSTIP="172.25.39.57"

echo "
=====================================================================================
		NGS pipeline for 454 Roche -> GAMMA AND LENTI RV
=====================================================================================
Date:	03rd December 2012 v0.1
		07th May 2013 v0.2
		27th August 2013 v0.2.1 -> reduced to accept 1 barcode per time and FastA
		04th April 2014 v1 -> BWA-MEM (ref to Illumina 3.0.4)
		12th December 2014 v3.0.7 -> BWA-MEM (ref to Illumina 3.0.8)
		25th January 2015 v3.0.8 -> new stats_summary, trimming optimized
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
MAXTHREADS="16"
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
mkdir ${TMPDIR}/iss
mkdir ${TMPDIR}/${OUTDIR}
# checking vars
RUN_ID=`date +"%Y%m%d%H%M%S"`
RUN_NAME="${DISEASE}|${PATIENT}|${POOL}|${TAG}"
PHIX_MAPPING="0"
BARCODE_MUX="0" # to add in case of barcoding demux
LTR_ID="NA"
LC_ID="NA"

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
	CREATE TABLE IF NOT EXISTS ${STATSDBNAME} (
	RUN_ID varchar(255) default null,
	RUN_NAME varchar(1000) default null,
	DISEASE varchar(255) default null,
	PATIENT varchar(255) default null,
	POOL varchar(255) default null,
	TAG varchar(255) default null,
	FASTA varchar(1000) default NULL,
	LTR_ID varchar(255) default null,
	LC_ID varchar(255) default null,
	PHIX_MAPPING int(30) DEFAULT NULL,
	BARCODE_MUX int(30) DEFAULT NULL,
	LTR_IDENTIFIED int(30) DEFAULT NULL,
	LTR_SEQAFTERTRIMMING int(30) DEFAULT NULL,
	LC_SEQAFTERTRIMMING int(30) DEFAULT NULL,
	LV_MAPPING_READS int(30) DEFAULT NULL,
	PLASMID_MAPPING_READS int(30) DEFAULT NULL,
	BWA_INPUT int(30) DEFAULT NULL,
	BWA_MAPPED int(30) DEFAULT NULL,
	FILTER_CIGARMD_TO_REMOVE int(30) DEFAULT NULL,
	FILTER_CIGARMD_REMAININGREADS int(30) DEFAULT NULL,
	FILTER_HQ_REMAININGREADS int(30) DEFAULT NULL,
	ISS_IDENTIFIED int(30) DEFAULT NULL,
	PRIMARY KEY (RUN_ID, TAG) 
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;
	" ${DBSCHEMA} ;


### ------------------ CONVERTING FASTA TO FASTQ ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Converting data FastA to FastQ"
if [[ ${FASTA} == *.gz ]]
	then
		zcat ${FASTA} | fasta_to_fastq -qv 9 > ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq
		BARCODE_MUX=`zcat ${FASTA} | grep ">" | wc -l | cut -d' ' -f1 `
	else
		# step 1 convert fasta into fastq, NB to fastq encoding!
		fasta_to_fastq -qv 9 ${FASTA} > ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq
		BARCODE_MUX=`cat ${FASTA} | grep ">" | wc -l | cut -d' ' -f1 `
fi


### ------------------ MARKING LTR READS ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Marking LTR"
# mark LTR reads and create list of only LTR reads (SAM accepted)
## ORIGINAL: # flexbar --source <(zcat ${TMPDIR}/${BASENAME}.fastq.gz) --target ${TMPDIR}/flexbar -f fastq-i1.8 -a ${LTR} --threads 8 -ae ANY -at 1 -ao 58 -aa -m 20 -q 10 --removal-tag ## this is for LENTI
flexbar2.5 --reads ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq --target ${TMPDIR}/reads/${BASENAME}.${TAG}.flexbar -f i1.8 -a ${LTR} --threads ${MAXTHREADS} -ae ANY -at 3.6 -ai -4 -ao 22 -m 2 -q 1 --removal-tags --max-uncalled 800 ## this is for LENTI

## HIV:: -at 3 -ai -10 -ao 22 -m 20 -q 10 --removal-tag

###flexbar --source <(zcat ${TMPDIR}/${BASENAME}.fastq.gz) --target ${TMPDIR}/flexbar -f fastq-i1.8 -a ${LTR} --threads 8 -ae ANY -at 1 -ao 53 -aa -m 20 -q 10 --removal-tag ## this is for GAMMA retro
grep "Flexbar" ${TMPDIR}/reads/${BASENAME}.${TAG}.flexbar.fastq > ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.only
sed -e 's/_Flexbar_removal//g' ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.only | sed -e 's/@//g' > ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.sam.list
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
flexbar2.5 --reads ${TMPDIR}/reads/${BASENAME}.${TAG}.fastq.gz --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR -f i1.8 -a ${LTR} --threads ${MAXTHREADS} -ae LEFT -at 3.6 -ai -4 -ao 22 -m 15 -q 1 --max-uncalled 800 -z GZ ## this is for HIV
### ------------------ end TRIMMING LTR ------------- ###

# flexbar --> you must then gzip files
gzip -f ${TMPDIR}/reads/${BASENAME}.*.fastq

### ------------------ TRIMMING LC ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Trimming LC"
#fastq-mcf ${LC} mld02.454.pool6.${TAG}.noLTR.fastq.gz -o mld02.454.pool6.${TAG}.noLTRLC.fastq.gz -m 12 -p 20 -l 30 -q 10 -P 33
## la precedente istruzione fallisce nel caso di più LC concatenate, esempio la read HHAUOBH02JFYTO
flexbar2.5 --reads ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR.fastq.gz --target ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC -f i1.8 -a ${LC} --threads ${MAXTHREADS} -ae RIGHT -at 4 -ao 8 -ai -4 -m 2 -q 5 --max-uncalled 800 -z GZ ;

# flexbar --> you must then gzip files
gzip -f ${TMPDIR}/reads/${BASENAME}.*.noLTRLC.fastq

LTR_IDENTIFIED=`wc -l ${TMPDIR}/reads/${BASENAME}.${TAG}.LTR.sam.list | cut -d' ' -f1 ` ;
LTR_SEQAFTERTRIMMING=$((`zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTR.fastq.gz | wc -l | cut -d' ' -f1 `/4)) ;
LC_SEQAFTERTRIMMING=$((`zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq.gz | wc -l `/4)) ;

    
##### ================ LV removal TAG based - start ======================== #####
OUTDIR_VECTOR_POOL="${TMPDIR}/lv";
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference VECTOR genome"
bwa-stable mem -k 14 -r 1 -T 15 -c 1 -R "@RG\tID:${TAG}\tSM:Vector\tCN:${DISEASE}.${PATIENT}.${POOL}.Andrea" -t ${MAXTHREADS} ${VECTORGENOME} <(zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq.gz) > ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sam

# create BAM and sort them
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter LV from here using only valid reads: mapped and primary)"
# samtools view -F 260 -q 5 -uS ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sam > ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.bam ;
# samtools sort ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.bam ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted;
samtools view -F 2308 -q 40 -uS ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sam | samtools sort - ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted;
samtools index ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bam;
samtools fillmd -b ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bam ${VECTORGENOME} > ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam
samtools index ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam ;

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filter by CIGAR and MD tag (only for starting position and CIGAR value) and then select only LV reads from BAM."
filter_by_cigar_bam --bam ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam -o ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list --minStartingMatches_fromNbp 2 --minStartingMatches 14 -p ${MAXTHREADS} -g ${VECTORCIGARGENOMEID} --minStartingBasesNoIndels 16 --endClipThreshold 100 --singleEnd

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

# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Statistics for LV BAM (both fastqc and samstat)"
# samstat ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam 

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Import data into DB from BAM (new table different from the standard one)" 
import_iss --bam ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam --singleEnd --locusRead alnstart -p ${MAXTHREADS} -g ${VECTORCIGARGENOMEID} --dbhost ${DBHOSTIP} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_lv --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG};

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Extract no LV reads from raw data"
zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.noLTRLC.fastq.gz | fqreverseextract_pureheader ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list | gzip -f > ${TMPDIR}/reads/${BASENAME}.${TAG}.no12.noLTRLC.noVector.fastq.gz 

# clean lv folder
rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sam;
rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.bam;
#rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam
rm ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bam ;
gzip -f ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.bfilter.toremovebycigar.list

LV_MAPPING_RESULTS=`samtools flagstat ${TMPDIR}/lv/${BASENAME}.${TAG}.noLTRLC.sorted.md.lvonly.bam ` ;
LV_MAPPING_READS=$((`echo ${LV_MAPPING_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;
#LV_MAPPING_READS=`samtools view ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.sorted.bam | cut -f1 | sort | uniq | wc -l | cut -d' ' -f1 `


#### ================ LV removal TAG based - end ======================== #####

##### ================ PLASMIDS quantification - start ======================== #####
INPUTFASTA="${TMPDIR}/reads/${BASENAME}.${TAG}.no12.noLTRLC.noVector.fastq.gz"
PLASMID="kana"
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Quantify and clean by ${PLASMID} reads"
# low seed value due to small genome
bwa-stable mem -k 14 -r 1 -T 15 -c 1 -R "@RG\tID:plasmid\tSM:${PLASMID}\tCN:Andrea" -t ${MAXTHREADS} /opt/genome/vector/plasmid/plasmids.${PLASMID}.fa ${INPUTFASTA} > ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.sam
# quality filter of seed-1 value
samtools view -F 2308 -q 20 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.sam | samtools sort - ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.sorted
rm ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.sam
# fastqc -o ${OUTDIR_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${BASEDIR}/${DISEASE}_${PATIENT}_${POOL}.LTR${b}.LC${k}.noLTRLC.${PLASMID}.merge.bam
samtools view ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.sorted.bam | cut -f1 | sort | uniq | gzip -f > ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.list.gz
rm ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.sorted.bam
PLASMID_MAPPING_READS=`zcat ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.list.gz | wc -l  | cut -d' ' -f1 `
# extract only NO PHIX reads
zcat ${INPUTFASTA} | fqreverseextract_pureheader <( zcat ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.noLV.${PLASMID}.list.gz ) | gzip -f > ${TMPDIR}/reads/${BASENAME}.${TAG}.no12.noLTRLC.noVector.noPlasmid.fastq.gz
##### ================ PLASMIDS quantification - end ======================== #####

### ------------------ ALIGNMENT ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Alignment to reference genome"
bwa-stable mem -k 18 -r 1 -M -T 15 -c 1 -R "@RG\tID:${DISEASE}.${PATIENT}.${POOL}.${TAG}\tSM:${TAG}\tCN:Andrea.${DISEASE}.${PATIENT}.${POOL}" -t ${MAXTHREADS} ${GENOME} <(zcat ${TMPDIR}/reads/${BASENAME}.${TAG}.no12.noLTRLC.noVector.noPlasmid.fastq.gz ) > ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam

# create BAM and sort them
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating BAM and indexes (filter from here the dataset using only valid reads: mapped and primary)"
samtools view -F 2308 -uS ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam | samtools sort - ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted;
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.bam ;
chmod -R 755 ${TMPDIR}/bam/*
#rm ${TMPDIR}/sam/${BASENAME}.${TAG}.noLTRLC.sam
rm ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.bam

### ------------------ RECALIBRATION ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Recreate/Recalibrate/Fill the MD tag"
samtools fillmd -b ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.bam ${GENOME} > ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam 
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam
rm ${TMPDIR}/bam/${BASENAME}.${TAG}.sorted.bam

BWA_RESULTS=`samtools flagstat ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam `
BWA_INPUT=$((`echo ${BWA_RESULTS} | grep ' in total ' | cut -d' ' -f1`)) ;
BWA_MAPPED=$((`echo ${BWA_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;

### ------------------ FILTERING ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filter by CIGAR and MD tag"
filter_by_cigar_bam --bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.bam -o ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list --minStartingMatches_fromNbp 2 --minStartingMatches 5 -p ${MAXTHREADS} -g ${CIGARGENOMEID} --minStartingBasesNoIndels 4 --compareSubOptimal --suboptimalThreshold ${SUBOPTIMALTHRESHOLD} --SAalnAction ignore --XSlikeTag XS --ASlikeTag AS --endClipThreshold 1000 --singleEnd 

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
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam

FILTER_CIGARMD_TO_REMOVE=`wc -l ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.toremovebycigar.list | cut -d' ' -f1 `
FILTER_CIGARMD_RESULTS=`samtools flagstat ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam ` ;
FILTER_CIGARMD_REMAININGREADS=$((`echo ${FILTER_CIGARMD_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;


echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering data (Bamtools)"
bamtools filter -in ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam -mapQuality ">=12" -out ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.hq.bam

# rename files and index them
mv ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.lq.bam
mv ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.hq.bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.lq.bam
samtools index ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam

FILTER_HQ_RESULTS=`samtools flagstat ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam ` ;
FILTER_HQ_REMAININGREADS=$((`echo ${FILTER_HQ_RESULTS} | grep '0 mapped ' | cut -d' ' -f15`)) ;

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
import_iss --bam ${TMPDIR}/bam/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bam --singleEnd --locusRead alnstart -p ${MAXTHREADS} -g ${CIGARGENOMEID} --dbhost ${DBHOSTIP} --dbschema ${DBSCHEMA} --dbtable ${DBTABLE}_refactored --dbuser andrea --dbpassword andrea -a ${ASSOCIATIONFILE} --tag ${TAG} --outfile ${TMPDIR}/iss/${DBSCHEMA}_${DBTABLE}_refactored.tsv;
# add comment to DB
mysql -u andrea --password=andrea -e "ALTER TABLE ${DBSCHEMA}.${DBTABLE}_refactored COMMENT = 'Pipe v.${PIPEVERSION}, last import at ${NOWIS}' ; " ;

### ------------------ ALIGNMENT STATS ------------- ###
# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Statistics for reliable merged BAM (both fastqc and samstat)"
# for k in $( ls ${TMPDIR}/bam/${BASENAME}.${TAG}*.bam ); do
# 	samstat $k;
# done

### ------------------ PIPELINE (to DB) STATS ------------- ###
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Pipeline Statistics"

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";


ISS_IDENTIFIED=`wc -l ${TMPDIR}/bed/${BASENAME}.${TAG}.noLTRLC.sorted.md.filter.iss.bed | cut -d' ' -f1 ` ;

echo "
::VARIABLE SUMMARY:: dlimiters<>
	RUN_ID=<${RUN_ID}>
	RUN_NAME=<${RUN_NAME}>
	DISEASE=<${DISEASE}>
	PATIENT=<${PATIENT}>
	POOL=<${POOL}>
	TAG=<${TAG}>
	FASTA=<${FASTA}>
	LTR_ID=<${LTR_ID}>
	LC_ID=<${LC_ID}>
	PHIX_MAPPING=<${PHIX_MAPPING}>
	BARCODE_MUX=<${BARCODE_MUX}>
	LTR_IDENTIFIED=<${LTR_IDENTIFIED}>
	LTR_SEQAFTERTRIMMING=<${LTR_SEQAFTERTRIMMING}>
	LC_SEQAFTERTRIMMING=<${LC_SEQAFTERTRIMMING}>
	LV_MAPPING_READS=<${LV_MAPPING_READS}>
	PLASMID_MAPPING_READS=<${PLASMID_MAPPING_READS}>
	BWA_INPUT=<${BWA_INPUT}>
	BWA_MAPPED=<${BWA_MAPPED}>
	FILTER_CIGARMD_TO_REMOVE=<${FILTER_CIGARMD_TO_REMOVE}>
	FILTER_CIGARMD_REMAININGREADS=<${FILTER_CIGARMD_REMAININGREADS}>
	FILTER_HQ_REMAININGREADS=<${FILTER_HQ_REMAININGREADS}>
	ISS_IDENTIFIED=<${ISS_IDENTIFIED}>
" ;
	
	echo "[TIGET] Import STATS into SUMMARY table"
mysql -u andrea --password=andrea -e "INSERT INTO ${DBSCHEMA}.${STATSDBNAME} 
(RUN_ID, RUN_NAME, DISEASE, PATIENT, POOL, TAG, FASTA, LTR_ID, LC_ID, PHIX_MAPPING, BARCODE_MUX, LTR_IDENTIFIED, LTR_SEQAFTERTRIMMING, LC_SEQAFTERTRIMMING, LV_MAPPING_READS, PLASMID_MAPPING_READS, BWA_INPUT, BWA_MAPPED, FILTER_CIGARMD_TO_REMOVE, FILTER_CIGARMD_REMAININGREADS, FILTER_HQ_REMAININGREADS, ISS_IDENTIFIED)
VALUES
('${RUN_ID}', '${RUN_NAME}', '${DISEASE}', '${PATIENT}', '${POOL}', '${TAG}', '${FASTA}', '${LTR_ID}', '${LC_ID}', '${PHIX_MAPPING}', '${BARCODE_MUX}', '${LTR_IDENTIFIED}', '${LTR_SEQAFTERTRIMMING}', '${LC_SEQAFTERTRIMMING}', '${LV_MAPPING_READS}', '${PLASMID_MAPPING_READS}', '${BWA_INPUT}', '${BWA_MAPPED}', '${FILTER_CIGARMD_TO_REMOVE}', '${FILTER_CIGARMD_REMAININGREADS}', '${FILTER_HQ_REMAININGREADS}', '${ISS_IDENTIFIED}')
" ${DBSCHEMA} ;

# add comment to DB
NOWIS=`date +'%y-%m-%d %H:%M:%S'`
mysql -u andrea --password=andrea -e "ALTER TABLE ${DBSCHEMA}.${DBTABLE} COMMENT = 'Pipe v.${PIPEVERSION}, from ${STARTTIME} to ${NOWIS}' ; " ;


### ------------------ COPY RESULTS AND CLEANING ------------- ###
rm -fr ${TMPDIR}/sam/*.sam
cp -r ${TMPDIR}/* ${SERVERWORKINGPATH}/${DISEASE}/${PATIENT}/${POOL}
# cp ${TMPDIR}/bam/${TAG}.sorted.md.*.ba* ${BASEPATH}/bam/
# cp ${TMPDIR}/${TAG}*.htm* ${BASEPATH}/bam/
# cp ${TMPDIR}/${OUTDIR}/${TAG}*.bed ${BASEPATH}/
# cp ${TMPDIR}/bam/${TAG}*.list ${BASEPATH}/


# a little final cleaning process
rm -fr ${TMPDIR}
