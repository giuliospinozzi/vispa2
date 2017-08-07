#!/bin/bash -x
### --- This program will: ----
## quantify and prune phix
## quantifY/identifY sequences with ltr
## quantify and prune plasmids
## log results (in append!!) to an output file
### ---------------
echo "
Andrea Calabria - Vector Integration Core TIGET
Giulio Spinozzi - Vector Integration Core TIGET
 Mail:     a.calabria@hsr.it
           spinozzi.giulio@hsr.it
 Version:  1    (2015 Jan 14)
           1.1  (2015 Jan 15, No PhiX removal)
"

usage()
{
	echo "This app cleans FASTA reads by PHIX, non LTR reads and PLASMIDS, and returns"
	echo "trimmed sequences in FASTA (gzip compressed, in the input dir)."
	echo
	echo "Usage: $0 [-i INPUT.fa] [-l LOGFILE]"
	echo
	echo "  [-i INPUT.fa]   - Input file FASTA (uncompressed)."
	echo "  [-l LOGFILE]    - Output log file, in append."
	echo "  [-m LTR]    - LTR fa"
	echo "  [-n LC]    - LC fa"
	echo
	exit 
}

while getopts ":i:l:m:n:h" Option
	do
	case $Option in
		i ) FASTA="$OPTARG" ;;
		l ) LOGFILE=$OPTARG ;;
		m ) LTR=$OPTARG ;;
		n ) LC=$OPTARG ;;
		h ) usage ;;
		* ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
	esac
done
shift $(($OPTIND - 1)) 

if [ -z "$FASTA" ]; then
	usage
fi

if [ ! -r "$FASTA" ]; then
	echo "Error: can't open input file ($1)." >&2
	exit 1
fi

# FASTA="${1}"
# DDIR="${2}" ; #"/home/andrea/test/shearsites"
# LOGFILE="${3}"

DDIR=`dirname ${FASTA}`
MAXTHREADS=10
PHIXGENOME="/opt/genome/control/phix174/bwa_7/phiX174.fa"
LTR=${LTR};
LC=${LC};
#LTR="/opt/applications/scripts/isatk/elements/sequences/LTR.32bp.fa"
#LC="/opt/applications/scripts/isatk/elements/sequences/LC.rc.fa"

RAW_READS=$((`cat ${FASTA} | grep ">" | wc -l | cut -d' ' -f1 `)) ;
BNAME_FASTA=`basename ${FASTA} | cut -d'.' -f1`;


# NOTE: PHIX quantification (here below) is not reliable because is based only on alignment quality (illumina may include phix/background bases to extend/complete the read up to 250/300 cycles, thus many reads map phix at the end and you will discard them by mistake)
# ##### ================ PHIX quantification ======================== #####
# echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Quantify and clean by PHIX reads"
# bwa-7.10 mem -k 16 -r 1 -T 15 -c 1 -R "@RG\tID:Test\tSM:PHIX\tCN:Andrea" -t ${MAXTHREADS} ${PHIXGENOME} ${FASTA} > ${DDIR}/${BNAME_FASTA}.phix.sam
# # quality filter of seed-1 value
# samtools view -F 2308 -q 40 -uS ${DDIR}/${BNAME_FASTA}.phix.sam | samtools sort - ${DDIR}/${BNAME_FASTA}.phix.sorted
# rm ${DDIR}/${BNAME_FASTA}.phix.sam
# # fastqc -o ${OUTDIR_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${BASEDIR}/${DISEASE}_${PATIENT}_${POOL}.LTR${b}.LC${k}.noLTRLC.${PLASMID}.merge.bam
# samtools view ${DDIR}/${BNAME_FASTA}.phix.sorted.bam | cut -f1 | sort | uniq | gzip -f > ${DDIR}/${BNAME_FASTA}.phix.list.gz
# rm ${DDIR}/${BNAME_FASTA}.*.sorted.bam

# PHIX_MAPPING_READS=`zcat ${DDIR}/${BNAME_FASTA}.phix.list.gz | wc -l  | cut -d' ' -f1 `
# # extract only NO PHIX reads
# fasta_to_fastq -q 9 ${FASTA} | fastq_to_fasta -n | fareverseextract_pureheader <( zcat ${DDIR}/${BNAME_FASTA}.phix.list.gz ) | gzip -f > ${DDIR}/${BNAME_FASTA}.nophix.fa.gz


##### ================ LTR identification and quantification ======================== #####
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] LTR: identify seqs and trim LTR and LC"
# INPUTFASTA="${DDIR}/${BNAME_FASTA}.nophix.fa.gz"
INPUTFASTA="${FASTA}"
# identify LTR
flexbar2.5 --reads ${INPUTFASTA} -f i1.8 -a ${LTR} --threads ${MAXTHREADS} -ae LEFT -at 3.6 -ai -4 -ao 22 -m 15 -q 1 --stdout-reads --removal-tags --max-uncalled 800 | grep Flexbar | cut -d' ' -f 1 | sed 's/>//g' | sed 's/_Flexbar_removal//g' | pigz --best -f -c > ${DDIR}/${BNAME_FASTA}.validLTR.list.gz
# extract only LTR identified reads
# # in case of compressed and trasformed (no reads split in multiple lines):: zcat ${INPUTFASTA} | faextract_pureheader <( zcat ${DDIR}/${BNAME_FASTA}.validLTR.list.gz ) | gzip -f > ${DDIR}/${BNAME_FASTA}.validLTR.fa.gz
fasta_to_fastq -q 9 ${FASTA} | fastq_to_fasta -n | faextract_pureheader <( zcat ${DDIR}/${BNAME_FASTA}.validLTR.list.gz ) | pigz --best -f -c > ${DDIR}/${BNAME_FASTA}.validLTR.fa.gz
# trim LTR
flexbar2.5 --reads ${DDIR}/${BNAME_FASTA}.validLTR.fa.gz --target ${DDIR}/${BNAME_FASTA}.trimLTR -f i1.8 -a ${LTR} --threads ${MAXTHREADS} -ae LEFT -at 3.6 -ai -4 -ao 22 -m 15 -q 1 -z GZ --max-uncalled 800 
# trim LC
flexbar2.5 --reads ${DDIR}/${BNAME_FASTA}.trimLTR.fasta.gz --target ${DDIR}/${BNAME_FASTA}.trimLTR.trimLC -f i1.8 -a ${LC} --threads ${MAXTHREADS} -ae RIGHT -at 4 -ao 8 -m 2 -q 5 -z GZ --max-uncalled 800 ;
# rename file
mv ${DDIR}/${BNAME_FASTA}.trimLTR.trimLC.fasta.gz ${DDIR}/${BNAME_FASTA}.trimLTR.trimLC.fa.gz
mv ${DDIR}/${BNAME_FASTA}.trimLTR.fasta.gz ${DDIR}/${BNAME_FASTA}.trimLTR.fa.gz
rm ${DDIR}/${BNAME_FASTA}.validLTR.fa.gz

LTR_IDENTIFIED_READS=`zcat ${DDIR}/${BNAME_FASTA}.validLTR.list.gz | wc -l | cut -d' ' -f1 `
READS_POST_TRIMMING=`zcat ${DDIR}/${BNAME_FASTA}.trimLTR.trimLC.fa.gz | grep ">" | wc -l | cut -d' ' -f1 `


##### ================ PLASMIDS quantification ======================== #####
INPUTFASTA="${DDIR}/${BNAME_FASTA}.trimLTR.trimLC.fa.gz"
PLASMID="kana"
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Quantify and clean by ${PLASMID} reads"
# low seed value due to small genome
bwa-stable mem -k 14 -r 1 -T 15 -c 1 -R "@RG\tID:plasmid\tSM:${PLASMID}\tCN:Andrea" -t ${MAXTHREADS} /opt/genome/vector/plasmid/plasmids.${PLASMID}.fa ${INPUTFASTA} > ${DDIR}/${BNAME_FASTA}.${PLASMID}.sam
# quality filter of seed-1 value
samtools view -F 2308 -q 20 -uS ${DDIR}/${BNAME_FASTA}.${PLASMID}.sam | samtools sort - ${DDIR}/${BNAME_FASTA}.${PLASMID}.sorted
rm ${DDIR}/${BNAME_FASTA}.${PLASMID}.sam
# fastqc -o ${OUTDIR_QUAL} --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f bam ${BASEDIR}/${DISEASE}_${PATIENT}_${POOL}.LTR${b}.LC${k}.noLTRLC.${PLASMID}.merge.bam
samtools view ${DDIR}/${BNAME_FASTA}.${PLASMID}.sorted.bam | cut -f1 | sort | uniq | pigz --best -f -c > ${DDIR}/${BNAME_FASTA}.${PLASMID}.list.gz
rm ${DDIR}/${BNAME_FASTA}.*.sorted.bam
PLASMID_MAPPING_READS=`zcat ${DDIR}/${BNAME_FASTA}.${PLASMID}.list.gz | wc -l  | cut -d' ' -f1 `
# extract only NO PHIX reads
zcat ${INPUTFASTA} | fareverseextract_pureheader <( zcat ${DDIR}/${BNAME_FASTA}.${PLASMID}.list.gz ) | pigz --best -f -c > ${DDIR}/${BNAME_FASTA}.trimLTR.trimLC.noplasmid.fa.gz

FINAL_READS=`zcat ${DDIR}/${BNAME_FASTA}.trimLTR.trimLC.noplasmid.fa.gz | grep ">" | wc -l | cut -d' ' -f1 `


##### ================ LOG and PRINT STATS ======================== #####
# echo -e "FASTA\tBNAME_FASTA\tRAW_READS\tPHIX_MAPPING_READS\tLTR_IDENTIFIED_READS\tREADS_POST_TRIMMING\tPLASMID_MAPPING_READS\tFINAL_READS" >> ${LOGFILE}
# echo -e "${FASTA}\t${BNAME_FASTA}\t${RAW_READS}\t${PHIX_MAPPING_READS}\t${LTR_IDENTIFIED_READS}\t${READS_POST_TRIMMING}\t${PLASMID_MAPPING_READS}\t${FINAL_READS}" >> ${LOGFILE}

echo -e "FASTA\tBNAME_FASTA\tRAW_READS\tLTR_IDENTIFIED_READS\tREADS_POST_TRIMMING\tPLASMID_MAPPING_READS\tFINAL_READS" >> ${LOGFILE}
echo -e "${FASTA}\t${BNAME_FASTA}\t${RAW_READS}\t${LTR_IDENTIFIED_READS}\t${READS_POST_TRIMMING}\t${PLASMID_MAPPING_READS}\t${FINAL_READS}" >> ${LOGFILE}


