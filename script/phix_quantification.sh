#!/bin/bash
source /etc/environment
source /etc/profile

echo "
  +--------------------------------------------------------+
  |                                                        |
  |                   PhiX quantification                  |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Giulio Spinozzi                             |
  |  Date:     Aug 2015                                    |
  |  Contact:  spinozzi.giulio@hsr.it                      |
  |  Version:  1.3.1.20 (BWA MEM) aka 3.0 -> latest 3.0.6  |
  |            3.0.5 Fixed bugs in CIGAR and filters       |
  |            3.0.6 Trimming optimized with LTR Analysis  |
  +--------------------------------------------------------+
"

usage()
{
	echo "This app quantify PHIX in a specific seq run"
	echo
	echo "Usage: $0"
	echo
	echo "  [-a R1.fastq.gz]  - R1 FASTQ (compressed)"
	echo "  [-b R2.fastq.gz]  - R2 FASTQ (compressed)"
	echo "  [-c phiX174.fa]   - phiX174 fa"
	echo "  [-d POOL]         - Sequencing POOL"
	echo "  [-e /tmp]         - Output DIR"
	echo
	exit 
}

while getopts ":a:b:c:d:e:h" Option
	do
	case $Option in
		a ) R1_FASTQ="$OPTARG" ;;
		b ) R2_FASTQ="$OPTARG" ;;
		c ) PHIXGENOME="$OPTARG" ;;
		d ) POOL="$OPTARG" ;;
		e ) TMPDIR="$OPTARG" ;;
		h ) usage ;;
		* ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
	esac
done
shift $(($OPTIND - 1)) 


# ##### ================ PHIX quantification START ======================== #####
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] PHIX Alignment to reference genome each single pair"
bwa-stable mem -k 16 -r 1 -M -T 15 -R "@RG\tID:Test\tSM:PHIX\tCN:TIGET" -t 16 ${PHIXGENOME} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} )> ${TMPDIR}/${POOL}.phix.PE.sam
samtools view -F 2308 -q 25 -f 35 -uS ${TMPDIR}/${POOL}.phix.PE.sam > ${TMPDIR}/${POOL}.phix.PE.bam;
rm ${TMPDIR}/${POOL}.phix.PE.sam
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering raw data"
samtools view ${TMPDIR}/${POOL}.phix.PE.bam | cut -f 1 | sort | uniq > ${TMPDIR}/PHIX.${POOL}.list
rm ${TMPDIR}/${POOL}.phix.PE.bam;

PHIX_MAPPING=`wc -l ${TMPDIR}/PHIX.${POOL}.list | cut -d' ' -f1 ` ;
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] PHIX Aligned in POOL (${POOL}): ${PHIX_MAPPING}"

# ##### ================ PHIX quantification END  ======================== #####