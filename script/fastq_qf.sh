#!/bin/bash

echo "
        +--------------------------------------------------------+
        |                                                        |
        |                  FastQ Quality Filter                  |
        |                                                        |
        +--------------------------------------------------------+
        |  Author:   Giulio Spinozzi, PhD student                |
        |  Date:     April 2016                              	 |
        |  Version:  1.0                                         |  
        |  Contact:  spinozzi.giulio@hsr.it                      |
        +--------------------------------------------------------+
"


##### ============================== RUN INFO =============================== #####
RUN_STARTED_AT=`date +"%Y-%m-%d %H:%M:%S"`;
RUN_ID="`whoami`"" ${RUN_STARTED_AT}";

TODAY=`date +"%Y%m%d%H%M"`;
#=================================================================================#


usage()
{
	echo "This app cleans FASTQ reads by low quality from R1 (<28 phred)"
	echo "low quality random barcodes from R2 - TAG (12n)"
	echo
	echo "Usage: $0 [-a R1.fastq.gz] [-b R2.fastq.gz] [-o output dir] [-t max threads] [-m method]"
	echo
	echo "  [-a R1.fastq.gz]   -   Input File: R1 FASTQ."
	echo "  [-b R2.fastq.gz]   -   Input File: R2 FASTQ."
	echo "  [-o output dir]    -   Output Directory"
	echo "  [-t max threads]   -   Maximum Number of Parallel Threads"
	echo "  [-m method]        -   Quality Filter on R1 only (lam) or R1 and R2 (slim)"
	echo
	exit 
}


while getopts ":a:b:o:t:m:h" Option
	do
	case $Option in
		a ) R1=$OPTARG ;;
		b ) R2=$OPTARG ;;
		o ) OUTDIR=$OPTARG ;;
		t ) MAXTHREADS=$OPTARG ;;
		m ) METHOD=$OPTARG ;;
		h ) usage ;;
		* ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
	esac
done
shift $(($OPTIND - 1)) 


##### ========================== PRELIMINARY CHECKS ========================= #####
## Arguments
if [ -z "$R1" ]; then
	usage
fi

if [ -z "$R2" ]; then
	usage
fi

if [ ! -r "$R1" ]; then
	echo "Error: can't open input file for R1."
	exit 1
fi

if [ ! -r "$R2" ]; then
	echo "Error: can't open input file for R2."
	exit 1
fi

if [ -z "$OUTDIR" ]; then
	usage
fi

if [ -z "$MAXTHREADS" ]; then
	MAXTHREADS=16
fi

if [ -z "$METHOD" ]; then
	METHOD="slim"
fi

## Space
LEFTSPACE=`df -P ${OUTDIR} | tail -1 | awk '{print $4}'`
if [ ${LEFTSPACE} -lt 100000000 ]; then 
    echo "
    
        *****************************************************************
        |                                                               |
        |   YOU DO NOT HAVE ENOUGH SPACE LEFT ON HARD DISK!!            |
        |                                                               |
        |   I WILL NOT PROCEED...                                       |
        |                                                               |       
        |   FREE SPACE BEFORE, LEAVING AT LEAST 100GB ON WD             |
        |                                                               |
        *****************************************************************
    
        "; 
    exit;
fi
#=================================================================================#


#---------------------------------------***---------------------------------------#
##### =============================== PROGRAM =============================== #####
echo "
---------------------------------------------------------------------------------
                       START PROCESSING AT: $RUN_STARTED_AT
---------------------------------------------------------------------------------
    " 

if [ ! -d "${OUTDIR}" ]; then
	mkdir ${OUTDIR};
fi


printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Setup\n\n"
R1_NAME="`basename ${R1}`";
R2_NAME="`basename ${R2}`";

#**** ============================= METHOD: SLiM ============================ ****#
if [ ${METHOD} = "slim" ]; then
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating R1 quality filtered reads list...\n"
	## remove low quality R1 reads
	trimmomatic SE -phred33 -threads $MAXTHREADS ${R1} "/dev/stdout" CROP:80 | fastq_quality_filter -q 28 -p 95 -Q 33 | fastq_to_fasta -Q 33 -n | fasta2csv | cut -d' ' -f1 > ${OUTDIR}/r1.quality.filtered.list
	printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating R2 quality filtered reads list...\n"
	## remove low quality R2 reads
	trimmomatic SE -phred33 -threads $MAXTHREADS ${R2} "/dev/stdout" CROP:12 | fastq_quality_filter -q 28 -p 100 -Q 33 | fastq_to_fasta -Q 33 -n | fasta2csv | cut -d' ' -f1 > ${OUTDIR}/r2.quality.filtered.list
	printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Intersecting R1 and R2 lists...\n\n"
	sort --parallel=5 ${OUTDIR}/r1.quality.filtered.list > ${OUTDIR}/r1.quality.filtered.sorted.list
	sort --parallel=5 ${OUTDIR}/r2.quality.filtered.list > ${OUTDIR}/r2.quality.filtered.sorted.list
	rm ${OUTDIR}/r1.quality.filtered.list ${OUTDIR}/r2.quality.filtered.list
	comm -12  <(cat ${OUTDIR}/r1.quality.filtered.sorted.list) <(cat ${OUTDIR}/r2.quality.filtered.sorted.list) > ${OUTDIR}/r1r2.quality.filtered.list
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Cleaning tmp files...\n\n"
	rm ${OUTDIR}/r1.quality.filtered.sorted.list ${OUTDIR}/r2.quality.filtered.sorted.list
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering R1 and R2 FastQ files...\n"
	zcat ${R1} | fqextract_pureheader <(cat ${OUTDIR}/r1r2.quality.filtered.list) | pigz --best -f -c > ${OUTDIR}/QF.${R1_NAME} &
	zcat ${R2} | fqextract_pureheader <(cat ${OUTDIR}/r1r2.quality.filtered.list) | pigz --best -f -c > ${OUTDIR}/QF.${R2_NAME} &
	wait
	rm ${OUTDIR}/r1r2.quality.filtered.list;
fi
#=================================================================================#


#**** ============================= METHOD: LAM ============================= ****#
if [ ${METHOD} = "lam" ]; then
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Creating R1 quality filtered reads list...\n"
	## remove low quality R1 reads
	trimmomatic SE -phred33 -threads $MAXTHREADS ${R1} "/dev/stdout" CROP:80 | fastq_quality_filter -q 28 -p 95 -Q 33 | fastq_to_fasta -Q 33 -n | fasta2csv | cut -d' ' -f1 > ${OUTDIR}/r1.quality.filtered.list
	printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Filtering R1 and R2 FastQ files...\n"
	zcat ${R1} | fqextract_pureheader <(cat ${OUTDIR}/r1.quality.filtered.list) | pigz --best -f -c > ${OUTDIR}/QF.${R1_NAME} &
	zcat ${R2} | fqextract_pureheader <(cat ${OUTDIR}/r1.quality.filtered.list) | pigz --best -f -c > ${OUTDIR}/QF.${R2_NAME} &
	wait
	rm ${OUTDIR}/r1.quality.filtered.list;
fi
#=================================================================================#


echo "
---------------------------------------------------------------------------------
                    END PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
---------------------------------------------------------------------------------
    " 
#=================================================================================#
#---------------------------------------***---------------------------------------#
