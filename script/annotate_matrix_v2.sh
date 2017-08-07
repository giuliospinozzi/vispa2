#!/bin/bash

echo "
        +--------------------------------------------------------+
        |                                                        |
        |               ISs Matrix Annotation V2                 |
        |                                                        |
        +--------------------------------------------------------+
        |  Author:   Giulio Spinozzi, PhD student                |
        |  Date:     April 2016                              	 |
        |  Version:  2.0                                         | 
        |  Contact:  spinozzi.giulio@hsr.it                      |
        +--------------------------------------------------------+
"


##### ============================== RUN INFO =============================== #####
RUN_STARTED_AT=`date +"%Y-%m-%d %H:%M:%S"`;
RUN_ID="`whoami`"" ${RUN_STARTED_AT}";

TODAY=`date +"%Y%m%d%H%M"`;
#=================================================================================#



##### ================================ ARGS ================================= #####
usage()
{
    echo "This app annotate ISs matrix file (csv, tsv)."
    echo
    echo "Usage: $0 [-m ISs.matrix.tsv] [-t type] [-g gtf_file.gtf] [-o output dir]"
    echo
    echo "  [-m ISs.matrix.tsv]   -   Input File: ISs Matrix (CSV, TSV)"
    echo "  [-t type]             -   Input Matrix Type: vispa or clustering"
    echo "  [-g gtf_file.gtf]     -   Input File: GTF File"
    echo "  [-o output dir]       -   Output Directory"
    echo
    exit 
}

while getopts ":g:m:o:t:h" Option
    do
    case $Option in
        g ) GTF=$OPTARG ;;
        m ) MATRIX=$OPTARG ;;
        o ) OUTDIR=$OPTARG ;;
        t ) TYPE=$OPTARG ;;
        h ) usage ;;
        * ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
    esac
done
shift $(($OPTIND - 1)) 
#=================================================================================#



##### ========================== PRELIMINARY CHECKS ========================= #####
## Arguments
if [ -z "$MATRIX" ]; then
    usage
fi

if [ -z "$GTF" ]; then
    usage
fi

if [ -z "$OUTDIR" ]; then
    usage
fi

if [ -z "$TYPE" ]; then
    usage
fi

if [ ! -f ${MATRIX} ]; then
    echo ""
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    echo " ${MATRIX} NOT EXISTS!!!!! "
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    exit 0
fi

if [ ! -f ${GTF} ]; then
    echo ""
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    echo " ${GTF} NOT EXISTS!!!!! "
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    exit 0
fi
#=================================================================================#


#---------------------------------------***---------------------------------------#
##### =============================== PROGRAM =============================== #####
echo "

---------------------------------------------------------------------------------
                    STARTING PROCESSING AT: $RUN_STARTED_AT
---------------------------------------------------------------------------------
    " 


##### =========================== SETUP PARAMETERS ========================== #####
mkdir ${OUTDIR};
#=================================================================================#

if [ ${TYPE} = "vispa" ]; then
    sed "s/\t0/\t/g" ${MATRIX} > ${MATRIX:0:-4}.no0.tsv;
    MATRIX="`basename ${MATRIX:0:-4}.no0.tsv`";

    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotation.\n"
    annotate_bed -a ${GTF} -b ${MATRIX} --matrix -o ${OUTDIR}/iss.annotation.tsv
    rm None;

    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Retrieving Annotation Data.\n\n"
    cat ${OUTDIR}/iss.annotation.tsv | cut -f7,11 > ${OUTDIR}/genes_data.tsv
    printf "GeneName\tGeneStrand\n" > ${OUTDIR}/genes.tsv
    cat ${OUTDIR}/genes_data.tsv >> ${OUTDIR}/genes.tsv

    printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Merging Annotation Data into ISs Matrix.\n"
    n_columns=$(head -n 1 ${MATRIX} | awk '{print NF}')
    cat ${MATRIX} | cut -f1,2,3 > ${OUTDIR}/matrix.part1.tsv
    cat ${MATRIX} | cut -f4-${n_columns} > ${OUTDIR}/matrix.part2.tsv
    MATRIX="`basename ${MATRIX:0:-4}`";
    paste -d '\t' ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/genes.tsv > ${OUTDIR}/matrix.part1.annotated.tsv
    paste -d '\t' ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/matrix.part2.tsv > ${OUTDIR}/${MATRIX}.annotated.tsv
    rm ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/matrix.part2.tsv ${OUTDIR}/genes.tsv ${OUTDIR}/genes_data.tsv ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/iss.annotation.tsv ${MATRIX:0:-4}.no0.tsv;

    echo "

    ---------------------------------------------------------------------------------
                          ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
    ---------------------------------------------------------------------------------
        " 
fi

if [ ${TYPE} = "clustering" ]; then
    sed 's/,$//g' ${MATRIX} > ${MATRIX:0:-4}.tmp.csv;
    sed 's|,|\t|g' ${MATRIX:0:-4}.tmp.csv > ${MATRIX:0:-4}.converted.tsv;
    rm ${MATRIX:0:-4}.tmp.csv;
    sed "s/\t0/\t/g" ${MATRIX:0:-4}.converted.tsv > ${MATRIX:0:-4}.no0.tsv;
    rm ${MATRIX:0:-4}.converted.tsv;
    MATRIX="`basename ${MATRIX:0:-4}.no0.tsv`";
    cat ${MATRIX} | awk -F$'\t' '{print $2"\t"$3"\t"$3}' | sed '1d' > ${MATRIX:0:-4}.bed

    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Annotation.\n"
    annotate_bed -a ${GTF} -b ${MATRIX:0:-4}.bed -o ${OUTDIR}/iss.annotation.tsv

    printf "\n<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Retrieving Annotation Data.\n\n"
    cat ${OUTDIR}/iss.annotation.tsv | cut -f7,11 > ${OUTDIR}/genes_data.tsv
    printf "GeneName\tGeneStrand\n" > ${OUTDIR}/genes.tsv
    cat ${OUTDIR}/genes_data.tsv >> ${OUTDIR}/genes.tsv

    printf "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Merging Annotation Data into ISs Matrix.\n"
    n_columns=$(head -n 1 ${MATRIX} | awk '{print NF}')
    cat ${MATRIX} | cut -f1,2,3,4 > ${OUTDIR}/matrix.part1.tsv
    cat ${MATRIX} | cut -f5-${n_columns} > ${OUTDIR}/matrix.part2.tsv
    MATRIX="`basename ${MATRIX:0:-4}`";
    paste -d '\t' ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/genes.tsv > ${OUTDIR}/matrix.part1.annotated.tsv
    paste -d '\t' ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/matrix.part2.tsv > ${OUTDIR}/${MATRIX}.annotated.tsv
    rm ${OUTDIR}/matrix.part1.tsv ${OUTDIR}/matrix.part2.tsv ${OUTDIR}/genes.tsv ${OUTDIR}/genes_data.tsv ${OUTDIR}/matrix.part1.annotated.tsv ${OUTDIR}/iss.annotation.tsv ${MATRIX:0:-4}.no0.tsv;

    echo "

    ---------------------------------------------------------------------------------
                          ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
    ---------------------------------------------------------------------------------
        " 
fi


#=================================================================================#
#---------------------------------------***---------------------------------------#