#!/bin/bash

#    DO STATS FROM FASTQ
#    Copyright (C) 2014  A. Calabria (a.calabria@hsr.it)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

usage()
{
	echo "FASTQ Stats, ensemble of stat programs"
	echo
	echo "Usage: $0 -r R1_fastq -p R2_fastq -l log_file -o output_dir"
	echo
	echo "  [-l]            - Log file."
	echo "  [-r FASTQ R1]   - R1 fastq file, absolute path."
	echo "  [-p FASTQ R2]   - R2 fastq file, absolute path."
	echo "  [-o OUT DIR]    - Output dir."
	echo " At the end, you will find stats data in the log file and summaries in the output folder (FastQC, SamStat, FastX). The TMP dir is /tmp/stats"
	exit 
}

while getopts ":l:r:p:o:h" Option
	do
	case $Option in
		l ) LOGFILE="$OPTARG" ;;
		r ) R1_FASTQ="$OPTARG" ;;
		p ) R2_FASTQ="$OPTARG" ;;
		o ) OUTDIR_QUAL="$OPTARG" ;;
		h ) usage ;;
		* ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
	esac
done
shift $(($OPTIND - 1)) 

# if [ -z "$R1_FASTQ" ]; then
# 	usage
# fi

# if [ ! -r "$R1_FASTQ" ]; then
# 	echo "Error: can't open input file ($1)." >&2
# 	exit 1
# fi

# LOGFILE="${01}"
# R1_FASTQ="${02}"
# R2_FASTQ="${03}"
# OUTDIR_QUAL="${04}"

PHIXGENOME="/opt/genome/control/phix174/bwa_7/phiX174.fa"
TMPDIR="/tmp/stats"
MAXTHREADS="10"

mkdir ${TMPDIR}
chmod 775 ${TMPDIR}

fastqc -o ${OUTDIR_QUAL} -t ${MAXTHREADS} -f fastq ${R1_FASTQ} ${R2_FASTQ}

BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' `;
fastx_quality_stats -i <(zcat ${R1_FASTQ}) -o ${OUTDIR_QUAL}/${BNAME_R1}.statslog -Q 33
fastq_quality_boxplot_graph.sh -i ${OUTDIR_QUAL}/${BNAME_R1}.statslog -o ${OUTDIR_QUAL}/${BNAME_R1}.perbasequality.png -t "${BNAME_R1}"
fastx_nucleotide_distribution_graph.sh -i ${OUTDIR_QUAL}/${BNAME_R1}.statslog -o ${OUTDIR_QUAL}/${BNAME_R1}.nucleocomposition.png -t "${BNAME_R1}"

BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' `;
fastx_quality_stats -i <(zcat ${R2_FASTQ}) -o ${OUTDIR_QUAL}/${BNAME_R2}.statslog -Q 33
fastq_quality_boxplot_graph.sh -i ${OUTDIR_QUAL}/${BNAME_R2}.statslog -o ${OUTDIR_QUAL}/${BNAME_R2}.perbasequality.png -t "${BNAME_R2}"
fastx_nucleotide_distribution_graph.sh -i ${OUTDIR_QUAL}/${BNAME_R2}.statslog -o ${OUTDIR_QUAL}/${BNAME_R2}.nucleocomposition.png -t "${BNAME_R2}"

# samstat ${R1_FASTQ}
# mv ${R1_FASTQ}.samstat.html .

# samstat ${R2_FASTQ}
# mv ${R2_FASTQ}.samstat.html .

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] PHIX Alignment to reference genome each single pair"
echo "run BWA: bwa-7.5 mem -r 1 -M -T 15 -R \"@RG\tID:Test\tSM:PHIX\tCN:Andrea\" -t ${MAXTHREADS} ${PHIXGENOME} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) > ${TMPDIR}/${BNAME_R1}.phix.PE.sam"
bwa-7.5 mem -r 1 -M -T 15 -R "@RG\tID:Test\tSM:PHIX\tCN:Andrea" -t ${MAXTHREADS} ${PHIXGENOME} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) > ${TMPDIR}/${BNAME_R1}.phix.PE.sam
echo "run samtools: samtools view -F 268 -q 30 -f 35 -uS ${TMPDIR}/${BNAME_R1}.phix.PE.sam | samtools sort - ${TMPDIR}/${BNAME_R1}.phix.PE.sorted ; "
samtools view -F 268 -q 30 -f 35 -uS ${TMPDIR}/${BNAME_R1}.phix.PE.sam | samtools sort - ${TMPDIR}/${BNAME_R1}.phix.PE.sorted ;
rm ${TMPDIR}/${BNAME_R1}.phix.PE.sam
echo "run bamtools"
# non necessario perche' fatto qui sopra... ma ok
bamtools filter -in ${TMPDIR}/${BNAME_R1}.phix.PE.sorted.bam -isMateMapped true -isPaired true -isProperPair true | bedtools bamtobed | cut -f4 | cut -d'/' -f1 | sort | uniq > ${TMPDIR}/${BNAME_R1}.PHIX.header.list

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Checking Variables ";
PHIX_MAPPING=`wc -l ${TMPDIR}/${BNAME_R1}.PHIX.header.list | cut -d' ' -f1 ` ;
OVERALL_READS=$((`zcat ${R1_FASTQ} | wc -l | cut -d' ' -f1 `/4)) ;

echo "--- REPORT ---
Input File: ${R1_FASTQ}
PhiX mapping: ${PHIX_MAPPING}
Overall Reads: ${OVERALL_READS}
" >> ${LOGFILE}
