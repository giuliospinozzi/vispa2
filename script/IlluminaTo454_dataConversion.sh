#!/bin/bash -x
source /etc/environment
source /etc/profile

echo "
  +--------------------------------------------------------+
  |                                                        |
  |      Illumina to 454 Data Conversion (filtered)        |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Andrea Calabria                             |
  |            Giulio Spinozzi                             |
  |  Date:     Oct 2015                                    |
  |  Contact:  andrea.calabria@hsr.it                      |
  |            spinozzi.giulio@hsr.it                      |
  +--------------------------------------------------------+
"

STARTTIME=`date +'%y-%m-%d %H:%M:%S'`;

## CONFIG start ### --------------------------------------------
MAXTHREADS='12';

R1_FASTQ="${1}";
R2_FASTQ="${2}";
DISEASE="${3}";
PATIENT="${4}";
POOL="${5}";
ASSOCIATIONFILE="${6}";
BARCODE_LTR="${7}";
BARCODE_LC="${8}";
TMPDIR="${9}" # cartella di destinazione temporanea (che non viene cancellata)
LTR="${10}";
LC="${11}";
SEQCONTAM="${12}";
FASTQ_QF="${13}"; # no_quality_filter, slim, lam
CONVERSION_TYPE="${14}"; # gtris, vispa
## CONFIG end ### ---------------------------------------------

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Preprocessing input variables (delimiters:<>)";
## print input variables (check for log utils)
INPUTVARNUM=0;
for INPUTVAR in "$@"; do
	let INPUTVARNUM++; 
	printf -v INPUTNUM '%02d' $INPUTVARNUM;
	echo "  => Input Variable: Order number = <${INPUTNUM}> ; Var Content = <${INPUTVAR}>";
done

### NB: qui i nomi dei barcode sono SEMPRE nella forma LTRNN.LCNN, ovvero separati da un punto. Se non dovesse andar bene per qualsiasi motivo allora cambiate questo carattere delimitatore.
TAGDELIMITER="."
TAGDELIMITER_2WRITE=""

FASTADIR="${TMPDIR}/fasta" # cartella di destinazione finale dei file fasta, per ora e' uguale a tmp dir + fasta
## CONFIG end ### --------------------------------------------

mkdir ${TMPDIR}/joint
mkdir ${TMPDIR}/quality
mkdir ${FASTADIR}

BCLTRf=`cut "${BARCODE_LTR}" -f1 | sed 's/LTR//g'`;
BCLTR=(${BCLTRf}) ;
BCLCf=`cut "${BARCODE_LC}" -f1 | sed 's/LC//g'`;
BCLC=(${BCLCf}) ;

# identify base names
R1_NAME="`basename ${R1_FASTQ}`";
R2_NAME="`basename ${R2_FASTQ}`";

ASSOBCLIST=(`cut -f1 ${ASSOCIATIONFILE} | sort | uniq `); ## barcode list (as first colun of the association file!) user defined!

## FastQC quality
fastqc -o ${TMPDIR}/quality --nogroup --contaminants ${SEQCONTAM} -t ${MAXTHREADS} -f fastq ${R1_FASTQ} ${R2_FASTQ}

## remove low quality R1-R2 reads with fastq_qf program
if ! [ ${FASTQ_QF} = "no_quality_filter" ]; then
	## remove low quality R1-R2 reads with fastq_qf program
	fastq_qf -a ${R1_FASTQ} -b ${R2_FASTQ} -o ${TMPDIR} -t ${MAXTHREADS} -m ${FASTQ_QF}
	R1_FASTQ="${TMPDIR}/QF.${R1_NAME}";
	R2_FASTQ="${TMPDIR}/QF.${R2_NAME}";
	rm ${TMPDIR}/${R1_NAME} ${TMPDIR}/${R2_NAME};
fi

## remove first N bases
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Remove first 12 bases"
trimmomatic SE -phred33 -threads ${MAXTHREADS} ${R1_FASTQ} ${TMPDIR}/r1.no12.fastq.gz HEADCROP:12
trimmomatic SE -phred33 -threads ${MAXTHREADS} ${R2_FASTQ} ${TMPDIR}/r2.no12.fastq.gz HEADCROP:12
rm ${R1_FASTQ} ${R2_FASTQ}

# 1. demux barcode R1
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Demux barcodes of R1 -> LTR only"
fastq-multx -m 1 -B ${BARCODE_LTR} ${TMPDIR}/r1.no12.fastq.gz ${TMPDIR}/r2.no12.fastq.gz -o ${TMPDIR}/r1.no12.%.fastq.gz -o ${TMPDIR}/r2.no12.%.fastq.gz
rm ${TMPDIR}/r1.no12.fastq.gz ${TMPDIR}/r2.no12.fastq.gz

# 2. demux barcode R2 (revert pairs!!)
echo "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Demux barcodex of R2 -> LC only"
for k in ${BCLTR[@]}; do
	fastq-multx -m 1 -B ${BARCODE_LC} ${TMPDIR}/r2.no12.LTR${k}.fastq.gz ${TMPDIR}/r1.no12.LTR${k}.fastq.gz -o ${TMPDIR}/r2.no12.LTR${k}${TAGDELIMITER}%.fastq.gz -o ${TMPDIR}/r1.no12.LTR${k}${TAGDELIMITER}%.fastq.gz
	rm ${TMPDIR}/r1.no12.LTR${k}.fastq.gz ${TMPDIR}/r2.no12.LTR${k}.fastq.gz
done

# a questo punto vi trovate con la combinazione N*N di file paired con all'interno il nome del barcode nella forma LTRXY.LCZK (XY e ZK sono numeri interi di massimo 2 cifre). da qui per ogni coppia, fate il join
# 3. join. crea 3 files: 1) join, 2-3) un1/un2
if [ ${CONVERSION_TYPE} = "vispa" ]; then
	for b in ${BCLTR[@]}; do
		for k in ${BCLC[@]}; do 
			TAG="LTR${b}${TAGDELIMITER}LC${k}"
			echo $TAG
			for a in ${ASSOBCLIST[@]}; do 
				if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
					echo -e "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Join overlapping paired-ends"
					# ${TMPDIR}/r1.no12.${TAG}.fastq.gz ${TMPDIR}/r2.no12.${TAG}.fastq.gz
					echo "fastq-join <(zcat ${TMPDIR}/r1.no12.${TAG}.fastq.gz) <(zcat ${TMPDIR}/r2.no12.${TAG}.fastq.gz) -o ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_%.fastq"
					fastq-join <(zcat ${TMPDIR}/r1.no12.${TAG}.fastq.gz) <(zcat ${TMPDIR}/r2.no12.${TAG}.fastq.gz) -o ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_%.fastq
					

					# a questo punto i file sono tutti nella forma con suffissi:
					# _join.fastq.gz -> le unite correttamente
					# _un1.fastq.gz -> quelle derivanti da R1 non unite (e similmente _un2.fastq.gz quelle da R2 non unite)
					# noi andremo ad usare SOLO le join e le un1, che dovremo isolare ed unire in un unico file fasta

					cat ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_join.fastq | fastq_to_fasta -n -Q 33 -o ${FASTADIR}/${TAG}.fa ;
					cat ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_un1.fastq | fastq_to_fasta -n -Q 33 >> ${FASTADIR}/${TAG}.fa ;
					# NB: controllate che l'append qui sopra funzioni correttamente
					FINAL_TAG="LTR${b}${TAGDELIMITER_2WRITE}LC${k}"
					mv ${FASTADIR}/${TAG}.fa ${FASTADIR}/${FINAL_TAG}.fa
				fi
			done
		done
	done
fi

if [ ${CONVERSION_TYPE} = "gtris" ]; then
	for b in ${BCLTR[@]}; do
		for k in ${BCLC[@]}; do
			TAG="LTR${b}${TAGDELIMITER}LC${k}"
			echo $TAG
			for a in ${ASSOBCLIST[@]}; do
				if [ ${TAG} == ${a} ]; then echo "--- This tag exists in the association file -> analyze it ---";
					echo -e "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Trimmomatic on R1-R2, paired-ends"
					trimmomatic PE -phred33 -threads ${MAXTHREADS} ${TMPDIR}/r1.no12.${TAG}.fastq.gz ${TMPDIR}/r2.no12.${TAG}.fastq.gz ${TMPDIR}/r1.no12.${TAG}.paired.fastq.gz ${TMPDIR}/r1.no12.${TAG}.unpaired.fastq.gz ${TMPDIR}/r2.no12.${TAG}.paired.fastq.gz ${TMPDIR}/r2.no12.${TAG}.unpaired.fastq.gz ILLUMINACLIP:/opt/applications/bin/trimmomatic/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
					rm ${TMPDIR}/r1.no12.${TAG}.fastq.gz ${TMPDIR}/r2.no12.${TAG}.fastq.gz
					echo -e "<`date +'%Y-%m-%d %H:%M:%S'`> [TIGET] Join overlapping paired-ends"
					echo "fastq-join <(zcat ${TMPDIR}/r1.no12.${TAG}.paired.fastq.gz) <(zcat ${TMPDIR}/r2.no12.${TAG}.paired.fastq.gz) -o ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_%.fastq"
					fastq-join <(zcat ${TMPDIR}/r1.no12.${TAG}.paired.fastq.gz) <(zcat ${TMPDIR}/r2.no12.${TAG}.paired.fastq.gz) -o ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_%.fastq
					rm ${TMPDIR}/r1.no12.${TAG}.paired.fastq.gz ${TMPDIR}/r2.no12.${TAG}.paired.fastq.gz


					# a questo punto i file sono tutti nella forma con suffissi:
					# _join.fastq.gz -> le unite correttamente
					# _un1.fastq.gz -> quelle derivanti da R1 non unite (e similmente _un2.fastq.gz quelle da R2 non unite)
					# noi andremo ad usare SOLO le join e le un1, che dovremo isolare ed unire in un unico file fasta

					cat ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_join.fastq | fastq_to_fasta -n -Q 33 -o ${FASTADIR}/${TAG}.fa ;
					rm ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_join.fastq
					cat ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_un1.fastq | fastq_to_fasta -n -Q 33 >> ${FASTADIR}/${TAG}.fa ;
					rm ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_un1.fastq ${TMPDIR}/joint/${PATIENT}_${POOL}_${TAG}_un2.fastq
					zcat ${TMPDIR}/r1.no12.${TAG}.unpaired.fastq.gz | fastq_to_fasta -n -Q 33 >> ${FASTADIR}/${TAG}.fa ;

					# NB: controllate che l'append qui sopra funzioni correttamente
					FINAL_TAG="LTR${b}${TAGDELIMITER_2WRITE}LC${k}"
					mv ${FASTADIR}/${TAG}.fa ${FASTADIR}/${FINAL_TAG}.fa

					/opt/applications/scripts/isatk/script/./fasta_cleaner_wLTRwLC.sh -i ${FASTADIR}/${FINAL_TAG}.fa -l ${FASTADIR}/fasta_cleaner.${FINAL_TAG}.log -m ${LTR} -n ${LC}

					rm ${FASTADIR}/${FINAL_TAG}.kana.list.gz ${FASTADIR}/${FINAL_TAG}.trimLTR.trimLC.fa.gz ${FASTADIR}/${FINAL_TAG}.validLTR.list.gz ${FASTADIR}/${FINAL_TAG}.trimLTR.fa.gz
					mv ${FASTADIR}/${FINAL_TAG}.trimLTR.trimLC.noplasmid.fa.gz ${TMPDIR}
				fi
			done
		done
	done
fi

rm -rf ${TMPDIR}/joint
rm -rf ${TMPDIR}/fasta
rm -rf ${TMPDIR}/r1.no12.* ${TMPDIR}/r2.no12.*