#!/bin/bash

source /etc/environment
source /etc/profile

echo "
  +--------------------------------------------------------+
  |                                                        |
  |       Data Formatting Pipeline:                        |
  |       - From MiSeq data to 454 data (FastA)            |
  |         formatted to run with the OLD TIGET pipeline.  |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Andrea Calabria                             |
  |  Date:     January 2013                                |
  |  Contact:  andrea.calabria@hsr.it                      |
  +--------------------------------------------------------+

  REQUIRED VARS and relative ORDER POSITION -> REMEMBER NO SPACES!!!!!!!!!
	1. disease id (MLD)
	2. patient id (MLD02)
	3. local server working path [/opt/NGS/results] -> it will be added DISEASE and PATIENT ID [/opt/NGS/results/MLD/MLD02/]
	4. R1
	5. R2
	6. pool id
	7. barcode list
	8. reference genome (NB: indexed by BWA, hg19:: /opt/NGS/index/human/hg19/bwa/hg19.nochr0)
	9. tmp dir
	
"

ADAPTER454_A="/home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454adapters.A.fa";
ADAPTER454_B="/home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454adapters.B.fa";
LTR="/home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LTR.fa";
LC="/home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LC.fa";

R1="${4}"
R2="${5}"
POOLNAME="${6}";
PATIENT="${2}";
GENOME="${8}"; # the right one
DISEASE="${1}";
NGSWORKINGPATH="${3}"; # WITHOUT LAST SLASH -> if not present, I will create it a new one with this name using mkdir
TMPDIR="${9}" ;
POOL="${POOLNAME}";
BARCODELIST="${7}";
BCF=`cat "$BARCODELIST" | cut -f1`; 
BARCODEARRAY=( ${BCF} ); # array of barcodes only
BASEDIR="${NGSWORKINGPATH}/${DISEASE}/${PATIENT}/miseqTo454";
OUTDIR="${BASEDIR}/${POOL}"
FASTADIR="${OUTDIR}/fasta"
FASTADIR_MERGE="${FASTADIR}/merge"
mkdir ${TMPDIR};
mkdir ${NGSWORKINGPATH}/${DISEASE} ;
mkdir ${NGSWORKINGPATH}/${DISEASE}/${PATIENT} ;
mkdir ${BASEDIR} ;
mkdir ${OUTDIR} ;
mkdir ${FASTADIR} ;
mkdir ${FASTADIR_MERGE} ;






# step 1: join. crea 3 files: 1) join, 2-3) un1/un2
echo -e "[TIGET] Join overlapping paired-ends"
fastq-join $R1 $R2 -o ${OUTDIR}/${PATIENT}_${POOL}_%.fastq
gzip ${OUTDIR}/${PATIENT}_${POOL}_*.fastq

########################################
### BRANCH 1: JOINED FILE MANAGEMENT ###
########################################
echo -e "[TIGET] BRANCH 1: from joined reads with consensus, find A adapter. If it is on the left, ok; if right then revert and complement the sequence to be compliant to orientation of 454 reads (always 5'-3'). Merge results in a fastq file of reliable reads."
# step 2: identifica A: estrai LEFT/RIGHT_TAIL, rimuovi A [non mettere l'opzione -a ${ADAPTER454_A} perche nel caso di sovrapposizione totale potresti raddoppiarlo]
echo -e "[TIGET] Identify adapter 454 A, left and right"
flexbar --source <(zcat ${OUTDIR}/${PATIENT}_${POOL}_join.fastq.gz) --target ${TMPDIR}/${PATIENT}_${POOL}_join.fb.left -f fastq-i1.8 -as "GCCTCCCTCGCGCCATCAG" --threads 8 -ae LEFT_TAIL -at 1 -ao 18 -aa -m 30 -q 5 --removal-tag
## MEMENTO: ricorda che flexbar NON e' in grado di fare il reverse complement della tua sequenza di ingresso, quindi la cerca AS IS! per questo devi metterlo tu in RC: CTGATGGCGCGAGGGAGGC
flexbar --source <(zcat ${OUTDIR}/${PATIENT}_${POOL}_join.fastq.gz) --target ${TMPDIR}/${PATIENT}_${POOL}_join.fb.right -f fastq-i1.8 -as "CTGATGGCGCGAGGGAGGC" --threads 8 -ae RIGHT_TAIL -at 1 -ao 18 -aa -m 30 -q 5 --removal-tag

# filtra il fastq di input con i soli header giusti (qui left):
# come prima cosa prendi i soli header con A (flag flexbar), LEFT!
echo -e "[TIGET] Handle left reads"
grep "Flexbar" ${TMPDIR}/${PATIENT}_${POOL}_join.fb.left.fastq | cut -d' ' -f1,2 | sed 's/@//g' > ${TMPDIR}/${PATIENT}_${POOL}_ltr.left.list
# come seconda cosa fai filtro fastq
zcat ${OUTDIR}/${PATIENT}_${POOL}_join.fastq.gz | fqextract ${TMPDIR}/${PATIENT}_${POOL}_ltr.left.list > ${TMPDIR}/${PATIENT}_${POOL}_join.ltr.left.fastq
# a questo punto hai reads NON ripulite!!! quindi originali su cui poi devi fare tutto (trimming e mux)

# step 3: fai il reverse complement delle RIGHT_TAIL
# filtra il fastq di input con i soli header da sistemare (qui roght):
# come prima cosa prendi i soli header con A (flag flexbar), RIGHT
echo -e "[TIGET] Handle right reads"
grep "Flexbar" ${TMPDIR}/${PATIENT}_${POOL}_join.fb.right.fastq | cut -d' ' -f1,2 | sed 's/@//g' > ${TMPDIR}/${PATIENT}_${POOL}_ltr.right.list
# ora dall'header list RIGHT elimina le entry che sono gia' state prese come LEFT -> evita la RIPETIZIONI!
comm -13 <(cat ${TMPDIR}/${PATIENT}_${POOL}_ltr.left.list | sort) <(cat ${TMPDIR}/${PATIENT}_${POOL}_ltr.right.list | sort) >  ${TMPDIR}/${PATIENT}_${POOL}_ltr.right.clean.list
# come seconda cosa fai filtro fastq
zcat ${OUTDIR}/${PATIENT}_${POOL}_join.fastq.gz | fqextract ${TMPDIR}/${PATIENT}_${POOL}_ltr.right.clean.list > ${TMPDIR}/${PATIENT}_${POOL}_join.ltr.right2rc.fastq
# ora fai RC delle sequenze fastq filtrate
fastx_reverse_complement -i ${TMPDIR}/${PATIENT}_${POOL}_join.ltr.right2rc.fastq -o ${TMPDIR}/${PATIENT}_${POOL}_join.ltr.right.fastq -Q 33
# a questo punto hai reads NON ripulite!!! quindi originali su cui poi devi fare tutto (trimming e mux)
# clean tmp file
rm ${TMPDIR}/${PATIENT}_${POOL}_join.ltr.right2rc.fastq

# step 3.1: unisci i dati ripuliti
echo -e "[TIGET] Merge valid left and right"
cat ${TMPDIR}/${PATIENT}_${POOL}_join.ltr.left.fastq > ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.fastq
cat ${TMPDIR}/${PATIENT}_${POOL}_join.ltr.right.fastq >> ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.fastq
gzip ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.fastq

# step 4: rimuovi A e B
echo -e "[TIGET] Remove A and B"
fastq-mcf ${ADAPTER454_A} ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.fastq.gz -o ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.noA.fastq.gz -m 19 -p 1 -l 30 -q 10 -P 33 -S
fastq-mcf ${ADAPTER454_B} ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.noA.fastq.gz -o ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.noAB.fastq.gz -m 15 -p 20 -l 30 -q 10 -P 33 -S

# step 5: barcode split
echo -e "[TIGET] Split (mux) by barcodes"
fastq-multx -m 0 -B ${BARCODELIST} ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.noAB.fastq.gz -o ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.noAB_%.fastq.gz

# step 6: fastq to fasta
echo -e "[TIGET] Create FastA from FastQ for each TAG"
for TAG in ${BARCODEARRAY[@]}; do 
	zcat ${OUTDIR}/${PATIENT}_${POOL}_join.ltr.noAB_${TAG}.fastq.gz | fastq_to_fasta -Q 33 -o ${FASTADIR}/${PATIENT}_${POOL}_join_${TAG}.fasta ;
done

#############################################
### BRANCH 2: PAIRED-ENDS FILE MANAGEMENT ###
#############################################
## handle: un1.fastq.gz/un2.fastq.gz
echo -e "[TIGET] BRANCH 2: from paired ends reads that do NOT overlap themself, find A and LTR, filter these reads, and then proceed as usual (trimming A, B, and barcoding)."
# work on both on R1 and R2
## detect LTR as first step -> these labeled reads will be potential real IS
echo -e "[TIGET] Detect adapter 454 A"
flexbar --source <(zcat ${OUTDIR}/${PATIENT}_${POOL}_un1.fastq.gz) --source2 <(zcat ${OUTDIR}/${PATIENT}_${POOL}_un1.fastq.gz) --target ${TMPDIR}/${PATIENT}_${POOL}_single.fb.454A -f fastq-i1.8 -as "GCCTCCCTCGCGCCATCAG" --threads 8 -ae LEFT_TAIL -at 1 -ao 18 -aa -m 30 -q 10 --removal-tag
# get only reads headers
grep "Flexbar" ${TMPDIR}/${PATIENT}_${POOL}_single.fb.454A_1.fastq | cut -d' ' -f1,2 | sed 's/@//g' > ${TMPDIR}/${PATIENT}_${POOL}_single.fb.454A.only.list
grep "Flexbar" ${TMPDIR}/${PATIENT}_${POOL}_single.fb.454A_2.fastq | cut -d' ' -f1,2 | sed 's/@//g' >> ${TMPDIR}/${PATIENT}_${POOL}_single.fb.454A.only.list

## detect LTR as first step -> these labeled reads will be potential real IS
echo "[TIGET] Detect LTR reads (Flexbar)"
####### ATTENZIONE!!!! CAMBIARE IL PARAMETRO DI VALUTAZIONE DINAMICA (-aa) IMPATTA PESANTEMENTE SUL NUMERO FINALE DI RICONOSCIMENTI
####--threads 8 -ae ANY -at 3 -ao 58 -aa -m 20 -q 10 --removal-tag
####--threads 8 -ae ANY -at 3 -ao 58 -m 20 -q 10 --removal-tag
flexbar --source <(zcat ${OUTDIR}/${PATIENT}_${POOL}_un1.fastq.gz) --source2 <(zcat ${OUTDIR}/${PATIENT}_${POOL}_un1.fastq.gz) --target ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr -f fastq-i1.8 -a ${LTR} --threads 8 -ae ANY -at 0.5 -ao 58 -aa -m 20 -q 10 --removal-tag
# get only reads headers
grep "Flexbar" ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr_1.fastq | cut -d' ' -f1,2 | sed 's/@//g' > ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr.only.list
grep "Flexbar" ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr_2.fastq | cut -d' ' -f1,2 | sed 's/@//g' >> ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr.only.list
#rm ${TMPDIR}/${PATIENT}_${POOL}_${STEPNUM_CURRENT}_flexbar*

### keep only THE intersection! because only in this way you know exactly which pair has got the ltr
echo -e "[TIGET] Select only intersection reads between A-list and LTR-list"
comm -12 <(cat ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr.only.list | sort) <(cat ${TMPDIR}/${PATIENT}_${POOL}_single.fb.454A.only.list | sort) > ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr.list

# now we can merge R1 and R2 fastq and then filter only LTR reads using the final list
echo -e "[TIGET] Merge paired reads into a single file"
cat ${OUTDIR}/${PATIENT}_${POOL}_un1.fastq.gz ${OUTDIR}/${PATIENT}_${POOL}_un2.fastq.gz > ${OUTDIR}/${PATIENT}_${POOL}_single.fastq.gz
zcat ${OUTDIR}/${PATIENT}_${POOL}_single.fastq.gz | fqextract ${TMPDIR}/${PATIENT}_${POOL}_single.fb.ltr.list > ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.fastq
gzip ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.fastq

# step 4: rimuovi A e B
echo -e "[TIGET] Remove A and B adapters"
fastq-mcf ${ADAPTER454_A} ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.fastq.gz -o ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.noA.fastq.gz -m 19 -p 1 -l 30 -q 10 -P 33 -S
fastq-mcf ${ADAPTER454_B} ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.noA.fastq.gz -o ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.noAB.fastq.gz -m 15 -p 20 -l 30 -q 10 -P 33 -S

# step 5: barcode split
echo -e "[TIGET] Mux barcodes"
fastq-multx -m 0 -B ${BARCODELIST} ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.noAB.fastq.gz -o ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.noAB_%.fastq.gz

# step 6: fastq to fasta
echo -e "[TIGET] Convert FastQ into FastA"
for TAG in ${BARCODEARRAY[@]}; do 
	zcat ${OUTDIR}/${PATIENT}_${POOL}_single.ltr.noAB_${TAG}.fastq.gz | fastq_to_fasta -Q 33 -o ${FASTADIR}/${PATIENT}_${POOL}_single_${TAG}.fasta ;
done


#########################################
### UN-BRANCH: COMBINE FASTA BARCODES ###
#########################################
echo -e "[TIGET] FINAL OPERATIONS: create a unique FastA file for each TAG."
for TAG in ${BARCODEARRAY[@]}; do 
	cat ${FASTADIR}/${PATIENT}_${POOL}_single_${TAG}.fasta ${FASTADIR}/${PATIENT}_${POOL}_join_${TAG}.fasta > ${FASTADIR_MERGE}/${PATIENT}_${POOL}_${TAG}.fasta ;
done

echo -e "[TIGET] Done."

# clean folder
rm ${OUTDIR}/*.gz

