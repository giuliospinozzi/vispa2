#!/bin/bash
source /etc/environment
source /etc/profile

PROJECT="VISPA2-test-CEM";
VISPA2_HOME="/opt/applications/scripts/vispa2"

echo "
        +--------------------------------------------------------+
        |                                                        |
        |                 VISPA2 (Sample Dataset)                |
        |                                                        |
        +--------------------------------------------------------+
        |  Author:   Giulio Spinozzi, PostDoc                    |
        |  Date:     Sept 2017                              	 |
        |  Version:  1.0                                         |  
        |  Contact:  spinozzi.giulio@hsr.it                      |
        +--------------------------------------------------------+

        >PROJECT: ${PROJECT}

"


TODAY=`date +"%Y%m%d%H%M%S"`;


##### ========================= Parameters to Check ========================= #####

R1="${VISPA2_HOME}/testing/sampledata/AssayValidation_R1.fastq.gz";
R2="${VISPA2_HOME}/testing/sampledata/AssayValidation_R2.fastq.gz";

# the file system tree folder will be as follows: DISEASE -> PATIENT -> {bam, bcmultx, iss, quality, bed, quantification} -> POOLNAME
DISEASE="VISPA2"; # project main name: for example the disease or the reference project ID
PATIENT="sampledata"; # place here a patient ID, or an experiment ID
POOLNAME="testcase"; # the sequencing pool ID 

GENOME="/opt/genome/human/hg19/index/bwa_7/hg19.fa"; ## hg19: /opt/genome/human/hg19/index/bwa_7/hg19.fa ; mm9: /opt/genome/mouse/mm9/index/bwa_7/mm9.fa ; mfa5: /opt/genome/monkey/mfa5/index/bwa_7/mfa5.fa
BARCODE_LTR="${VISPA2_HOME}/elements/barcode/barcode.LTR.48.list";
BARCODE_LC="${VISPA2_HOME}/elements/barcode/barcode.LC.48.list";

ASSOCIATIONFILE="${VISPA2_HOME}/testing/sampledata/AssayValidation.tsv";  # https://www.dropbox.com/s/6slo4l8g9i8sk3s/association%20file%20template.xlsx
LTR="${VISPA2_HOME}/elements/sequences/LTR.32bp.fa"; # LTR in forward
LTR_rc="${VISPA2_HOME}/elements/sequences/LTR.32bp.rev.fa"; # LTR in reverse complement
LC_fwd="${VISPA2_HOME}/elements/sequences/LC.assayvalidation.fwd.fa"; # Linker Cassette in forward
LC_rev="${VISPA2_HOME}/elements/sequences/LC.assayvalidation.rc.fa"; # Linker Cassette in reverse

DBHOSTID="local"; ## lkeywords, to keed as it is for localhost. Otherwise please change also the program import_iss_from_bed
DBTARGETSCHEMA="review_VISPA2";
DBTARGETTABLE="testcase";

REPEATS="repeats_no";
GATKREFGENOME="/opt/genome/human/hg19/index/bwa_7/hg19.fa"; # riportare lo stesso file del genoma GENOME:: hg19: /opt/genome/human/hg19/index/bwa_7/hg19.fa ; mm9: /opt/genome/mouse/mm9/index/bwa_7/mm9.fa ; mfa5: /opt/genome/monkey/mfa5/index/bwa_7/mfa5.fa hiv: /opt/genome/hiv/hiv_hxb2cg/bwa_7/hiv.fa
CIGARGENOMEID="hg19" ; # Reference genome ID: choose among {hg19 | mm9 | mfa5}
VECTORCIGARGENOMEID="lv"; ## This is the vector reference name (id) used to remove vector sequences. Choose among: {lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada, hiv}
LVGENOME="/opt/genome/vector/lv/bwa_7/lv.backbone.fa"; # Change it ONLY if you want to quantify and remove other vectors or inserted sequences. Alternatives in the GEMINI folder /opt/genome/vector/lv/bwa_7/: {lv.backbone.fa, lv.backbone.hpgk.arsa.wprem.fa, lv.backbone.wasp.was.wprem.fa, lv.plasmid.amp.fa, lv.plasmid.kana.fa}. HIV: /opt/genome/hiv/hiv_hxb2cg/bwa_7/hiv.fa
PHIXGENOME="/opt/genome/control/phix174/bwa_7/phiX174.fa";
CONTAMINANTDB="${VISPA2_HOME}/elements/sequences/UniVec_Tiget.fa"; # con lo stesso path, le alternative sono: {UniVec_Tiget_Gamma.fa, UniVec.2013.TIGET.LV.ARSA.fa, UniVec.2013.TIGET.LV.WAS.fa}

## Computational Parameters
CPUN="`cat /proc/cpuinfo | grep "model name" | wc -l`";
MAXTHREADS=6;
SUBOPTIMALTHRESHOLD='40';
FASTQ_QF="lam"; # FASTQ Quality Filter Methods: slim (QF on R1 80bp and R2 TAGs) or lam (QF only on R1 80bp)

NGSWORKINGPATH="/opt/NGS/results"; # WITHOUT LAST SLASH -> if not present, I will create it a new one with this name using mkdir
REMOVE_TMP_DIR="remove_tmp_yes"; # remove tmp dirs? remove_tmp_yes

#=================================================================================#


##### ========================== Fixed Parameters =========================== #####
## Creating Dirs
TMPDIR="/opt/NGS/pipetmpdir/${TODAY}";
mkdir -p ${TMPDIR};
#=================================================================================#


#---------------------------------------***---------------------------------------#
##### =============================== PROGRAM =============================== #####
## Space
LEFTSPACE=`df -P ${TMPDIR} | tail -1 | awk '{print $4}'`
if [ ! -r "$R1" ]; then
	echo "Error: can't open input file for R1."
	exit 1
fi
if [ ! -r "$R2" ]; then
	echo "Error: can't open input file for R2."
	exit 1
fi
if [ ! -r "$ASSOCIATIONFILE" ]; then
	echo "Error: can't open input ASSOCIATION FILE."
	exit 1
fi
echo "
STARTING PROCESSING AT:
"
date
vispa2 ${DISEASE} ${PATIENT} ${NGSWORKINGPATH} ${R1} ${R2} ${POOLNAME} ${BARCODE_LTR} ${BARCODE_LC} ${GENOME} ${TMPDIR} ${ASSOCIATIONFILE} ${DBHOSTID} ${DBTARGETSCHEMA} ${DBTARGETTABLE} ${PHIXGENOME} ${LVGENOME} ${CONTAMINANTDB} ${MAXTHREADS} ${GATKREFGENOME} ${CIGARGENOMEID} ${VECTORCIGARGENOMEID} ${SUBOPTIMALTHRESHOLD} ${REMOVE_TMP_DIR} ${LTR} ${LTR_rc} ${LC_fwd} ${LC_rev} ${FASTQ_QF} ${REPEATS}
echo "
FINISHED PROCESSING AT:
"
date;
#=================================================================================#
#---------------------------------------***---------------------------------------#


# echo "
#         +--------------------------------------------------------+
#         |                                                        |
#         |                  VISPA2 (Create Matrix)                |
#         |                                                        |
#         +--------------------------------------------------------+
#         |  Author:   Giulio Spinozzi, PostDoc                    |
#         |  Date:     Sept 2017                                   |
#         |  Version:  1.0                                         |  
#         |  Contact:  spinozzi.giulio@hsr.it                      |
#         +--------------------------------------------------------+

#         >PROJECT: ${PROJECT}

# "


TODAY=`date +"%Y%m%d%H%M%S"`;


create_matrix --dbDataset "${DBTARGETSCHEMA}.${DBTARGETTABLE}" --columns tissue,sample,treatment,vector,enzyme --IS_method classic --bp_rule 7 --tsv --no_xlsx

# echo "
#         +--------------------------------------------------------+
#         |                                                        |
#         |                VISPA2 (Annotate Matrix)                |
#         |                                                        |
#         +--------------------------------------------------------+
#         |  Author:   Giulio Spinozzi, PostDoc                    |
#         |  Date:     Sept 2017                                   |
#         |  Version:  1.0                                         |  
#         |  Contact:  spinozzi.giulio@hsr.it                      |
#         +--------------------------------------------------------+

#         >PROJECT: ${PROJECT}

# "


TODAY=`date +"%Y%m%d%H%M%S"`;

annotate_matrix -m IS_matrix_classic_strand_specific_method_review_VISPA2_testcase.tsv -g /opt/applications/scripts/vispa2/annotation/ucsc.hg19.refSeq.gtf -o . -t vispa

