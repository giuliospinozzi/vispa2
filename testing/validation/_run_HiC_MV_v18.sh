#!/bin/bash
source /etc/environment
source /etc/profile

#   +--------------------------------------------------------+
#   |                                                        |
#   |       Illumina NEW Pipeline: data -> Merged BAM        |
#   |       - no 454 adapters, illumina-based only           |
#   |                                                        |
#   +--------------------------------------------------------+
#   |  Author:   Andrea Calabria                             |
#   |  Date:     January 2013                                |
#   |  Contact:  andrea.calabria@hsr.it                      |
#   |  Version:  1.3.1.10                                    |
#   +--------------------------------------------------------+
# 
#   REQUIRED VARS and relative ORDER POSITION -> REMEMBER NO SPACES!!!!!!!!!
# 	1. disease id (MLD)
# 	2. patient id (MLD02)
# 	3. local server working path [/opt/NGS/results] -> it will be added DISEASE and PATIENT ID [/opt/NGS/results/MLD/MLD02/]
# 	4. R1
# 	5. R2
# 	6. pool id
# 	7. barcode list LTR
# 	8. barcode list LC
# 	9. reference genome (NB: indexed by BWA, hg19:: /opt/NGS/index/human/hg19/bwa/hg19.nochr0)
# 	10. tmp dir
# 	11. association file
# 	12. db host id
# 	13. db target schema
# 	14. db target table
# 	15. PhiX genome (bwa index, i.e.: PHIX=/opt/NGS/index/phix/phix174/bwa/phiX174.fasta)
# 	16. Vector genome (i.e.: LV 751: LVGENOME=/opt/NGS/index/virus/lv751wg/bwa/lv.fa)
# 	17. Contaminant FastA DB (i.e: /opt/NGS/vectors/UniVec_Tiget.fa)
# 	18. Number of MAXIMUM processes to apply to each step (max=number of CPU)
# 	19. GATK reference genome (i.e. /opt/NGS/index/human/hg19/bwa/hg19.nochr0.fa)
# 	20. Genome ID for CIGAR-MD program. NB: this option collects only whole genome so far, not single regions (even the program allows it!). Example: hg19, mm9.
# 	21. Stringent regions BED file (for example: /opt/NGS/annotation/hg19/ucsc.simplerepeat.hg19.bed or /opt/NGS/annotation/mm9/ucsc.simplerepeat.mm9.bed)
# 	

#########################################################################################################################################################

######################################################
#                                                    #
#                    MV HiC                          #
#                                                    #
######################################################

##### POOL 1 ######

TODAY=`date +"%Y%m%d%H%M"`;
R1="/tmp/mv1/SA55_S1_L001_R1_001.fastq.gz" ;
R2="/tmp/mv1/SA55_S1_L001_R2_001.fastq.gz" ;
POOLNAME="MV1_p18lc";
PATIENT="General";
GENOME="/opt/genome/mouse/mm9/index/bwa_7/mm9.fa"; # the right one
DISEASE="Monica";
BARCODE_LTR="/opt/applications/scripts/isatk/elements/barcode/barcode.hic.mv.ltr.list" ;
BARCODE_LC="/opt/applications/scripts/isatk/elements/barcode/barcode.hic.mv.lc.list" ;
NGSWORKINGPATH="/tmp/mv1/results"; # WITHOUT LAST SLASH -> if not present, I will create it a new one with this name using mkdir
TMPDIR="/tmp/${TODAY}" ;
LOGF="/tmp/mv1/${TODAY}.${DISEASE}.${PATIENT}.${POOLNAME}.log"
mkdir ${TMPDIR};
ASSOCIATIONFILE="/opt/applications/scripts/isatk/elements/association/asso.hic.mv.tsv"
DBHOSTID="local"
DBTARGETSCHEMA="sequence_hic"
DBTARGETTABLE="mv1_p18lc"
PHIXGENOME="/opt/genome/control/phix174/bwa_7/phiX174.fa"
LVGENOME="/opt/genome/vector/lv751wg/bwa_7/lv.fa"
CONTAMINANTDB="/opt/applications/scripts/isatk/elements/sequences/UniVec.2013.TIGET.LV.ARSA.fa";
GATKREFGENOME="/opt/genome/mouse/mm9/index/bwa_7/mm9.fa"
CIGARGENOMEID="mm9" ;
# available CPUs
CPUN="`cat /proc/cpuinfo | grep "model name" | wc -l`";
MAXTHREADS=$((CPUN-2))
MAXTHREADS="8"
STRINGETREGIONSBED="/opt/NGS/annotation/mm9/ucsc.simplerepeat.mm9.bed"


##########

LEFTSPACE=`df | grep /dev/md0 | cut -d' ' -f18`
if [ ${LEFTSPACE} -lt 15000000 ]; then 
	echo "
	
    *****************************************************************
    |                                                               |
    |   YOU DO NOT HAVE ENOUGH SPACE LEFT ON HARD DISC!!            |
    |                                                               |
    |   I WILL NOT PROCEED...                                       |
    |                                                               |		
    |   FREE SPACE BEFORE, LEAVING AT LEAST 10GB ON MAIN DISC       |
    |                                                               |
    *****************************************************************
	
		"; 
else
	echo "
	STARTING PROCESSING AT:
	"  > ${LOGF};
	date >> ${LOGF}
	 /opt/applications/scripts/isatk/pipeline/illumina/./NewIlluminaMiSeq.pipeline.v1.3.1.18_gemini.sh ${DISEASE} ${PATIENT} ${NGSWORKINGPATH} ${R1} ${R2} ${POOLNAME} ${BARCODE_LTR} ${BARCODE_LC} ${GENOME} ${TMPDIR} ${ASSOCIATIONFILE} ${DBHOSTID} ${DBTARGETSCHEMA} ${DBTARGETTABLE} ${PHIXGENOME} ${LVGENOME} ${CONTAMINANTDB} ${MAXTHREADS} ${GATKREFGENOME} ${CIGARGENOMEID} ${STRINGETREGIONSBED} 2>&1 >> ${LOGF} ;
	echo "
	FINISHED PROCESSING AT:
	"  >> ${LOGF};
	date >> ${LOGF};
	# gzip ${LOGF}
fi
#rm -fr ${TMPDIR};

#########################################################################################################################################################
