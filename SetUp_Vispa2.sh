#!/bin/bash

RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'


######################  LINKS  ######################

	#mother link
	link_ucsc='http://hgdownload.soe.ucsc.edu/goldenPath'
	#sons links
	link_hg19=${link_ucsc}'/hg19/bigZips/chromFa.tar.gz'
	link_hg38=${link_ucsc}'/hg38/bigZips/hg38.fa.gz'
	link_mm10=${link_ucsc}'/mm10/bigZips/chromFa.tar.gz'
	link_mm9=${link_ucsc}'/mm9/bigZips/chromFa.tar.gz'

	link_trimmomatic='http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip'
	link_vispa2='https://bitbucket.org/andreacalabria/vispa2'

usage()
{
	echo
	echo "Set up: VISPA 2 PE"
	echo
	echo "Usage: $0 [-s START] [-i GENOME.tar.gz]"
	echo
	echo "  [-s START] - Use the flag -s to initialize the program. Ex. $0 -s start"
	echo "  [-i GENOME]  - Use the flag -i to download a Genome. Which GENOME would you like to Download? You MUST use keyword like:"
  printf "	${GREEN}hg19${NC} --> human genome 19 from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz \n"
  printf "	${GREEN}hg38${NC} --> human genome 38 from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \n"
  printf "	${GREEN}mm9${NC} --> mouse 9 from http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz \n"
  printf "	${GREEN}mm10${NC} --> mouse 10 from http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz \n" 
	echo
	exit 
}

# parse options
while getopts ":si:h" Option
	do
	case $Option in
		s ) START="$OPTARG" ;;
		i ) GENOME="$OPTARG" ;;
		h ) usage ;;
		* ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;
	esac
done
shift $(($OPTIND - 1)) #shift $((OPTIND-1)) removes all the options that have been parsed by getopts from 
#the parameters list, and so after that point, $1 will refer to the first non-option argument passed to the script.

# check input GENOME

if [ "$GENOME" == "hg19" ]; then
	to_Download=$link_hg19
	species="human"
	folder_name="hg19"
elif [ "$GENOME" == "hg38" ]; then
	to_Download=$link_hg38
	species="human"
	folder_name="hg38"
elif [ "$GENOME" == "mm10" ]; then
	to_Download=$link_mm10
	folder_name="mm10"
	species="mouse"
elif [ "$GENOME" == "mm9" ]; then
	to_Download=$link_mm9
	folder_name="mm9"
	species="mouse"
fi

######################  VARS and relative ORDER POSITION ######################

	file_name=`echo $to_Download | cut -d'/' -f 7` #[ex chromFa.tar.gz]
	out_name=`echo $to_Download | cut -d'/' -f 5`	#[ex hg38]
	# Create destination folder var
	APPLICATIONS="/opt/applications"
	# Create destination folder var
	BIN=${APPLICATIONS}"/bin"
	# Create destination folder var
	PATHGENOME="/opt/genome"
	# Create destination folder var
	if [ -n "$GENOME" ]; then
		FOLDERGENOME=${PATHGENOME}"/"${species}"/"${folder_name}
	else 
		FOLDERGENOME=${PATHGENOME}"/species/"	
	fi

###################### DEPENDENCIES ######################


	echo "STEP 1: DEPENDENCIES"
	echo ""
	echo "Updating your system...Insert your password to update"
	sudo apt-get update
	#create folder if not exist
	sudo mkdir -p ${APPLICATIONS}/scripts
	cd ${APPLICATIONS}																											#pwd ----> /opt/applications/

	echo ""
	printf "${YELLOW}##### Installing MySQL #####${NC}\n"
		sudo apt-get install mysql-server
		#printf "${RED}REMEMBER THE PASSWORD FOR MYSQL, YOU WILL NEED LATER...${NC}\n" 
	echo ""
	printf "${YELLOW}##### Installing pigz #####${NC}\n"
		sudo apt-get install pigz
	echo ""
	printf "${YELLOW}##### Installing parallel #####${NC}\n"
		sudo apt-get install parallel
	echo ""
	printf "${YELLOW}##### Installing fastqc #####${NC}\n"
		sudo apt-get install fastqc
	echo ""
	printf "${YELLOW}##### Installing bwa #####${NC}\n"
		sudo apt-get install bwa
		sudo ln -s /usr/bin/bwa /usr/bin/bwa-stable
	echo ""
	printf "${YELLOW}##### Installing samtools #####${NC}\n"
		sudo apt-get install samtools
	echo ""
	printf "${YELLOW}##### Installing trimmomatic #####${NC}\n"
		curl -O $link_trimmomatic
		unzip Trimmomatic-0.36.zip
		cd Trimmomatic-0.36																										#pwd ----> /opt/applications/trimmomatic-0.36
		mkdir -p ${BIN}
		echo '#!/bin/bash' > ${BIN}/trimmomatic
		echo 'java -jar '${APPLICATIONS}'/Trimmomatic-0.36/trimmomatic-0.36.jar $@' >> ${BIN}/trimmomatic
		chmod +x ${BIN}/trimmomatic
	echo ""
	printf "${YELLOW}##### Installing ea-utils #####${NC}\n"
		sudo apt-get install ea-utils
	echo ""
	printf "${YELLOW}##### Installing flexbar #####${NC}\n"
		sudo apt-get install flexbar
	echo ""
	printf "${YELLOW}##### Installing bamtools #####${NC}\n"
		sudo apt-get install bamtools
	echo ""
	printf "${YELLOW}##### Installing picard-tools #####${NC}\n"
		sudo apt-get install picard-tools
		echo '#!/bin/bash' > ${BIN}/FilterSamReads
		echo 'picard-tools FilterSamReads $@' >> ${BIN}/FilterSamReads
		chmod +x ${BIN}/FilterSamReads
		echo '#!/bin/bash' > ${BIN}/MergeSamFiles
		echo 'picard-tools MergeSamFiles $@' >> ${BIN}/MergeSamFiles
		chmod +x ${BIN}/MergeSamFiles
	echo ""
	printf "${YELLOW}##### Installing bedtools #####${NC}\n"
		sudo apt-get install bedtools
	echo ""
	printf "${YELLOW}##### Installing curl #####${NC}\n"
		sudo apt-get install curl
	echo ""
	printf "${YELLOW}##### Installing mercurial #####${NC}\n"
		sudo apt-get install mercurial
	echo ""
	printf "${YELLOW}##### Installing pip #####${NC}\n"
		sudo apt-get install python-pip
		sudo pip install --upgrade pip
	echo ""
	printf "${YELLOW}##### Installing libmysqlclient-dev #####${NC}\n"
		sudo apt-get install libmysqlclient-dev
	echo ""
	printf "${YELLOW}##### CLONING VISPA2 IN "$APPLICATIONS"${NC}\n"
	#cloning Vispa2
		cd $APPLICATIONS/scripts
		sudo hg clone $link_vispa2

###################### CONFIGURATION ######################

	echo ""
	echo "STEP 2: CONFIGURATION"
	echo ""

	printf "${YELLOW}@@@@ Folder creation --> ${BIN}${NC}\n"
		sudo mkdir -p ${BIN}
	printf "${YELLOW}@@@@ Adjust permission --> ${APPLICATIONS}${NC}\n"
		sudo chmod -R 777 ${APPLICATIONS}
	printf "${YELLOW}@@@@ Folder creation --> ${PATHGENOME}${NC}\n"
		sudo mkdir -p ${PATHGENOME}
	printf "${YELLOW}@@@@ Adjust permission --> ${PATHGENOME}${NC}\n" 
		sudo chmod -R 777 ${PATHGENOME}
	printf "${YELLOW}@@@@ Folder creation --> /opt/NGS/results${NC}\n"
		sudo mkdir -p /opt/NGS/results
	printf "${YELLOW}@@@@ Adjust permission --> /opt/NGS${NC}\n"
		sudo chmod -R 777 /opt/NGS

###################### Download Genome ######################

	echo ""
	echo "STEP 3: Download Genome"
	echo ""

	
	if [ -z $to_Download ]; then
		echo "No Genome will be downloaded, SKIP to step 4..."
	else
		printf "Folders creation --> "${FOLDERGENOME}"/annotation"
		mkdir -p ${FOLDERGENOME}/annotation
		echo ""
		echo "Downloading from UCSC the FASTA file and annotation file in GTF and saved in the this path --> "${FOLDERGENOME}
		echo ""
		echo "Getting fileGenome FROM "$to_Download
		cd ${FOLDERGENOME}																													 	#pwd ----> /opt/genome/species/type
		curl -OL $to_Download
		echo ""
		echo "File Extraction of "$file_name
		tar -xvzf $file_name
		echo ""
		echo "File Concatenation..."
		
		if [ $species == "human" ]; then	
			for i in {1..22}; do
				cat chr$i.fa >> $out_name.fa
			done

		elif [ $species == "mouse" ]; then
			for i in {1..19}; do
				cat chr$i.fa >> $out_name.fa
			done

		fi

		cat chrX.fa >> $out_name.fa
		cat chrY.fa >> $out_name.fa
		cat chrM.fa >> $out_name.fa
		
		echo "...."

		printf "${GREEN}@@@@ Concatenation file complete!${NC}\n"
		echo ""
		echo "Removing extracted chr files...."
		mv $file_name ../
		rm chr*
		mv ../$file_name ${FOLDERGENOME}
		echo ""
		cd ${FOLDERGENOME}/annotation																											#pwd ----> /opt/genome/species/type/annotation
		
		if [[ $out_name == "hg19" ]]; then
			cp /opt/applications/scripts/vispa2/annotation/ucsc.$out_name.refSeq.gtf.tar.gz .
			tar -xvzf ${FOLDERGENOME}/annotation/ucsc.$out_name.refSeq.gtf.tar.gz
		else
			cp /opt/applications/scripts/vispa2/annotation/ucsc.$out_name.refSeq.gtf.gz .
			gunzip ${FOLDERGENOME}/annotation/ucsc.$out_name.refSeq.gtf.gz
		fi

		printf "${YELLOW}@@@@ Folders creation --> "${FOLDERGENOME}"/index/bwa_7${NC}\n"
		mkdir -p ${FOLDERGENOME}/index/bwa_7

		cd ${FOLDERGENOME}																														#pwd ----> /opt/genome/species/type/
		echo "moving "$out_name.fa" in "${FOLDERGENOME}"/index/bwa_7"
		mv $out_name.fa index/bwa_7
		echo ""
		printf "${NC}Indexing of $out_name.fa...this operation require a lot of time. Better if you go to take some coffee.\n"
		cd ${FOLDERGENOME}/index/bwa_7																											#pwd ----> /opt/genome/species/type/index/bwa_7

		bwa-stable index -a bwtsw ${FOLDERGENOME}/index/bwa_7/$out_name.fa  #BWA
		samtools faidx $out_name.fa 										#SAMTOOLS
		picard-tools CreateSequenceDictionary R=$out_name.fa O=$out_name.dict 	#PICARD

		printf "${GREEN}Indexing complete!${NC}\n"
	fi
	
	printf "${YELLOW} @@@@ Folders creation --> "${PATHGENOME}"/vector${NC}\n"
	mkdir -p ${PATHGENOME}/vector
	printf "${YELLOW} @@@@ Folders creation --> "${PATHGENOME}"/control/phix${NC}\n"
	mkdir -p ${PATHGENOME}/control/phix

###################### DATABASE CONFIGURATION ######################

	#VISPA2 uses MySQL as DBMS.
	echo "Install the dev lib and create two users (all privileges and readonly)"
	echo "	
		GRANT ALL PRIVILEGES ON *.* To 'vispa2'@'localhost' IDENTIFIED BY 'vispa2';
		GRANT SELECT ON *.* TO 'readonly'@'localhost' IDENTIFIED BY 'readonlypswd';
		MYSQL...
	"
	mysql -uroot -p -e "

		GRANT ALL PRIVILEGES ON *.* To 'vispa2'@'localhost' IDENTIFIED BY 'vispa2';

		GRANT SELECT ON *.* TO 'readonly'@'localhost' IDENTIFIED BY 'readonlypswd';

	"

###################### SOFTWARE CONFIGURATION ######################

	echo ""
	echo "STEP 4: SOFTWARE CONFIGURATION"
	echo ""

	echo "Soft link creation..."
	sudo ln -s /usr/bin/bwa /usr/bin/bwa-stable

	echo 'export LD_LIBRARY_PATH=/path/FlexbarDir:$LD_LIBRARY_PATH' >> ~/.bashrc
	source ~/.bashrc
	
	sudo ln -s /usr/bin/flexbar /usr/bin/flexbar2.5
	sudo ln -s ${BIN}/trimmomatic /usr/bin/trimmomatic
	sudo ln -s ${BIN}/FilterSamReads /usr/bin/FilterSamReads
	sudo ln -s ${BIN}/MergeSamReads /usr/bin/MergeSamReads

	# Python
	# The following packages are required:

	echo ""
	echo "STEP 5: PYTHON CONFIGURATION"
	echo ""

	sudo -H pip install MySQL-python
	sudo -H pip install pysam==0.7.7
	sudo -H pip install biopython
	sudo -H pip install HTSeq
	sudo -H pip install rpy2
	sudo -H pip install scipy
	sudo -H pip install numpy
	sudo -H pip install matplotlib
	sudo -H pip install xlsxwriter
	sudo -H pip install pandas

	cd $APPLICATIONS 																													#pwd ----> /opt/applications
	
	echo "Python program link creation in /usr/bin"
	echo "" 

	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/import_iss"
	sudo ln -s /opt/applications/scripts/vispa2/script/import_iss.py /usr/bin/import_iss
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/fqreverseextract.pureheader"
	sudo ln -s /opt/applications/scripts/vispa2/script/fqreverseextract.pureheader.py /usr/bin/fqreverseextract.pureheader
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/fqextract.pureheader"
	sudo ln -s /opt/applications/scripts/vispa2/script/fqextract.pureheader.py /usr/bin/fqextract.pureheader
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/rev_extract_header"
	sudo ln -s /opt/applications/scripts/vispa2/script/rev_extract_header.py /usr/bin/rev_extract_header
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/extract_header"
	sudo ln -s /opt/applications/scripts/vispa2/script/extract_header.py /usr/bin/extract_header
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/filter_by_cigar_bam"
	sudo ln -s /opt/applications/scripts/vispa2/script/filter_by_cigar_bam.py /usr/bin/filter_by_cigar_bam
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/filter_by_mate"
	sudo ln -s /opt/applications/scripts/vispa2/script/filter_by_mate.py /usr/bin/filter_by_mate
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/isa_importrediss_frombed"
	sudo ln -s /opt/applications/scripts/vispa2/script/dbimport_redundantiss_from_bed.v2.py /usr/bin/isa_importrediss_frombed
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/fastq_qf"
	sudo chmod +x /opt/applications/scripts/vispa2/script/fastq_qf.sh
	sudo ln –s /opt/applications/scripts/vispa2/script/fastq_qf.sh /usr/bin/fastq_qf
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/fasta_to_csv"
	sudo chmod +x /opt/applications/scripts/vispa2/script/fasta_to_csv.rb
	sudo ln –s /opt/applications/scripts/vispa2/script/fasta_to_csv.rb /usr/bin/fasta_to_csv


echo "Set up Completed!"
