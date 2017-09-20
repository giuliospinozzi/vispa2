#!/bin/bash
######################  COLOR  ######################

	RED='\033[0;31m' 
	NC='\033[0m' # No Color
	GREEN='\033[0;32m'
	YELLOW='\033[1;33m'

######################  MY FUNCTION  ######################

	function downloadFromUCSC() {
		
		if [[ $1 == "hg38" ]]; then
			filename="hg38.fa.gz"
		elif [[ $1 == "hg19" || $1 == "mm10" || $1 == "mm9" || $1 == "mm8" ]] ; then
			filename="chromFa.tar.gz"
		else
			filename="chromFa.zip"
		fi

		finalLink="http://hgdownload.soe.ucsc.edu/goldenPath/$1/bigZips/$filename"

		echo $finalLink
	}

	function checkURL() {

		#check if the link exist
		status=`curl -s --head -w %{http_code} $1 -o /dev/null`
		
		if [[ $status == "200" ]]; then
			 echo "This URL Exist"
		else
		    echo "This URL Not Exist"
		    echo "Check the option -i used. Visit the link http://hgdownload.soe.ucsc.edu/goldenPath/ to search the right Keyword like hg19-mm10-mm9..."
		    exit -1
		fi
	}

######################  LINKS  ######################

	link_trimmomatic='http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip'
	link_vispa2='https://bitbucket.org/andreacalabria/vispa2'

######################  HELP  ######################

	usage()
	{
		echo
		echo "Set up VISPA-2 PE"
		echo
		echo "#Summary:"
		echo "	This program is an introduction for VISPA2"
		echo "	You will be able to set up your machine to allow work Vispa2 correctly"
		echo "	Will downloaded all the dependencies required; $0 installs packages to their own directory [/opt/applications/] and then symlinks their files into /usr/bin"
		echo "	Furthermore, can decide to download or not genomes present in the list below."
		echo "	If you prefer to configure manually your machine visit the link: https://bitbucket.org/andreacalabria/vispa2/wiki/VISPA2-PairedEnd"
		echo
		echo "#Usage: $0 [-s SPECIES] [-i GENOME.tar.gz]"
		echo
		echo "  [-s SPECIES] - Use the flag -x to indicate the species of the genome to download. If you use this flag, [-i is required]"
		echo "  [-i GENOME]  - Use the flag -i to download a Genome. If you use this flag, [-s is required]. 
			 Which GENOME would you like to Download? You MUST use keyword like:"
		echo
	  printf "	${GREEN}hg19${NC} --> human genome 19 from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz \n"
	  printf "	${GREEN}hg38${NC} --> human genome 38 from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \n"
	  printf "	${GREEN}mm9${NC} --> mouse 9 from http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz \n"
	  printf "	${GREEN}mm10${NC} --> mouse 10 from http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz \n" 
		echo
		exit 
	}

######################  MAIN  ######################

	# parse options
	while getopts ":s:i:h" Option
		do
		case $Option in
			s ) SPECIES="$OPTARG" ;;
			i ) GENOME="$OPTARG" ;;
			h ) usage ;;
			* ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;
		esac
	done
	shift $(($OPTIND - 1)) #shift $((OPTIND-1)) removes all the options that have been parsed by getopts from 
	#the parameters list, and so after that point, $1 will refer to the first non-option argument passed to the script.

	#check input
	if [[ $GENOME != "" ]]; then
		if [[ -z $SPECIES ]]; then
			echo "option '-i' requires option '-s'. use '-h' for usage information."; exit -1 ;
		fi
	fi

	if [[ -n $SPECIES ]]; then
		if [[ -z $GENOME ]]; then
			echo "option '-s' requires option '-i'. use '-h' for usage information."; exit -1 ;
		fi
	fi
	
	if [[ -n "$GENOME" ]]; then

		to_Download=`downloadFromUCSC $GENOME`
		checkURL $to_Download
		file_name=`echo $to_Download | cut -d'/' -f 7` #[ex chromFa.tar.gz]
		out_name=`echo $to_Download | cut -d'/' -f 5`  #[ex hg38]
	
	fi
	
######################  VARS and relative ORDER POSITION ######################
	
	# Create destination folder var
	APPLICATIONS="/opt/applications"
	# Create destination folder var
	BIN=${APPLICATIONS}"/bin"
	# Create destination folder var
	PATHGENOME="/opt/genome"
	# Create destination folder var
	if [ -n "$GENOME" ]; then
		FOLDERGENOME=${PATHGENOME}"/"$SPECIES"/"$GENOME 
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
	printf "${YELLOW}##### Installing apache2 #####${NC}\n"
		sudo apt install apache2
	echo ""
	printf "${YELLOW}##### php libapache2-mod-php php-mcrypt php-mysql #####${NC}\n"
		sudo apt-get install php libapache2-mod-php php-mcrypt php-mysql
	echo ""
	printf "${YELLOW}##### CLONING VISPA2 IN "$APPLICATIONS"/scripts${NC}\n"
	#cloning Vispa2
		cd $APPLICATIONS/scripts
		sudo hg clone $link_vispa2
	echo ""
	printf "${YELLOW}##### integration_analysis IN "$APPLICATIONS"/scripts${NC}\n"
	#cloning integration_analysis
		sudo hg clone https://bitbucket.org/tigetbioinformatics/integration_analysis

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
		sudo mkdir -p ${FOLDERGENOME}/annotation
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
		
		chr_numer=`ls | grep -E "[[:digit:]].fa" | wc -l`
		
		for i in $(seq 1 $chr_number); do
			cat chr$i.fa >> $out_name.fa
		done

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
		sudo mkdir -p ${FOLDERGENOME}/index/bwa_7

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
	sudo mkdir -p ${PATHGENOME}/vector
	printf "${YELLOW} @@@@ Folders creation --> "${PATHGENOME}"/control/phix${NC}\n"
	sudo mkdir -p ${PATHGENOME}/control/phix

###################### DATABASE CONFIGURATION ######################

	#VISPA2 uses MySQL as DBMS.
	echo "Install the dev lib and create three users (all privileges and readonly and all PRIVILEGES andrea)"
	echo "	
		GRANT ALL PRIVILEGES ON *.* To 'vispa2'@'localhost' IDENTIFIED BY 'vispa2';
		GRANT SELECT ON *.* TO 'readonly'@'localhost' IDENTIFIED BY 'readonlypswd';
		GRANT SELECT ON *.* TO 'andrea'@'localhost' IDENTIFIED BY 'andrea';
		MYSQL...
	"
	mysql -uroot -p -e "

		GRANT ALL PRIVILEGES ON *.* To 'vispa2'@'localhost' IDENTIFIED BY 'vispa2';

		GRANT SELECT ON *.* TO 'readonly'@'localhost' IDENTIFIED BY 'readonlypswd';

		GRANT ALL PRIVILEGES ON *.* TO 'andrea'@'localhost' IDENTIFIED BY 'andrea';

	"

###################### SOFTWARE CONFIGURATION ######################

	echo ""
	echo "STEP 4: SOFTWARE CONFIGURATION"
	echo ""

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

	sudo ln -s /opt/applications/scripts/vispa2/pipeline/illumina/VISPA2.IlluminaMiSeq.pipeline.sh /usr/bin/vispa2
	sudo ln -s /opt/applications/scripts/vispa2/script/annotate_matrix_v2.sh /usr/bin/annotate_matrix
	sudo ln -s /opt/applications/scripts/vispa2/script/annotate_bed.py /usr/bin/annotate_bed
	sudo ln -s /opt/applications/scripts/integration_analysis/src/Integration_Analysis_main.py /usr/bin/create_matrix


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
	sudo ln -s /opt/applications/scripts/vispa2/script/fastq_qf.sh /usr/bin/fastq_qf
	echo "link from /opt/applications/scripts/vispa2/script in /usr/bin/fasta_to_csv"
	sudo chmod +x /opt/applications/scripts/vispa2/script/fasta_to_csv.rb
	sudo ln -s /opt/applications/scripts/vispa2/script/fasta_to_csv.rb /usr/bin/fasta_to_csv

echo "Set up Completed!"