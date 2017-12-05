#!/bin/bash
######################  COLOR  ######################

	RED='\033[0;31m' 
	NC='\033[0m' # No Color
	GREEN='\033[0;32m'
	YELLOW='\033[1;33m'
	CYAN='\033[0;36m'
	BLUE='\033[0;34m'
	PURPLE='\033[0;35m'

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

		if [[ `dpkg -s curl 2> /dev/null | wc -l` > 0 ]]; then
			status=`curl -s --head -w %{http_code} $1 -o /dev/null`
			#check if the link exist
			if [[ $status == "200" ]]; then
				echo "This URL Exist"
			else
				echo "This URL Not Exist"
				echo "Check the option -i used. Visit the link http://hgdownload.soe.ucsc.edu/goldenPath/ to search the right Keyword like hg19-mm10-mm9..."
				exit -1
			fi
		else
			printf "[ ${RED}CURL NOT FOUND${NC} ] curl not installed\n"
			echo Download...
			sudo apt-get install curl
			checkURL $to_Download
		fi
	}

	function installDepend () {
		sudo apt-get install $1
	}

	function installPyPack () {
		sudo -H pip install $1
	}

	function len () {
        echo $(wc -w <<< "$@")
	}

	function checkAllPacks () {

		echo "@  Checking if all pack are installed"
		list_dependencies=(mysql-server pigz parallel curl bwa samtools ea-utils flexbar fastx-toolkit bamtools picard-tools bedtools mercurial python-pip libmysqlclient-dev apache2 ruby php libapache2-mod-php php-mcrypt php-mysql dos2unix)
		for i in ${list_dependencies[@]}; do
			if [[ `dpkg -s $i 2> /dev/null | wc -l` > 0 ]]; then
				printf "[ ${GREEN}OK${NC} ] $i installed correctly\n"
			else
				printf "[ ${RED}MISSING${NC} ] $i not installed correctly\n"
				installDepend "$i"
			fi
		done
	}

######################  HELP  ######################

	usage()
	{
		echo
		echo "Setting up VISPA2"
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

		if [ ${GENOME:0:2} = "hg" ]; then
			SPECIES="human";
		else
			if [ ${GENOME:0:2} = "mm" ]; then
			SPECIES="mouse";
			fi
		fi
		
		file_name=`echo $to_Download | cut -d'/' -f 7` #[ex chromFa.tar.gz]
		out_name=`echo $to_Download | cut -d'/' -f 5`  #[ex hg38]
	
	fi
	
######################  VARS and relative ORDER POSITION ######################
	
	VISPA2="/opt/applications/scripts/vispa2"
	INTEGRATION="/opt/applications/scripts/integration_analysis"
	SCRIPT="${VISPA2}/script"
	APPLICATIONS="/opt/applications"
	BIN=${APPLICATIONS}"/bin"
	PATHGENOME="/opt/genome"
	
	if [ -n "$GENOME" ]; then
		FOLDERGENOME=${PATHGENOME}"/"$SPECIES"/"$GENOME 
	fi

######################  LINKS  ######################
	
	__fastqc="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip"
	__trimmomatic='http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip'
	__vispa2='https://bitbucket.org/andreacalabria/vispa2'
	__version_trim="Trimmomatic-0.36.zip"
	__version_fastqc="fastqc_v0.11.5.zip"
	__integration_analysis="https://bitbucket.org/tigetbioinformatics/integration_analysis"

###################### DEPENDENCIES ######################

	list_dependencies=(mysql-server fastqc trimmomatic pigz parallel curl flexbar fastx-toolkit bamtools picard-tools bedtools mercurial python-pip libmysqlclient-dev apache2 ruby php libapache2-mod-php php-mcrypt php-mysql dos2unix $__vispa2 $__integration_analysis bwa samtools ea-utils)

	echo "
        +--------------------------------------------------------+
        |                                                        |
        |                 CONFIGURATION (Vispa2)                 |
        |                                                        |
        +--------------------------------------------------------+
        |  Collaborator:  Adriano De Marino, PhD Student I year  |
        |  Date:     Sept 2017                                   |
        |  Version:  1.2                                         |  
        |  Contact:  demarino.adriano@hsr.it (HSR-TIGET)         |
        +--------------------------------------------------------+

    "
	echo "STEP 1: DEPENDENCIES"
	echo "Installing missing dependencies"
	echo
	
	echo "Updating your system...Insert your password to update"
	sudo apt-get update
	sudo apt-get upgrade
	#create folder if not exist
	sudo mkdir -p ${APPLICATIONS}/scripts
	sudo chmod -R 777 ${APPLICATIONS}
	cd ${APPLICATIONS}																		#pwd ----> /opt/applications/

	for i in "${list_dependencies[@]}"; do 
		
		cd ${APPLICATIONS}/scripts 

		if [[ $i == "fastqc" ]]; then

			cd ${APPLICATIONS}
			if [[( -a /etc/fastqc/fastqc && -a /usr/bin/fastqc )]]; then
				printf "[ ${GREEN}OK${NC} ] Fastqc already exists\n"
				sudo chmod +x /etc/fastqc/fastqc > /dev/null
				continue
			else
				printf "[ ${RED}NO${NC} ] Fastqc not found, Installing ${i}\n"
				sudo apt-get install curl 2> /dev/null
				sudo curl -OL $__fastqc 
				unzip fastqc_v0.11.5.zip > /dev/null
				sudo mv FastQC fastqc 2> /dev/null
				sudo mv fastqc /etc 2> /dev/null
				sudo chmod +x /etc/fastqc/fastqc
				sudo ln -sf /etc/fastqc/fastqc /usr/bin/fastqc 2> /dev/null
			fi	
		elif [[ $i == "trimmomatic" ]]; then

			cd ${APPLICATIONS}
			if [[( -a ${APPLICATIONS}/bin/trimmomatic && -a /usr/bin/trimmomatic )]]; then
				printf "[ ${GREEN}OK${NC} ] Trimmomatic already exists\n"
				continue
			else
				printf "[ ${RED}NO${NC} ] Trimmomatic not found, Installing ${i}\n"
				sudo apt-get install curl
				sudo curl -OL $__trimmomatic > /dev/null
				unzip Trimmomatic-0.36.zip > /dev/null
				cd Trimmomatic-0.36																										#pwd ----> /opt/applications/trimmomatic-0.36
				mkdir -p ${BIN}
				echo "#!/bin/bash" > ${BIN}/trimmomatic
				echo "java -jar '${APPLICATIONS}'/Trimmomatic-0.36/trimmomatic-0.36.jar \$@" >> ${BIN}/trimmomatic
				chmod +x ${BIN}/trimmomatic
			fi
		elif [[ $i == $__vispa2 ]]; then
			if [[ `dpkg -s mercurial 2> /dev/null | wc -l ` > 0 ]]; then
				printf "[ ${GREEN}OK${NC} ] Mercurial already exists\n"
				printf "[ ${YELLOW}CLONING${NC} ] repository vispa2 in ${VISPA2}\n"
				hg clone $__vispa2 > /dev/null
			else
				printf "[ ${RED}NO${NC} ] Mercurial not found, Installing Mercurial to download repository\n"			
				installDepend "mercurial" > /dev/null
				printf "[ ${YELLOW}CLONING${NC} ] repository vispa2 in ${VISPA2}\n"
				hg clone $__vispa2 > /dev/null
			fi
		elif [[ $i == $__integration_analysis ]]; then
			if [[ `dpkg -s mercurial 2> /dev/null | wc -l ` > 0 ]]; then
				printf "[ ${GREEN}OK${NC} ] Mercurial already exists\n"
				printf "[ ${YELLOW}CLONING${NC} ] repository integration_analysis in ${INTEGRATION}\n"
				hg clone -b "2.1-seqTracker" $__integration_analysis 2> /dev/null
			else
				printf "[ ${RED}NO${NC} ] Mercurial not found, Installing Mercurial to download repository\n"
				printf "[ ${YELLOW}CLONING${NC} ] repository integration_analysis in ${INTEGRATION}\n"			
				installDepend "mercurial" > /dev/null
				hg clone -b "2.1-seqTracker" $__integration_analysis 2> /dev/null
			fi
		elif [[ $i == "bedtools" ]]; then
			cd ${APPLICATIONS}
			sudo apt-get purge bedtools
			curl -OL https://github.com/arq5x/bedtools2/releases/download/v2.18.0/bedtools-2.18.0.tar.gz
			tar -xvzf bedtools-2.18.0.tar.gz
			cd bedtools-2.18.0
			sudo make > /dev/null
			for i in `ls bin/`; do sudo ln -sf ${APPLICATIONS}/bedtools-2.18.0/bin/$i /usr/bin/$i; done
		else
			if [[ `dpkg -s $i 2> /dev/null | wc -l ` > 0 ]]; then
				printf "[ ${GREEN}OK${NC} ] $i already exists\n"
			else
				printf "[ ${RED}NO${NC} ] $i not found, Installing ${i}\n"
				installDepend "$i"
			fi
		fi
	done

	echo
	checkAllPacks

	#DOS/Mac to Unix and vice versa text file format converter all script in integration_analysis/src/
	echo
	cd integration_analysis/src/
	for e in `ls`; do 
		sudo dos2unix $e 2> /dev/null
		printf "[ ${PURPLE}converting${NC} ] converting file ${e} to Unix format\n" 
	done 

	echo
	cd ${APPLICATIONS}/scripts 

	sudo ln -sf /usr/bin/bwa /usr/bin/bwa-stable 2> /dev/null
	printf "[ ${CYAN}link${NC} ] /usr/bin/bwa /usr/bin/bwa-stable\n"
	echo '#!/bin/bash' > ${BIN}/FilterSamReads
	echo 'picard-tools FilterSamReads $@' >> ${BIN}/FilterSamReads
	chmod +x ${BIN}/FilterSamReads 2> /dev/null
	echo '#!/bin/bash' > ${BIN}/MergeSamFiles
	echo 'picard-tools MergeSamFiles $@' >> ${BIN}/MergeSamFiles
	chmod +x ${BIN}/MergeSamFiles 2> /dev/null


###################### CONFIGURATION ######################

	echo ""
	echo "STEP 2: CONFIGURATION"
	echo ""

	printf "[ ${BLUE}new folder${NC} ] Folder creation --> ${BIN}\n"
		sudo mkdir -p ${BIN} 2> /dev/null
		sudo chmod -R 777 ${APPLICATIONS} 2> /dev/null
	printf "[ ${BLUE}new folder${NC} ] Folder creation --> ${PATHGENOME}\n"
		sudo mkdir -p ${PATHGENOME} 2> /dev/null
		sudo chmod -R 777 ${PATHGENOME} 2> /dev/null
	printf "[ ${BLUE}new folder${NC} ] Folder creation --> /opt/NGS/results\n"
		sudo mkdir -p /opt/NGS/results 2> /dev/null
		sudo chmod -R 777 /opt/NGS 2> /dev/null

###################### Download Genome ######################

	echo ""
	echo "STEP 3: Download Genome"
	echo ""

	#file annotation step
	cd ${VISPA2}/annotation
	tar -xvzf ucsc.hg19.refSeq.gtf.tar.gz
	sudo mv hg19.refGene.TIGET.gtf ucsc.hg19.refSeq.gtf
	
	if [ -z $to_Download ]; then
		echo "No Genome will be downloaded, SKIP to step 4..."
		echo
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
		sudo tar -xvzf $file_name
		echo ""
		echo "File Concatenation..."
		
		chr_number=`ls -vI "*_*" | grep -E "chr[[:digit:]]" | wc -l`
		
		for i in $(seq 1 $chr_number); do
			cat chr$i.fa >> $out_name.fa
		done

		cat chrX.fa >> $out_name.fa
		cat chrY.fa >> $out_name.fa
		cat chrM.fa >> $out_name.fa
		
		echo "...."

		printf "[ ${GREEN}CONTATENATION COMPLETE${NC} ]\n"
		echo ""
		echo "Removing extracted chr files...."
		sudo mv $file_name ../
		sudo rm chr*
		sudo mv ../$file_name ${FOLDERGENOME}
		echo ""
		cd ${FOLDERGENOME}/annotation																											#pwd ----> /opt/genome/species/type/annotation
		
		if [[ $out_name == "hg19" ]]; then
			cp ${VISPA2}/annotation/ucsc.hg19.refSeq.gtf.tar.gz .
			tar -xvzf ${FOLDERGENOME}/annotation/ucsc.hg19.refSeq.gtf.tar.gz
			sudo mv hg19.refGene.TIGET.gtf ucsc.hg19.refSeq.gtf
		else
			cp ${VISPA2}/annotation/ucsc.$out_name.refSeq.gtf.gz .
			gunzip ${FOLDERGENOME}/annotation/ucsc.$out_name.refSeq.gtf.gz
		fi

		printf "[ ${BLUE}new folder${NC} ] Folders creation --> "${FOLDERGENOME}"/index/bwa_7\n"
		sudo mkdir -p ${FOLDERGENOME}/index/bwa_7 2> /dev/null

		cd ${FOLDERGENOME}																														#pwd ----> /opt/genome/species/type/
		echo "moving "$out_name.fa" in "${FOLDERGENOME}"/index/bwa_7"
		mv $out_name.fa index/bwa_7
		echo ""
		printf "[ ${PURPLE}INDEXING${NC} ] Indexing of $out_name.fa...this operation require a lot of time. Better if you go to take some coffee.\n"
		cd ${FOLDERGENOME}/index/bwa_7																											#pwd ----> /opt/genome/species/type/index/bwa_7

		bwa-stable index -a bwtsw ${FOLDERGENOME}/index/bwa_7/$out_name.fa  #BWA
		samtools faidx $out_name.fa 										#SAMTOOLS
		picard-tools CreateSequenceDictionary R=$out_name.fa O=$out_name.dict 	#PICARD

		printf "[ ${GREEN}INDEXING COMPLETE${NC} ]\n"
	fi
	
	printf "[ ${BLUE}new folder${NC} ] Folders creation --> "${PATHGENOME}"/vector${NC}\n"
	sudo mkdir -p ${PATHGENOME}/vector 2> /dev/null
	printf "[ ${BLUE}new folder${NC} ] Folders creation --> "${PATHGENOME}"/control/phix${NC}\n"
	sudo mkdir -p ${PATHGENOME}/control/phix 2> /dev/null

	sudo mkdir -p /opt/genome/control/phix174/bwa_7/
	sudo mkdir -p /opt/genome/vector/lv/bwa_7/
	sudo mkdir -p /opt/genome/vector/alu/
	sudo cp -r /opt/applications/scripts/vispa2/genomes/phix/bwa7/* /opt/genome/control/phix174/bwa_7/
	sudo cp -r /opt/applications/scripts/vispa2/genomes/lv/bwa_7/* /opt/genome/vector/lv/bwa_7/
	sudo cp -r /opt/applications/scripts/vispa2/genomes/alu/* /opt/genome/vector/alu/

	cd /opt/genome/vector/lv/bwa_7/

	for k in $( ls *.fa ) ; do 
		INAME=`basename $k | sed 's/.fa//g'`; 
		echo $INAME;
		sudo rm ${INAME}.dict
		bwa-stable index $k; 
		samtools faidx $k;
		picard-tools CreateSequenceDictionary R=${k} O=${INAME}.dict; 
	done
###################### DATABASE CONFIGURATION ######################

	#VISPA2 uses MySQL as DBMS.
	echo
	echo "Create three users (all privileges and readonly and all PRIVILEGES andrea)"
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

	if grep -Rq "LD_LIBRARY_PATH=/path/FlexbarDir:" /home/`whoami`/.bashrc ; then 
		printf "[ ${CYAN}link${NC} ] Flexbar path founded in ./bashrc\n"
	else 
		printf "[ ${CYAN}link${NC} ] Flexbar path not found ./bashrc\n"
		echo 'export LD_LIBRARY_PATH=/path/FlexbarDir:$LD_LIBRARY_PATH' >> ~/.bashrc
		source ~/.bashrc ;
	fi
	
	printf "[ ${CYAN}link${NC} ] /usr/bin/flexbar /usr/bin/flexbar2.5 \n"
	sudo ln -sf /usr/bin/flexbar /usr/bin/flexbar2.5 2> /dev/null
	printf "[ ${CYAN}link${NC} ] ${BIN}/trimmomatic /usr/bin/trimmomatic \n"
	sudo ln -sf ${BIN}/trimmomatic /usr/bin/trimmomatic 2> /dev/null
	printf "[ ${CYAN}link${NC} ] ${BIN}/FilterSamReads /usr/bin/FilterSamReads \n"
	sudo ln -sf ${BIN}/FilterSamReads /usr/bin/FilterSamReads 2> /dev/null
	printf "[ ${CYAN}link${NC} ] ${BIN}/MergeSamFiles /usr/bin/MergeSamFiles \n"
	sudo ln -sf ${BIN}/MergeSamFiles /usr/bin/MergeSamFiles 2> /dev/null

	# Python
	# The following packages are required:

	echo ""
	echo "STEP 5: PYTHON CONFIGURATION"
	echo ""

	listPypack=(MySQL-python pybedtools biopython HTSeq rpy2==2.7.8 scipy numpy matplotlib xlsxwriter pandas pysam==0.7.7)
	for i in ${listPypack[@]}; do
		printf "[ ${GREEN}add python pack${NC} ] python pack ${i}\n" 
		installPyPack "$i"
	done

	cd $APPLICATIONS 																													#pwd ----> /opt/applications
	
	echo "" 
	echo "Python program link creation in /usr/bin"
	echo "" 

	printf "[ ${CYAN}link${NC} ] ${VISPA2}pipeline/illumina/VISPA2.IlluminaMiSeq.pipeline.sh /usr/bin/vispa2\n"
	sudo ln -sf ${VISPA2}/pipeline/illumina/VISPA2.IlluminaMiSeq.pipeline.sh /usr/bin/vispa2 2> /dev/null
	printf "[ ${CYAN}link${NC} ] ${SCRIPT}/annotate_matrix_v2.sh /usr/bin/annotate_matrix\n"
	sudo ln -sf ${SCRIPT}/annotate_matrix_v2.sh /usr/bin/annotate_matrix 2> /dev/null
	printf "[ ${CYAN}link${NC} ] ${SCRIPT}/annotate_bed.py /usr/bin/annotate_bed\n"
	sudo ln -sf ${SCRIPT}/annotate_bed.py /usr/bin/annotate_bed 2> /dev/null
	printf "[ ${CYAN}link${NC} ] ${INTEGRATION}/src/Integration_Analysis.py /usr/bin/create_matrix\n"
	sudo ln -sf ${INTEGRATION}/src/Integration_Analysis.py /usr/bin/create_matrix 2> /dev/null


	printf "[ ${CYAN}link${NC} ] in /usr/bin/import_iss\n"
	sudo ln -sf ${SCRIPT}/import_iss.py /usr/bin/import_iss 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/fqreverseextract.pureheader\n"
	sudo ln -sf ${SCRIPT}/fqreverseextract.pureheader.py /usr/bin/fqreverseextract_pureheader 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/fqextract.pureheader\n"
	sudo ln -sf ${SCRIPT}/fqextract.pureheader.py /usr/bin/fqextract_pureheader 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/rev_extract_header\n"
	sudo ln -sf ${SCRIPT}/rev_extract_header.py /usr/bin/rev_extract_header 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/extract_header\n"
	sudo ln -sf ${SCRIPT}/extract_header.py /usr/bin/extract_header 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/filter_by_cigar_bam\n"
	sudo ln -sf ${SCRIPT}/filter_by_cigar_bam.py /usr/bin/filter_by_cigar_bam 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/filter_by_mate\n"
	sudo ln -sf ${SCRIPT}/filter_by_mate.py /usr/bin/filter_by_mate 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/isa_importrediss_frombed\n"
	sudo ln -sf ${SCRIPT}/dbimport_redundantiss_from_bed.v2.py /usr/bin/isa_importrediss_frombed 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/fastq_qf\n"
	sudo chmod +x ${SCRIPT}/fastq_qf.sh
	sudo ln -sf ${SCRIPT}/fastq_qf.sh /usr/bin/fastq_qf 2> /dev/null
	printf "[ ${CYAN}link${NC} ] in /usr/bin/fasta_to_csv\n"
	sudo chmod +x ${SCRIPT}/fasta_to_csv.rb
	sudo ln -sf ${SCRIPT}/fasta_to_csv.rb /usr/bin/fasta_to_csv 2> /dev/null
	sudo ln -sf ${SCRIPT}/fasta_to_csv.rb /usr/bin/fasta2csv 2> /dev/null
	
	cd /opt
	sudo chmod -R 775 ./*   

printf "[ ${GREEN}SET UP COMPLETED!${NC} ] your machine is ready to vispa2!\n"