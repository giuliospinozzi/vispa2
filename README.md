# README #

VISPA2: A Scalable Pipeline for High-Throughput Identification and Annotation of Vector Integration Sites.

Authors: Giulio Spinozzi, Andrea Calabria.

Online web site for demo purpose: http://openserver.itb.cnr.it/vispa


## What is this repository for? ##

* Quick summary
Tool for the analysis of retroviral vector integration sites.

* Version
Versione 2.


## How do I get set up? ##

The full research article reports a short manual for tool installation and configuration.

Further details for the configuration and running are updated here.

### Summary of set up ###
The Bash version of VISPA2 is reachable at: https://bitbucket.org/andreacalabria/vispa2

It can process only paired–end Illumina sequencing reads. The single end mode is present (with also the paired-end) in the web version, at: http://openserver.itb.cnr.it/vispa/ 

### Dependencies ###
Before proceeding with the next steps make sure you install the software first, making sure these are already in the path.

* [MySQL server] (https://dev.mysql.com/downloads/mysql/)
* [Pigz] (https://zlib.net/pigz/)
* [GNU Parallel] (https://www.gnu.org/software/parallel/)
* [Fastqc] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bwa] (http://bio-bwa.sourceforge.net/bwa.shtml)
* [Samtools] (http://www.htslib.org/doc/samtools.html)
* [Trimmomatic] (www.usadellab.org/cms/?page=trimmomatic)
* [EA-Utils] (https://expressionanalysis.github.io/ea-utils/)
* [flexbar] (https://github.com/seqan/flexbar)
* [bamtools] (https://github.com/pezmaster31/bamtools)
* [FilterSamReads/MergeSamFiles] (https://broadinstitute.github.io/picard/)
* [Bedtools] (https://github.com/arq5x/bedtools2/releases)

### Configuration
#### Directories ####
In /opt/applications/bin you should install all third-party software like picard, ea-utils… while in /opt/NGS/results the final result folder (bam, bed…)

```
sudo mkdir /opt/applications

sudo chmod -R 777 /opt/applications

mkdir /opt/applications/bin

sudo mkdir /opt/genome

sudo chmod -R 777 /opt/genome

sudo mkdir /opt/NGS

sudo chmod -R 777 /opt/NGS

mkdir /opt/NGS/results
```


#### Genomes ####
You can download from UCSC the FASTA file (index) and annotation file in GTF and save in the following paths:

```
mkdir /opt/genome/human

mkdir /opt/genome/human/hg19/

mkdir /opt/genome/human/hg19/annotation

mkdir /opt/genome/human/hg19/index

mkdir /opt/genome/human/hg19/index/bwa_7

mkdir /opt/genome/vector
...
mkdir /opt/genome/control

mkdir /opt/genome/control/phix
...
```
In /opt/genome/human/hg19/index/bwa_7 must be inserted the genome (FASTA) and its indexes. The indexes can be built in this way:

```
bwa index -a bwtsw REF.fa

samtools faidx REF.fa

java -jar /opt/applications/bin/picard/picard-tools-1.79/CreateSequenceDictionary.jar R=CE.cns.fa O=CE.cns.dict
```

The same for vector genomes.



### Database configuration
VISPA2 uses MySQL as DBMS.

Install the dev lib and create two users (all privileges and readonly)

```
sudo apt-get install libmysqlclient-dev

mysql -uroot -p

mysql> GRANT ALL PRIVILEGES ON *.* To 'vispa2'@'localhost' IDENTIFIED BY 'vispa2';

mysql> GRANT SELECT ON *.* TO 'readonly'@'localhost' IDENTIFIED BY 'readonlypswd';
```


### Software configuration
VISPA2 requires specific links to run the programs, here the creation.

#### BWA
```
sudo ln -s /opt/applications/bin/bwa/bwa-0.7.15/bwa /usr/bin/bwa-stable
```

#### Flexbar (Version 2.5)
Put this line in you path:
```
export LD_LIBRARY_PATH=/path/FlexbarDir:$LD_LIBRARY_PATH
```
And the specific link to this version:
```
sudo ln -s /opt/applications/bin/flexbar/flexbar_v2.5/flexbar /usr/bin/flexbar2.5
```

#### Trimmomatic
```
sudo ln -s /opt/applications/scripts/isatk/utils/trimmomatic.sh /usr/bin/trimmomatic
```
The bash script, trimmomatic.sh, contains:
```
#!/bin/bash

java -jar /opt/applications/bin/trimmomatic/trimmomatic-0.36.jar $@

```

#### Picard
```
sudo ln -s /opt/applications/scripts/isatk/utils/FilterSamReads.sh /usr/bin/FilterSamReads
```
The bash script, FilterSamReads.sh, contains:
```
#!/bin/bash

picard FilterSamReads $@
```
```
sudo ln -s /opt/applications/scripts/isatk/utils/MergeSamReads.sh /usr/bin/MergeSamReads
```
The bash script, MergeSamReads.sh, contains:
```
#!/bin/bash

picard FilterSamReads $@
```


#### Python
The following packages are required:
```
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
```

#### Repositoy configuration
Cloning
```
cd /opt/applications

hg clone –b ‘v3’ https://bitbucket.org/andreacalabria/vispa2
```
Linking the programs in /usr/local/bin
```
sudo ln -s /opt/applications/scripts/isatk/script/import_iss.py /usr/bin/import_iss

sudo ln -s /opt/applications/scripts/isatk/script/fqreverseextract.pureheader.py /usr/bin/fqreverseextract.pureheader

sudo ln -s /opt/applications/scripts/isatk/script/fqextract.pureheader.py /usr/bin/fqextract.pureheader

sudo ln -s /opt/applications/scripts/isatk/script/rev_extract_header.py /usr/bin/rev_extract_header

sudo ln -s /opt/applications/scripts/isatk/script/extract_header.py /usr/bin/extract_header

sudo ln -s /opt/applications/scripts/isatk/script/filter_by_cigar_bam.py /usr/bin/filter_by_cigar_bam

sudo ln -s /opt/applications/scripts/isatk/script/filter_by_mate.py /usr/bin/filter_by_mate

sudo ln -s /opt/applications/scripts/isatk/script/dbimport_redundantiss_from_bed.v2.py /usr/bin/isa_importrediss_frombed

sudo ln -s /opt/applications/scripts/isatk/script/annotate_matrix_v2.sh /usr/bin/annotate_matrix

sudo ln –s /opt/applications/scripts/isatk/script/fastq_qf.sh /usr/bin/fastq_qf

sudo ln –s /opt/applications/scripts/isatk/script/fasta_to_csv.rb /usr/bin/fasta_to_csv
```

### How to run tests
### Deployment instructions


### Who do I talk to? ###

* Repo owner or admin
Giulio Spinozzi (spinozzi.giulio@hsr.it)
Andrea Calabria (calabria.andrea@hsr.it)

