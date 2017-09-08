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
* [MySQL server] (https://dev.mysql.com/downloads/mysql/)
* [Fastqc] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bwa] (http://bio-bwa.sourceforge.net/bwa.shtml)
* [Samtools] (http://www.htslib.org/doc/samtools.html)
* [Trimmomatic] (www.usadellab.org/cms/?page=trimmomatic)
* [fastq-multx] (https://expressionanalysis.github.io/ea-utils/)
* [flexbar] (https://github.com/seqan/flexbar)
* [bamtools] (https://github.com/pezmaster31/bamtools)
* [FilterSamReads/MergeSamFiles] (https://broadinstitute.github.io/picard/)
* [Bedtools] (https://github.com/arq5x/bedtools2/releases)

### Configuration
#### Directories ####
In /opt/applications/bin you should install all third-party software like picard, ea-utils… while in /opt/NGS/results the final result folder (bam, bed…)

```
sudo mkdir /opt/applications
```
```
sudo chmod -R 777 /opt/applications
```
```
mkdir /opt/applications/bin
```
```
sudo mkdir /opt/genome
```
```
sudo chmod -R 777 /opt/genome
```
```
sudo mkdir /opt/NGS
```
```
sudo chmod -R 777 /opt/NGS
```
```
mkdir /opt/NGS/results
```


#### Genomes ####
You can download from UCSC the FASTA file (index) and annotation file in GTF and save in the following paths:

```
mkdir /opt/genome/human
```
```
mkdir /opt/genome/human/hg19/
```
```
mkdir /opt/genome/human/hg19/annotation
```
```
mkdir /opt/genome/human/hg19/index
```
```
mkdir /opt/genome/human/hg19/index/bwa_7
```
```
mkdir /opt/genome/vector
```
In /opt/genome/human/hg19/index/bwa_7 must be inserted the genome (FASTA) and its indexes. The indexes can be built in this way:

```
bwa index -a bwtsw REF.fa
```
```
samtools faidx REF.fa
```
```
java -jar /opt/applications/bin/picard/picard-tools-1.79/CreateSequenceDictionary.jar R=CE.cns.fa O=CE.cns.dict
```

The same for vector genomes.



### Database configuration
VISPA2 uses MySQL as DBMS.

Install the dev lib and create two users (all privileges and readonly)

```
sudo apt-get install libmysqlclient-dev
```
```
mysql -uroot -p
```
```
mysql> GRANT ALL PRIVILEGES ON *.* To 'vispa2'@'localhost' IDENTIFIED BY 'vispa2';
```
```
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

### How to run tests
### Deployment instructions


### Who do I talk to? ###

* Repo owner or admin
Giulio Spinozzi (spinozzi.giulio@hsr.it)
Andrea Calabria (calabria.andrea@hsr.it)

