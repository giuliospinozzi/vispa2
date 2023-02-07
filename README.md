# README #

VISPA2: A Scalable Pipeline for High-Throughput Identification and Annotation of Vector Integration Sites.

Authors: Giulio Spinozzi, Andrea Calabria.

Online web site for demo purpose: http://openserver.itb.cnr.it/vispa

Scientific Paper at BMC Bioinformatics: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1937-9


## What is this repository for? ##

* Quick summary
Tool for the analysis of retroviral vector Integration Sites (IS).

* Version
Vector Integration Site Parallel Analysis, Version 2.


## How do I get set up? ##

The full research article describes the software and its procecures.

An easy configuration script will support users in the process of installation of VISPA2. Use the following command line statement that will guide you through the installation and set up (as root user):

If you have a Linux System (Debian based) with Python 2.7 and R 3 you can use the configuration script below:

```
cd vispa2
# get help from the configuration file
./config_vispa2.sh -h
# run the installation and configuration of the required tools and download the genome(s), 
# here hg19 (this option will also index the reference genome). 
# The script exploits reference genome common names, as reported in UCSC web 
# site http://hgdownload.soe.ucsc.edu/downloads.html#mouse. 
# Please, enable internet network connections to UCSC web site.
./config_vispa2.sh -s human -i hg19
```

If you have a newer system please use Conda to install some important packages and comment in the configuration script the relative lines:

```
conda create -n vispa2

conda install -c bioconda rpy2

conda install -c conda-forge/label/cf202003 r-soniclength=1.4.4

conda install -c conda-forge/label/cf202003 pandas=0.17.1

conda install -c bioconda matplotlib

conda install -c conda-forge editdistance mysql-python xlsxwriter

conda install -c bioconda pybedtools pysam=0.7.7 htseq bedtools=2.20.1

```
Change py interpreter in each py script used (/opt/anaconda3/envs/vispa2/bin/python2)



Details for the configuration and running:

* [CreateMatrix](https://github.com/giuliospinozzi/integration_analysis), the program to generate the final matrix file of annotated IS sites.

* [Installation Test](use vispa2_sampleData.sh in testing folder) for automated test.



### Test the usage ###

To run a sample dataset please use the template and data in testing folder.

### Who do I talk to? ###

* Repo owner or admin
Giulio Spinozzi (spinozzi.giulio@hsr.it), Andrea Calabria (calabria.andrea@hsr.it)

