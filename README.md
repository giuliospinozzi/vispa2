# README #

VISPA2: A Scalable Pipeline for High-Throughput Identification and Annotation of Vector Integration Sites.

Authors: Giulio Spinozzi, Andrea Calabria.

Online web site for demo purpose: http://openserver.itb.cnr.it/vispa


## What is this repository for? ##

* Quick summary
Tool for the analysis of retroviral vector Integration Sites (IS).

* Version
Vector Integration Site Parallel Analysis, Version 2.


## How do I get set up? ##

The full research article reports a short manual for tool installation and configuration.
An easy configuration script will support the installation of the pipeline VISPA2. Use the following command line statement that will guide you through the installation and set up (as root user):

```
cd VISPA2
# get help from the configuration file
./SetUp_Vispa2.sh -h
# run the installation and configuration of the required tools and download the genome(s), here hg19 and mm10 (this option will also index the reference genome). The script exploits reference genomee common names, as reported in UCSC web site http://hgdownload.soe.ucsc.edu/downloads.html#mouse. Please, enable internet network connections to UCSC web site.
./SetUp_Vispa2.sh -s -g hg19,mm10
```

Details for the configuration and running:

* [Paired-end mode](https://bitbucket.org/andreacalabria/vispa2/wiki/VISPA2-PairedEnd) (Illumina) of VISPA2

* [Single-end mode](https://bitbucket.org/andreacalabria/vispa2/wiki/VISPA2-SingleEnd) (454-like) of VISPA2

* [CreateMatrix](https://bitbucket.org/andreacalabria/vispa2/wiki/VISPA2-IS_Matrix), the program to generate the final matrix file of annotated IS sites.

### Who do I talk to? ###

* Repo owner or admin
Giulio Spinozzi (spinozzi.giulio@hsr.it), Andrea Calabria (calabria.andrea@hsr.it)

