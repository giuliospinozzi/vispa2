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

The full research article describes the software and its procecures.

An easy configuration script will support users in the process of installation of VISPA2. Use the following command line statement that will guide you through the installation and set up (as root user):

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

<!---
Details for the configuration and running:

* [Paired-end mode](https://bitbucket.org/andreacalabria/vispa2/wiki/VISPA2-PairedEnd) (Illumina) of VISPA2

* [Single-end mode](https://bitbucket.org/andreacalabria/vispa2/wiki/VISPA2-SingleEnd) (454-like) of VISPA2

* [CreateMatrix](https://bitbucket.org/andreacalabria/vispa2/wiki/VISPA2-IS_Matrix), the program to generate the final matrix file of annotated IS sites.

* [Installation Test](https://bitbucket.org/andreacalabria/vispa2/wiki/Test%20VISPA2%20installation%20with%20a%20sample%20dataset) for automated test.

* [How to add new reference genomes](https://bitbucket.org/andreacalabria/vispa2/wiki/How%20to%20add%20a%20new%20reference%20genome)
-->

To run a sample dataset please use the template and data in testing folder.

### Who do I talk to? ###

* Repo owner or admin
Giulio Spinozzi (spinozzi.giulio@hsr.it), Andrea Calabria (calabria.andrea@hsr.it)

