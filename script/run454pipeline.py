#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, os, re, argparse, csv, time, random, tempfile
from operator import itemgetter, attrgetter
from time import gmtime, strftime

"""
Python program to run the pipeline for 454 data, give:
1. the new version of the 454 pipeline single tag
2. the association file as a whole file in which last columns represents paths and file names

Example:
[08/2013]
./run454pipeline.py -a /home/andrea/Dropbox/TIGET/Workspace/Andrea/MLD/Association/mld02.test.tsv -p /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454pipeline_v1.2.1.sh --fastaRootDir /storage/d1/archive/MLD_WAS_Science2013_Fasta_and_WAS1004/454data --tmpdir /tmp2 --genome /opt/NGS/index/human/hg19/bwa/hg19.nochr0 --disease MLD --patient MLD02 --workingdir /opt/NGS/results --dbproperties localhost.andrea.mld02_pool1_t1 --ltrfile /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LTR.fa --lcfile /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LC.fa --cigargenome hg19 

[04/2014]
./run454pipeline.py -a /home/andrea/Dropbox/TIGET/Workspace/Andrea/MLD/Association/mld02.test.tsv -p /home/andrea/Dropbox/TIGET/Workbench/isatk/pipeline/454/454pipeline.lv.gemini.sh --fastaRootDir /storage/d1/archive/MLD_WAS_Science2013_Fasta_and_WAS1004/454data --tmpdir /tmp2 --genome /opt/NGS/index/human/hg19/bwa/hg19.nochr0 --disease MLD --patient MLD02 --workingdir /opt/NGS/results --dbproperties localhost.andrea.mld02_pool1_t1 --ltrfile /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LTR.fa --lcfile /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LC.fa --cigargenome hg19 --vectorcigargenome lv --suboptimalthreshold 50

"""

header = """

+--------------------------------------------------+
 Author:    Andrea Calabria
 Date:      August 2013 - April 2014
 Contact:   andrea.calabria@hsr.it
 Revision:  0.1 (08/2013) rerun Science
            0.2 (04/2014) upgraded to BWA-MEM pipe
+--------------------------------------------------+

 Description
  - 454 pipeline at august 2013 requires these variables:
  /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454pipeline_v1.2.1.sh ${DISEASE} ${PATIENT} ${NGSWORKINGPATH} ${FASTA} ${POOL} ${BARCODELIST} ${GENOME} ${TMPDIR} ${ASSOCIATIONFILE} ${DBHOSTID} ${DBSCHEMA} ${DBTABLE} ${LTR} ${LC} ${CIGARGENOMEID} ;
	echo "
  - current 454 pipeline requires:

  
 Note:
  - You must know how many input parameters the pipeline requires and which version it is
  
 TODO:
  - 
  
 Steps

""" 

description = "This application will run the new 454 pipeline as a whole."

usage_example = """
Examples of usage:
 (1) 

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('-a', '--associationFile', dest="associationFile", help="The association file. No default option.", action="store", required=True)
parser.add_argument('-p', '--pipeline', dest="pipeline", help="Pipeline path. No default option.", action="store", required=True)
parser.add_argument('--associationFieldsOrder', dest="associationFieldsOrder", help="CSV string decsribing how to re-write the association file that the pipeline needs to reads. Put the field list in terms of numbers (1-based) with respect to the associationFile order. Default option, coming from Fabrizio's order: 4,4,7,5,8,2,10,12,9,6,3 (former 4,4,7,5,8,2,10,12,9,1)", action="store", default="4,4,7,5,8,2,10,12,9,6,3")
parser.add_argument('--poolField', dest="poolField", help="Pool field wrt association file, 1-based number. Default option: 16.", action="store", default='16')
parser.add_argument('--fastaRootDir', dest="fastaRootDir", help="Root directory of the FastA file: this root will be added to the path reported in the association file.", action="store", required=True)
parser.add_argument('--tmpdir', dest="tmpdir", help="Tmp directory. Default: '/tmp", action="store", default="/tmp")
parser.add_argument('--genome', dest="genome", help="Pipeline option. No default option.", action="store", required=True)
parser.add_argument('--disease', dest="disease", help="Pipeline option. No default option.", action="store", required=True)
parser.add_argument('--patient', dest="patient", help="Pipeline option. No default option.", action="store", required=True)
parser.add_argument('--workingdir', dest="workingdir", help="Pipeline option. No default option.", action="store", required=True)
parser.add_argument('--dbproperties', dest="dbproperties", help="Pipeline option. No default option. Write CSV string (comma separated) with: host,schema,table", action="store", required=True)
parser.add_argument('--ltrfile', dest="ltrfile", help="Pipeline option. No default option.", action="store", required=True)
parser.add_argument('--lcfile', dest="lcfile", help="Pipeline option. No default option.", action="store", required=True)
parser.add_argument('--cigargenome', dest="cigargenome", help="Pipeline option. No default option. Value examples: hg19, mm9, mfa5, ce, lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada.", action="store", required=True)
parser.add_argument('--vectorcigargenome', dest="vectorcigargenome", help="Pipeline option. No default option. lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada.", action="store", required=True)
parser.add_argument('--suboptimalthreshold', dest="suboptimalthreshold", help="Pipeline option. Default: 50.", action="store", required=False, default="50")
#parser.add_argument('--logfile', dest="logfile", help="Pipeline option. No default option.", action="store", required=True)

args = parser.parse_args()

#########################################################################################
####### GLOBAL VARS
#########################################################################################
timeprefix = strftime("%Y%m%d%H%M%S", gmtime())

#########################################################################################
### MY FUNCTIONS
#########################################################################################

def checkArgs(args):
	"""
	Check file path and do required pre-operations.
	"""
	if not os.path.isfile(args.associationFile) or not os.path.isfile(args.genome) or not os.path.isfile(args.ltrfile) or not os.path.isfile(args.lcfile) or not os.path.isdir(args.workingdir):
		print "\n[AP]\tError while reading files: no valid paths.\n\tExiting...\n" 
		sys.exit()
	if not os.path.isdir(args.tmpdir):
		os.makedirs(args.tmpdir)

def getAssociationInfo(assofile, delimiter = "\t"):
	"""
	Acquire data.
	K = last 2 columns of the file
	"""
	dictionary = {} # k = (folder, file fasta), v = array of values from the row
	with open(assofile, 'rb') as inf:
		reader = csv.reader(inf, delimiter = delimiter)
		for row in reader:
			if len(row) > 2:
				dictionary[(row[-2], row[-1])] = row
	return dictionary

def writePoolAssociationFile(tmpfile, rowarray, order):
	"""
	
	"""
	fo = open(tmpfile, 'w')
	content = []
	for index in order.split(','):
		content.append(rowarray[int(index)-1])
	fo.write('\t'.join(x for x in content))
	fo.close()
	
	
#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	# first check args and file paths
	print "[AP]\tChecking inputs."
	checkArgs(args)
	print "[AP]\tAcquire association file data."
	dict_pool = getAssociationInfo(args.associationFile)
	print "[AP]\tAssign values to vars."
	for (poolfolder, fastaname), rowdata in dict_pool.iteritems(): # for each pool and file, set vars and run the pipeline
		poolassofile = os.path.join(tempfile.gettempdir(), timeprefix + str(random.randint(10, 99)) + "." + fastaname + ".tsv")
		writePoolAssociationFile(poolassofile, rowdata, args.associationFieldsOrder)
		print "[AP]\tRun pipeline for each pool."
		## /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454pipeline_v1.2.1.sh ${DISEASE} ${PATIENT} ${NGSWORKINGPATH} ${FASTA} ${POOL} ${BARCODELIST} ${GENOME} ${TMPDIR} ${ASSOCIATIONFILE} ${DBHOSTID} ${DBSCHEMA} ${DBTABLE} ${LTR} ${LC} ${CIGARGENOMEID} ## >> ${LOGF} ;
		os.system("%(pipeline)s %(disease)s %(patient)s %(workingdir)s %(fasta)s %(pool)s %(barcodelist)s %(genome)s %(tmpdir)s %(poolassofile)s %(dbhost)s %(dbschema)s %(dbtable)s %(ltrfile)s %(lcfile)s %(cigargenome)s %(vectorcigargenome)s %(suboptimalthreshold)s " 
			%{	'pipeline': args.pipeline, #0
				'disease': args.disease, #1
				'patient': str(rowdata[0]), #args.patient, #2
				'pool':	rowdata[int(args.poolField)-1].split('/')[-1].strip(), #3
				'genome': args.genome, #4
				'barcodelist': 'barcodelist', #5
				'workingdir': args.workingdir, #6
				'fasta': os.path.join(args.fastaRootDir, os.path.join(poolfolder, fastaname) ), #7
				'tmpdir': args.tmpdir, #8
				'poolassofile': poolassofile, #9
				'dbhost': args.dbproperties.split('.')[0], #10
				'dbschema': args.dbproperties.split('.')[1], #11
				'dbtable': args.dbproperties.split('.')[2], #12
				'ltrfile': args.ltrfile, #13
				'lcfile': args.lcfile, #14
				'cigargenome': args.cigargenome, #15
				'vectorcigargenome': args.vectorcigargenome, #16
				'suboptimalthreshold': args.suboptimalthreshold, #17
				}
			)

# sentinel
if __name__ == "__main__":
    main()
