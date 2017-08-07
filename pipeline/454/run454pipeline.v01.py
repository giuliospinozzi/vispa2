#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, os, re, argparse, csv, time
from operator import itemgetter, attrgetter
from time import gmtime, strftime

"""
Python program to run the pipeline for 454 data, give:
1. the new version of the 454 pipeline single tag
2. the association file as a whole file in which last columns represents paths and file names

Example:
./run454pipeline.py -a /home/andrea/Dropbox/TIGET/Workspace/Andrea/MLD/Association/mld02.test.tsv -p /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454pipeline_v1.2.1.sh --fastaRootDir /storage/d1/archive/MLD_WAS_Science2013_Fasta_and_WAS1004/454data --tmpdir /tmp2 --genome /opt/NGS/index/human/hg19/bwa/hg19.nochr0 --disease MLD --patient MLD02 --workingdir /opt/NGS/results --dbproperties localhost.andrea.mld02_pool1_t1 --ltrfile /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LTR.fa --lcfile /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/LC.fa --cigargenome hg19 --gatkgenome /opt/genome/human/hg19/index/bwa/hg19.nochr0.fa 

"""

header = """

+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		August 2013
 Contact:	andrea.calabria@hsr.it
 Revision:	0.1
+--------------------------------------------------+

 Description
  - 454 pipeline at august 2013 requires these variables:
  /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454pipeline_v1.2.1.sh ${DISEASE} ${PATIENT} ${NGSWORKINGPATH} ${FASTA} ${POOL} ${BARCODELIST} ${GENOME} ${TMPDIR} ${ASSOCIATIONFILE} ${DBHOSTID} ${DBSCHEMA} ${DBTABLE} ${LTR} ${LC} ${CIGARGENOMEID} ;
	echo "
  
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
parser.add_argument('--associationFieldsOrder', dest="associationFieldsOrder", help="CSV string decsribing how to re-write the association file that the pipeline needs to reads. Put the field list in terms of numbers (1-based) with respect to the associationFile order. Default option, coming from Fabrizio's order: 4,4,7,5,8,2,10,12,9,1", action="store", default="4,4,7,5,8,2,10,12,9,1")
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
parser.add_argument('--gatkgenome', dest="gatkgenome", help="Pipeline option. GATK genome, No default option. Value examples: /opt/genome/human/hg19/index/bwa/hg19.nochr0.fa.", action="store", required=True)
parser.add_argument('--cigargenome', dest="cigargenome", help="Pipeline option. No default option. Value examples: hg19, mm9.", action="store", required=True)
parser.add_argument('--contaminantDb', dest="contaminantDb", help="Pipeline option. No default option. For example: /opt/NGS/vectors/UniVec.2013.TIGET.LV.ARSA.fa | UniVec.2013.TIGET.LV.WAS.fa | UniVec_Tiget_Gamma.fa.", action="store", required=True)
parser.add_argument('--stringentFilterRegions', dest="stringentFilterRegions", help="Pipeline option. No default option. For example:  /opt/NGS/annotation/hg19/ucsc.simplerepeat.hg19.bed or /opt/NGS/annotation/mm9/ucsc.simplerepeat.mm9.bed.", action="store", required=False)
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
	if not os.path.isfile(args.associationFile) or not os.path.isfile(args.genome) or not os.path.isfile(args.ltrfile) or not os.path.isfile(args.lcfile) or not os.path.isdir(args.workingdir) or not os.path.isfile(args.contaminantDb) or not os.path.isfile(args.gatkgenome):
		print "\n[AP]\tError while reading files: no valid paths.\n\tExiting...\n" 
		sys.exit()

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
		poolassofile = os.path.join(args.tmpdir, timeprefix + "." + fastaname + ".tsv")
		writePoolAssociationFile(poolassofile, rowdata, args.associationFieldsOrder)
		print "[AP]\tRun pipeline for each pool."
		## /home/andrea/Dropbox/TIGET/Workbench/ISAtools/NGS/454pipeline_v1.2.1.sh ${DISEASE} ${PATIENT} ${NGSWORKINGPATH} ${FASTA} ${POOL} ${BARCODELIST} ${GENOME} ${TMPDIR} ${ASSOCIATIONFILE} ${DBHOSTID} ${DBSCHEMA} ${DBTABLE} ${LTR} ${LC} ${CIGARGENOMEID} ## >> ${LOGF} ;
		os.system("%(pipeline)s %(disease)s %(patient)s %(workingdir)s %(fasta)s %(pool)s FOO_barcodelist %(genome)s %(tmpdir)s %(poolassofile)s %(dbhost)s %(dbschema)s %(dbtable)s %(ltrfile)s %(lcfile)s %(cigargenome)s %(contaminantDb)s %(stringentFilterRegions)s %(gatkgenome)s " 
			%{	'pipeline': args.pipeline,
				'disease': args.disease,
				'patient': args.patient,
				'pool':	rowdata[int(args.poolField)-1].split('/')[-1].strip(),
				'genome': args.genome,
				'workingdir': args.workingdir,
				'fasta': os.path.join(args.fastaRootDir, os.path.join(poolfolder, fastaname) ),
				'tmpdir': args.tmpdir,
				'poolassofile': poolassofile,
				'dbhost': args.dbproperties.split(',')[0],
				'dbschema': args.dbproperties.split(',')[1],
				'dbtable': args.dbproperties.split(',')[2],
				'ltrfile': args.ltrfile,
				'lcfile': args.lcfile,
				'cigargenome': args.cigargenome,
				'contaminantDb': args.contaminantDb,
				'stringentFilterRegions': args.stringentFilterRegions,
				'gatkgenome': args.gatkgenome,
				}
			)

# sentinel
if __name__ == "__main__":
    main()
