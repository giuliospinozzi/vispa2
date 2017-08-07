#!/usr/bin/python
# -*- coding: utf-8 -*-
from subprocess import call
import sys, os, re, argparse, csv
from time import gmtime, strftime

header = """

+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		May 2013
 Contact:	andrea.calabria@hsr.it
 Revision:	0.1
+--------------------------------------------------+

 Description
  - This program automates the annotation in this way:
  	Preconditions:
  	 1. given a GTF file (-> if missed, just do it yourself: mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "select * from ensGene;" mm9 | cut -f2- | genePredToGtf file stdin ensGene.gtf )(or in TIGET alien server: mysql --user=annotation --password=annotation --host=172.25.39.2 -A -N -e "select * from refGene;" ucsc_hg19 | cut -f2- | genePredToGtf file stdin refGene.gtf )( genePredToGtf:: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf ) and the input BED
  	 2. BEDTools globally installed
  	 3. data come from UCSC -> so that we are always aligned to genomic reference (that we downloaded from UCSC)
  	Steps:
  	 1. run annotation with system call to bedtools
  	 2. parse output moving the SOURCE field into DETAILS field and replacing it as GENE_NAME
  - Matrix input in a standard matrix file format
  	 
 TODO:
  - Nothing yet
 
""" 

description = "This application annotates input BED file"

usage_example = """
Examples of usage:
 (1) Annotate BED given refGene DB as GTF file: 
    <appname> 

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ Annotate BED file ] \n", description = description)
parser.add_argument('-b', '--inBed', dest="inBed", help="Input BED file. No default option.", action="store", required=True)
parser.add_argument('-a', '--annotation', dest="annotation", help="Annotation file. Couls be GTF, GFF or BED file format, compatible with BEDtools. No default option.", action="store", required=True)
parser.add_argument('-o', '--outBed', dest="outBed", help="Output file name of the BED file. No default option.", action="store", required=True)
parser.add_argument('--bedtoolsOptions', dest="bedtoolsOptions", help="bedtools closest options. Write a csv string with desired options as you were typing them to run bedtools, without the input and output options (of course!) that is the -a and -b parameters. I.e.: Output file name of the BED file. Default: '-d,-D,b,-t,first'.", action="store", required=False, default="-d,-D,b,-t,first")
parser.add_argument('--matrix', dest="matrix", action="store_true", help="Is the input file directly the matrix file? If yes, activate this option and the aprogram will automatically generate the annotated BED file from the matrix. Default: 'False'; alternative: 'True' or just activate the option.", default=False)

args = parser.parse_args()

########################################################################
####### GLOBAL VARS
########################################################################

#######################################################################
####### FUNCTIONS
########################################################################

def checkArgs(args):
	"""
	Check file path
	"""
	if not os.path.isfile(args.inBed) and not os.path.isfile(args.annotation):
		print "\n[AP]\tError while reading files: no valid paths.\n\tExit\n"
		sys.exit()

def runBedTools(inbed, outbed, annotation, optionstring):
	"""
	Run BEDtools as desired
	"""
	whattorun = ["bedtools", "closest"] + optionstring.split(',') + ['-a', inbed, '-b', annotation, ' > ', outbed]
# 	whattorun = ["bedtools", "closest"] + optionstring.split(',') + ['-a', inbed, '-b', annotation, ]
# 	whattorun = ["closestBed"] + optionstring.split(',') + ['-a', inbed, '-b', annotation, ]
# 	bedout = open(outbed, 'w')
# 	ok = call(whattorun, shell=True, stdout=bedout)
#	print ok, ' '.join(str(x) for x in whattorun)
# 	if ok is not 0:
# 		print "[AP]\tERROR:\tProblems invoking BEDtools, check it please (options included)."
# 		sys.exit()
# 	bedout.close()

	os.system(' '.join(str(x) for x in whattorun))

def rewriteBed(inbed, outbed, delimiter = "\t", gene_name_position = 7):
	"""
	Rewrite BED, delimiter is spacer for CSV.
	gene_name_position -> 0-based position in which to write gene name.
	"""
	# initialize final out file
	bedout = csv.writer(open(outbed, 'wb'), delimiter=delimiter, quotechar=None)
	# read input and process it
	with open(inbed, 'rb') as inf:
		reader = csv.reader(inf, delimiter = delimiter)
		for row in reader:
			gene_name = "NOT FOUND!!!!!!!"
			for col in row:
				if col.startswith("gene_"):
					#print col
					annotation = col.split(';')
				 	gene_name = annotation[0].strip().split()[1].strip()
			if gene_name is not "NOT FOUND!!!!!!!":
				## first remove not significant columns, so far the '.' and chr repetition
				colchrtoremove_found, colchrtoremove = findColNumber('chr', row, scan_from = 4)
				if colchrtoremove_found:
					row.pop(colchrtoremove)
				#coldottoremove_found, coldottoremove = findColNumber('.', row, scan_from = 4)
				#if coldottoremove_found:
				#	row.pop(coldottoremove)
				## insert gene annotation at given position
				row.insert(gene_name_position-1, gene_name.replace('"',''))
				bedout.writerow( row )
			else:
				print "[AP]\tERROR:\tProblems in annotation: gene_name not found in this row:", row
				sys.exit()

def findColNumber(pattern, content_array, scan_from = 4):
	"""
	Fiven a string as pattern (not reg exp, only string) find the corresponding column number (0-based) through 'startswith' and return number and Bool val.
	scan_from -> start looking at col number from this slice. it is useful for BED, keeping first N cols but removing others
	"""
	found = False
	colnumber = None
	for k in range(scan_from, len(content_array)):
		if content_array[k].startswith(pattern):
			colnumber = k
			found = True
	return found, colnumber

def createBedFromMatrix(matrixfile, ):
	"""
	Invoke AWK to format the standard MATRIX file. If this file format changes, pelase change this function!!!
	"""
	filetoreturn = None
	os.system("awk '{print \"chr\"$1\"\t\"$2\"\t\"$2}' '%s' | tail -n +2 > %s" %(matrixfile, filetoreturn) )
	return filetoreturn
	

#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	# first check args and file paths
	checkArgs(args)
	# create the tmp file
	tmpfile = os.path.join("/tmp", os.path.basename(args.outBed) + "." + strftime("%Y%m%d%H%M%S", gmtime()),)
	# run bedtools
	print "[AP]\tRun BEDtools."
	# if input file is BED (not the matrix file)
	if not args.matrix:
		# run bed tools
		statusok = runBedTools(args.inBed, tmpfile, args.annotation, args.bedtoolsOptions)
	else:
		# run bed tools
		inBed = createBedFromMatrix(args.inBed)
		statusok = runBedTools(inBed, tmpfile, args.annotation, args.bedtoolsOptions)
	# rewrite BED
	print "[AP]\tRewrite BED file and produce output file."
	rewriteBed(tmpfile, args.outBed, delimiter = "\t", gene_name_position = 7)
	# clean up files
	os.remove(tmpfile)

# sentinel
if __name__ == "__main__":
    main()
