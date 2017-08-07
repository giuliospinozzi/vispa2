#!/usr/bin/python
# -*- coding: utf-8 -*-

from subprocess import call
import math
import sys, os, re, argparse, csv, time, random
from time import gmtime, strftime

header = """
+--------------------------------------------------+
   ***           LTR ANALYSIS                ***
   ***   DO SIMULATED DATA FROM REAL CASES   ***
+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		November 2014
 Contact:	andrea.calabria@hsr.it
 Revision:	0.1
 Note:		
+--------------------------------------------------+

Steps:
 - parse the PREFIX.info.csv file generated by the program
   'fastaquery_from_bam' and look for the fields (order is
   important):
    #sample_read_id
    #sequence_aligned
    #occurence
    #abundance
    #chr
    #sample_read_start
    #sample_read_end
   NB: the program will fail if sample_read_end is not the 
       las column and if the first 4 ones are the same (in
       the same order).
 - create the GROUND TRUTH file and the FASTA/Q file for
   the simulation/optimization program
""" 

description = "Parse readscount file and produce TSV data file to be plotted."

usage_example = """
Examples of usage:
 (1) ./app -i TEST.LTRgenome.rawR1.alm.sorted.md.alignedltr.fa.info.csv -o /home/andrea/trimming_analysis/SIMFOMINFO_output --genomicSeq AGAGATCAAGTCTCACTATGTTGCCCAGGTTGGTCTCGAACTCTTGGGCTCAAGCAATCCTCTCACCTCAGCCTCCCAAAGTGCTGGGATTACAGACATGAGCCACCCTGCTCGGCTAGAATT --adapterSeq ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA --maxAcceptedDistanceFrom3end 3 --overallSimulatedReads 20000

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('-i', '--infofile', dest="infofile", help="Input file: the info file generated by the program fastaquery_from_bam.", action="store", required=True)
parser.add_argument('-o', '--outPrefix', dest="outPrefix", help="The output prefix for the files (with full path!): (1) ground truth, (2) FASTA of sequences, (3) FASTQ of sequences.", action="store", required=True)
parser.add_argument('--genomicSeq', dest="genomicSeq", help="Genomic reference sequence used in the ground truth file. No default option. Write the sequence only, not FASTA format", action="store", required=True)
parser.add_argument('--adapterSeq', dest="adapterSeq", help="Adapter reference sequence used in the ground truth file, usually it is the LTR sequence. Write the sequence only. No default option.", action="store", required=True)
# parser.add_argument('--adapterFasta', dest="adapterFasta", help="FASTA file of the adapter reference sequence, usually the LTR sequence. Write absolute paths. No default option.", action="store", required=True)
parser.add_argument('--maxAcceptedDistanceFrom3end', dest="maxAcceptedDistanceFrom3end", help="Maximum distance of accepted aligned reads with respect to the LTR sequence length. Default: 3. This value 3 means that in case of an LTR of 32 bp, only the sequences aligned up to the bp 29 (included) are accepted as true, else as false.", action="store", default=3)
parser.add_argument('--overallSimulatedReads', dest="overallSimulatedReads", help="The overall number of simulated reads to return based on the input file (this number will be only a guide for the output file, since using percentages the program will approximate the count to the upper closer integer that means that at least one occurrence for each input read will be returned. Default: 20k, write an integer number.", action="store", required=False, default=20000)

args = parser.parse_args()


#########################################################################################
####### GLOBAL VARS
#########################################################################################


#########################################################################################
### MY FUNCTIONS
#########################################################################################

def checkArgs(args):
	"""
	Check file path and do required pre-operations.
	"""
	if not os.path.isfile(args.infofile):
		print "\n[AP]\tError while reading files: no valid paths.\n\tInput file:\n\tReadsCount: %s\n\n\tExiting...\n" %(args.infofile)
		sys.exit()
	if not os.path.isdir(os.path.split(args.outPrefix)[0]):
		print "\n[AP]\tError: invalid prefix path (destination dir not existing %s)\n\tExiting...\n" %(args.outPrefix)
		sys.exit()
	
def acquireInfo(infofile, ltrseq, genseq, maxAcceptedDistanceFrom3end, overallSimulatedReads):
	"""
	Given the info file, return a tuple array of the sequences to return.
	"""
	seq_tuplist = [] # list of tuples (header, full sequence, ltr sequence, label [0,1])
	with open(infofile, 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			if not row[0].startswith('#'): # avoid the header
				# print row
				ltr_end = row[-1]
				label = 0 # default reject
				header = row[0]
				ltr_sequence = row[1]
				full_sequence = ltr_sequence + genseq
				number_of_replica = int(math.ceil(float(row[3])*overallSimulatedReads/100))
				if len(ltrseq)-int(ltr_end) <= int(maxAcceptedDistanceFrom3end): # to change the label in TRUE
					label = 1
				for replica in range(0,number_of_replica):
					# (header, full sequence, ltr sequence, label [0,1])
					seq_tuplist.append((header+":%d" %(replica), full_sequence, ltr_sequence, label))
	return seq_tuplist

def writeFiles(prefix, seq_tuplist):
	"""
	Given the input seq_tuplist and the prefix of the output files to write, write the fasta, fastq and the ground truth file.
	"""
	fasta_file = prefix + ".fa"
	fastq_file = prefix + ".fq"
	gt_file = prefix + ".gt"
	# gt
	filetarget = open(gt_file, 'w') # init output file
	filetarget.write("#header\tfull sequence\tltr sequence\tlabel\n")
	for read_id, full_sequence, ltr_sequence, label in seq_tuplist:
		filetarget.write("%s\t%s\t%s\t%s\n" %(read_id, full_sequence, ltr_sequence, label))
	filetarget.close() # finalize file
	# fasta_file
	filetarget = open(fasta_file, 'w') # init output file
	for read_id, full_sequence, ltr_sequence, label in seq_tuplist:
		filetarget.write(">%s\n%s\n" %(read_id, full_sequence))
	filetarget.close() # finalize file
	# fastq_file
	filetarget = open(fastq_file, 'w') # init output file
	for read_id, full_sequence, ltr_sequence, label in seq_tuplist:
		filetarget.write("@%s\n%s\n+\n%s\n" %(read_id, full_sequence, '9'*len(full_sequence)))
	filetarget.close() # finalize file
	


#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	start_time = time.time() # get read of starting time
	# first check args and file paths
	print "[AP]\tYour settings:"
	for k, v in vars(args).iteritems():
		print "\t\t%s::\t\t%s" %(k, v)
	
	print "[AP]\tChecking inputs."
	checkArgs(args)

	print "[AP]\tNow acquiring reads from INFO file."
	seq_towrite = acquireInfo(args.infofile, args.adapterSeq, args.genomicSeq, int(args.maxAcceptedDistanceFrom3end), int(args.overallSimulatedReads))

	print "[AP]\tWrite FASTA (.fa), FASTQ (.fq) and GROUND TRUTH (.gt) files with prefix", args.outPrefix
	writeFiles(args.outPrefix, seq_towrite)

	elapsed_time = time.time() - start_time # get read of finishing time
	print "\n[AP]\tTask Finished, closing.\n\tElapsed time: %.2f [seconds]\n" %(elapsed_time)
	

# sentinel
if __name__ == "__main__":
    main()



