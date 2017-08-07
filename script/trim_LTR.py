#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, os, re, argparse, csv, time
import pysam
from operator import itemgetter, attrgetter
from string import Template
#from Bio import SeqIO
from time import gmtime, strftime
from multiprocessing import Process, Lock, Queue
import multiprocessing, math

header = """
+--------------------------------------------------+
        ***   TRIM LTR SEQUENCES   ***
+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		November 2014
 Contact:	andrea.calabria@hsr.it
 Revision:	0.1
 Note:		Trim LTR from aligned sequences and 
 			return only valid headers.
+--------------------------------------------------+
 
 Note:
  - given the alinment of the sequences to the reference LTR,
    (BAM file format), this application will identify the se-
    quences with a valid LTR and trim them all, discarding the
    sequences with an invalid recongnized LTR.

""" 

description = "This application trims the sequences from the mapped LTR (LTR Analysis)."

usage_example = """
Examples of usage:
 (1) 

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('-b', '--bam', dest="bamfile", help="BAM file to process. No default option.", action="store", required=True)
parser.add_argument('-o', '--outfilename', dest="outfilename", help="Output file name. The output will be a list of headers to be removed. No default option.", action="store", required=True)
parser.add_argument('--ltrpair', dest="ltrpair", help="Select the LTR pair. In case of paired-end reads you will have both r1 and r2, else only r1. Choose between r1 and r2. Default: r1.", action="store", default="r1")
parser.add_argument('-r', '--regions', dest="regions", help="Regions to process as CSV string in the format: CHR:START-END. Default value is for HG19: 'LTR:1-32'.", action="store", default="LTR:1-32")
parser.add_argument('-p', '--processes', dest="processes", help="Number of processes to use. Default 2.", action="store", default=2)
parser.add_argument('--minOverlap', dest="minOverlap", help="Minimum overpal of the aligned sequence to the LTR. Default: 10.", action="store", default=10)
parser.add_argument('--maxDistanceFromLastBase', dest="maxDistanceFromLastBase", help="Maximal accepted distance in bp from last base of the LTR (boundaries of the LTR region defined in the option --regions). Default 4.", action="store", default=4)
args = parser.parse_args()

#########################################################################################
####### GLOBAL VARS
#########################################################################################

cigar_flag = {
		'M': 0, #BAM_CMATCH 
		'I': 1, #BAM_CINS 
		'D': 2, #BAM_CDEL 
		'N': 3, #BAM_CREF_SKIP 
		'S': 4, #BAM_CSOFT_CLIP 
		'H': 5, #BAM_CHARD_CLIP 
		'P': 6, #BAM_CPAD 
		'=': 7, #BAM_CEQUAL 
		'X': 8, #BAM_CDIFF 
		} 
cigar_number = {
		0: 'M', #BAM_CMATCH
		1: 'I', #BAM_CINS
		2: 'D', #BAM_CDEL
		3: 'N', #BAM_CREF_SKIP
		4: 'S', #BAM_CSOFT_CLIP
		5: 'H', #BAM_CHARD_CLIP
		6: 'P', #BAM_CPAD
		7: '=', #BAM_CEQUAL
		8: 'X', #BAM_CDIFF
		}

#########################################################################################
### MY FUNCTIONS
#########################################################################################

def checkArgs(args):
	"""
	Check file path and do required pre-operations.
	"""
	if not os.path.isfile(args.bamfile):
		print "\n[AP]\tError while reading files: no valid paths.\n\tInput files:\n\tBAM: %s\n\n\tExiting...\n" %(args.bamfile)
		sys.exit()

def sortBAM(bamfile, sorted_bamfile = None):
	"""
	Sort BAM using external tools -> samtools
	Returns new BAM file name.
	"""
	if sorted_bamfile is None:
		sorted_bamfile = os.path.abspath(bamfile) + "tmpsorted_" + os.path.basename(bamfile) 
	os.system( "samtools sort %s %(s)" %(bamfile, sorted_bamfile) )
	return sorted_bamfile

def indexBAM(bamfile, ):
	"""
	Index BAM using external tools -> samtools
	"""
	os.system( "samtools index %s" %(bamfile) )

def revertString(instr):
	"""
	Given a string in input, revert it: example: AAB5 -> 5BAA
	"""
	return instr[::-1]

def getRegions(regions_str):
	"""
	Given a string CSV of regions in the form CHR:START-END, return the list of tuples of regions
	"""
	# build an array of regions
	regions_list = [] # array of tuples
	for reg in regions_str.split(';'):
		if reg.strip().split(':')[0].strip().isdigit():
			chr = "chr" + reg.strip().split(':')[0].strip()
			start = int(reg.strip().split(':')[1].strip().split('-')[0])
			end = int(reg.strip().split(':')[1].strip().split('-')[1])
			regions_list.append( (chr, start, end) )
		else:
			chr = reg.strip().split(':')[0].strip()
			start = int(reg.strip().split(':')[1].strip().split('-')[0])
			end = int(reg.strip().split(':')[1].strip().split('-')[1])
			regions_list.append( (chr, start, end) )
	return regions_list

def worker(todo_regions, bam, out_q, ): 
	""" 
	The worker function, invoked in a process. 
		region: array of tuples with coordinate to slice chromosome
		bam: the input bam to splice
		out_q: queue of the process
	"""
	product = {} # k = header, v = dict:: "sample", "target", "header", "var_type", "var_start", "var_end"
	for region in todo_regions:
		chr, start, end = region # now start and end are the boundaries of putative targets!
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t...processing region (chr, start, end)", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		
		print "[AP]\t\t\tAcquiring dictionary of reads in the region (chr, start, enc)\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		for bam_alignment in mybam: # for each alignment in the BAM file (both R1 and R2!!!)
			# get orientation and check for correct position (get end from orientation and start)
			orientation = '+' # default, if fwd
			aln_start = bam_alignment.pos # default, if fwd
			aln_end = bam_alignment.pos + bam_alignment.alen # default, if fwd
			if bam_alignment.is_reverse:
				orientation = '-'
				# aln_end = bam_alignment.pos
				# aln_start = bam_alignment.pos + bam_alignment.alen
			# now acquire the DB of reads
			if product.has_key(bam_alignment.qname):
				if not bam_alignment.is_read2: # all but not read 2: this means read 1 and potentially all singletones from r1
					product[bam_alignment.qname]['r1'] = {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
													}
				else:
					product[bam_alignment.qname]['r2'] = {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
															} 
			else: # dict does not contain the key
				if not bam_alignment.is_read2: # all but not read 2: this means read 1 and potentially all singletones from r1
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in ELSE key exists and IF NOT R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
					product[bam_alignment.qname] = {'r1': {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
														}
													}
				else: # the r2 reads
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in ELSE key exists and ELSE R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
					#print bam_alignment.qname
					product[bam_alignment.qname] = {'r2': {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
														}
													}
		print "[AP]\t\t\tDone, region completed (chr, start, end)\t\t\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" 
	#print len(outdata), outdata[0]
	# queue results at the end of this process
	out_q.put(product)

def single_worker(todo_regions, bam): 
	""" 
	The worker function, invoked in a process. 
		region: array of tuples with coordinate to slice chromosome
		bam: the input bam to splice
		out_q: queue of the process
	"""
	product = {} # k = header, v = dict:: "sample", "target", "header", "var_type", "var_start", "var_end"
	for region in todo_regions:
		chr, start, end = region # now start and end are the boundaries of putative targets!
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t...processing region (chr, start, end)", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		
		print "[AP]\t\t\tAcquiring dictionary of reads in the region (chr, start, enc)\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		for bam_alignment in mybam: # for each alignment in the BAM file (both R1 and R2!!!)
			# get orientation and check for correct position (get end from orientation and start)
			orientation = '+' # default, if fwd
			aln_start = bam_alignment.pos # default, if fwd
			aln_end = bam_alignment.pos + bam_alignment.alen # default, if fwd
			if bam_alignment.is_reverse:
				orientation = '-'
				# aln_end = bam_alignment.pos
				# aln_start = bam_alignment.pos + bam_alignment.alen
			# now acquire the DB of reads
			if product.has_key(bam_alignment.qname):
				if not bam_alignment.is_read2: # all but not read 2: this means read 1 and potentially all singletones from r1
					product[bam_alignment.qname]['r1'] = {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
													}
				else:
					product[bam_alignment.qname]['r2'] = {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
															} 
			else: # dict does not contain the key
				if not bam_alignment.is_read2: # all but not read 2: this means read 1 and potentially all singletones from r1
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in ELSE key exists and IF NOT R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
					product[bam_alignment.qname] = {'r1': {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
														}
													}
				else: # the r2 reads
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in ELSE key exists and ELSE R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
					#print bam_alignment.qname
					product[bam_alignment.qname] = {'r2': {	'start': aln_start, 
															'end': aln_end, 
															'strand': orientation , 
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'query': bam_alignment.query,
														}
													}
		print "[AP]\t\t\tDone, region completed (chr, start, end)\t\t\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" 
	#print len(outdata), outdata[0]
	# queue results at the end of this process
	return product


def doMSA(input_fasta):
	"""
	Using the FASTA input file and the reference LTR string (!!, not FASTA format), do the MSA: clustalw and mview pipeline.
	"""
	print "[AP]\tCreating MSA -> Running CLUSTALW"
	out_msa = os.path.splitext(input_fasta)[0]+".msa"
	clustal_cmd = "clustalw2 -ALIGN -TYPE=DNA -OUTFILE=%s -OUTPUT=FASTA -OUTORDER=INPUT -INFILE=%s" %(out_msa , input_fasta, )
	os.system(clustal_cmd)
	print "[AP]\tCreating HTML files from MSA -> Running MVIEW"
	out_html = os.path.splitext(input_fasta)[0]+".msa.html"
	os.system("mview -in pearson -ruler on -html head -css on -coloring any -colormap CLUSTAL %s > %s " %(out_msa, out_html) )


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
	selected_pair = args.ltrpair # from which read do you have the LTR in the correct orientation?
	print "[AP]\tChecking inputs."
	checkArgs(args)
	# if both selections are active (or not), assign the correct output file names
	file_to_writeOccurrenceFasta = None
	file_to_writeOccurrenceFasta_top = None
	file_to_writeFullFasta = None
	if args.writeOccurrenceFasta and args.writeFullFasta:
		file_to_writeFullFasta = os.path.join(os.path.splitext(args.outfilename)[0] + ".all" + os.path.splitext(args.outfilename)[1])
		file_to_writeOccurrenceFasta = os.path.join(os.path.splitext(args.outfilename)[0] + ".unique" + os.path.splitext(args.outfilename)[1])
		file_to_writeOccurrenceFasta_top = os.path.join(os.path.splitext(args.outfilename)[0] + ".top" + os.path.splitext(args.outfilename)[1])
	elif args.writeOccurrenceFasta:
		file_to_writeOccurrenceFasta = args.outfilename
		file_to_writeOccurrenceFasta_top = os.path.join(os.path.splitext(args.outfilename)[0] + ".top" + os.path.splitext(args.outfilename)[1])
	elif args.writeFullFasta:
		file_to_writeFullFasta = args.outfilename
	print "[AP]\tIndexing BAM (this is required if you did not do it)."
	indexing = pysam.index(args.bamfile) # indexing is BLANK
	# acquire BED and BAM files: INPUT data
	print "[AP]\tAcquiring data, both BED and BAM."
	bam = pysam.Samfile(args.bamfile) # the bam file to process
	print "[AP]\tAcquiring region list to use while slicing BAM."
	regions = getRegions(args.regions)

	print "[AP]\tNow looping over BAM alignments Splitting run in processes..."
	##### -------- MULTIPROCESS APPROACH -------- ########
	# out_q = Queue() # init process queue
	# nprocs = int(args.processes)
	# chunksize = int(math.ceil(len(regions) / float(nprocs))) # find the chunk wrt processes for running regions
	# procs = [] # array of processes
	# # run processes in a range of chunks
	# for i in range(nprocs):
	# 	p = multiprocessing.Process(target=worker, args=(regions[chunksize * i:chunksize * (i + 1)], bam, out_q )) # prototype: todo_regions, bam, bed_dictionary, minStartingMatches, out_q
	# 	procs.append(p)
	# 	p.start()
	# # Collect all results into a single result list.
	# print "[AP]\tGetting results from different processes."
	# resultdata = [] # this is what will be written into the file
	# for i in range(nprocs):
	# 	resultdata.append(out_q.get())
	# # Wait for all worker processes to finish
	# print "[AP]\tWaiting that all processes end."
	# for p in procs:
	# 	p.join()
	# union/merge of dictionaries
	# bam_dictionary = {} # general result dictionary
	# for d in resultdata:
	# 	bam_dictionary.update(d)
	# bam_dictionary = {} # general result dictionary
	# for d in resultdata:
	# 	bam_dictionary.update(d)
	##### -------- SINGLE CPU APPROACH -------- ########
	bam_dictionary = single_worker(regions, bam) # prototype: todo_regions, bam, bed_dictionary, minStartingMatches, out_q
	## - then proceed
	
	# write a fasta file with all sequences?
	if file_to_writeFullFasta is not None:
		# write output results from joined array of strings
		print "[AP]\tWriting results into output file", file_to_writeFullFasta
		# selected_pair = 'r1' # which pair contains the LTR in the correct orientation. NOW THIS IS AN OPTION
		filetarget = open(file_to_writeFullFasta, 'w') # init output file
		filetarget.write(">REFERENCE_LTR_SEQUENCE\n%s\n" %(args.referenceLTRseq)) # put as first sequence the reference LTR
		for read_name, pairs in bam_dictionary.iteritems():
			filetarget.write(">%s\n%s\n" %(read_name, pairs[selected_pair]['query']))
		filetarget.close() # finalize file
		# MSA?
		if args.doMSA:
			doMSA(file_to_writeOccurrenceFasta)

	# write a fasta file with the sample sequences with occurrences (sorted by abundance) and the top ones (with abundance_threshold at 0.1)?
	if file_to_writeOccurrenceFasta is not None:
		# now create a file of sequence occurrences (with sample ID)
		print "[AP]\tWriting results of sequence occurrences (aligned query) into output file", file_to_writeOccurrenceFasta
		sequence_dictionary = {} # k = sequence, v = [occurrence, sample_read_id]
		overall_aligned_sequences = 0
		for read_name, pairs in bam_dictionary.iteritems():
			if sequence_dictionary.has_key(pairs[selected_pair]['query']):
				sequence_dictionary[pairs[selected_pair]['query']][0] += 1
				overall_aligned_sequences += 1
			else:
				sequence_dictionary[pairs[selected_pair]['query']] = [1, read_name]
				overall_aligned_sequences += 1
		# acquire a tuple array to use it for sorting purposes
		sequence_ta = [] # tuple array of data that can be sorted
		for sequence_aligned, [occurence, sample_read_id] in sequence_dictionary.iteritems():
			sequence_ta.append(tuple([sample_read_id, sequence_aligned, occurence, 100*float(occurence)/overall_aligned_sequences]))
		sorted_sequence_ta = sorted(sequence_ta, key=itemgetter(2), reverse=True) # sort sequence by occurrences
		filetarget = open(file_to_writeOccurrenceFasta, 'w') # init output file
		filetarget.write(">REFERENCE_LTR_SEQUENCE\n%s\n" %(args.referenceLTRseq))
		### ---- DICTIONARY BASED FILE WRITER
		# for sequence_aligned, [occurence, sample_read_id] in sequence_dictionary.iteritems():
		# 	## uncomment the following line to produce a general CSV tab file
		# 	# filetarget.write("%s\t%s\t%s\n" %(sample_read_id, sequence_aligned, occurence))
		# 	filetarget.write(">%d|%.6f|%s\n%s\n" %(occurence, 100*float(occurence)/overall_aligned_sequences, sample_read_id, sequence_aligned, ))
		### ---- TUPLE ARRAY BASED FILE WRITER
		for sample_read_id, sequence_aligned, occurence, abundance in sorted_sequence_ta:
			## uncomment the following line to produce a general CSV tab file
			# filetarget.write("%s\t%s\t%s\n" %(sample_read_id, sequence_aligned, occurence))
			filetarget.write(">%d|%.6f|%s\n%s\n" %(occurence, abundance, sample_read_id, sequence_aligned, ))
		filetarget.close() # finalize file
		# use the most abundant ones only
		abundance_threshold = 0.1 # do not write data if they are <X%, this is a percentage value!
		filetarget = open(file_to_writeOccurrenceFasta_top, 'w') # init output file
		filetarget.write(">REFERENCE_LTR_SEQUENCE\n%s\n" %(args.referenceLTRseq))
		### ---- TUPLE ARRAY BASED FILE WRITER
		for sample_read_id, sequence_aligned, occurence, abundance in sorted_sequence_ta:
			## uncomment the following line to produce a general CSV tab file
			# filetarget.write("%s\t%s\t%s\n" %(sample_read_id, sequence_aligned, occurence))
			if abundance >= abundance_threshold:
				filetarget.write(">%d|%.6f|%s\n%s\n" %(occurence, abundance, sample_read_id, sequence_aligned, ))
		filetarget.close() # finalize file
		# MSA?
		if args.doMSA:
			doMSA(file_to_writeOccurrenceFasta_top)

	elapsed_time = time.time() - start_time # get read of finishing time
	print "\n[AP]\tTask Finished, closing.\n\tElapsed time: %.2f [seconds]\n" %(elapsed_time)
	

# sentinel
if __name__ == "__main__":
    main()
