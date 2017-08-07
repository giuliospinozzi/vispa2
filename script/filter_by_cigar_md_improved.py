#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, os, re, argparse, csv
import HTSeq, pybedtools, pysam
from operator import itemgetter, attrgetter
from string import Template
#from Bio import SeqIO
from time import gmtime, strftime

##### NOTE #####
# 1. pybedtools importa bam e bed e processa stringa bene. inoltre scrive BED e fa intersezione up/down-stream con features (provalo!!!!)
#    x = pybedtools.example_bedtool('/opt/NGS/results/qLAMdNEF/test/bam/pool5plus6_t1/LTR22.LC74.sorted.rel.pg.dupflag.tagfilter.bam')
# 2. HTSeq permette lettura BAM e tutti gli altri formati e fa processare stringhe: http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html#tour
#		http://www-huber.embl.de/users/anders/HTSeq/doc/alignments.html#HTSeq.Alignment
#    bam_reader = HTSeq.BAM_Reader('/opt/NGS/results/qLAMdNEF/test/bam/pool5plus6_t1/LTR22.LC74.sorted.rel.pg.dupflag.tagfilter.bam')
#	 for q in bam_reader:
#		print q.read, q.cigar, q.optional_field('MD')
# 	 for a in HTSeq.itertools.islice( bam_reader, 5 ):
# 		print len(a.cigar), a.cigar, a.optional_field('MD'), a.aQual, a.cigar[0].type, a.cigar[0].size, a.cigar[0].ref_iv, a.cigar[0].query_from, a.cigar[0].query_to, a.cigar[0].check(), a.iv.strand, a.iv.start, a.iv.start_d, a.read.name, a.read.seq, a.read.qual, a.pe_which, a.iv.strand
		

### Idea di base: scrivere un programma in cui dai in input un BAM (allineamenti globali), un BED (dove e' specificato l'LTR) e filtri le reads del BED sfruttando il CIGAR e il tag MD. Come prima cosa guardi lo strand, poi valuti il CIGAR. Poi valuti il tag MD. Questo tag (una volta stabilito lo strand, e quindi se iniziare da 5' o 3') puo' essere valutato tramite una REG ex data come input al programma, una opzione; in particolare 2 input, una per il fwd e una per il rev. Cosi' rendi il programma flessibile e versatile.
"""
Flusso:
bam = read BAM
for k in BED
	# controlla nome
	if bam.read.name == k.header:
		# se e' paired end, devi capire quale mate ha LTR; else vai avanti
		if bam.paired_end:
			filtra solo quelle che hanno LTR -> in base al BED sai orientamento e quindi prendi quella nell'orientamento del BED (o del tipo di esperimento, per esempio sempre il first mate [pe_which])
		latoOK: Acquisisci orientamento per sapere dove valutare in modo stringente lo score (se a inizio o fine sequenza)
		Cigar passa filtro? -> ha delezioni/inserzioni/altro non M nel latoOK? se si, elimina riga, else procedi (o printa ID sequenza)
		MD passa filtro? -> ha problemi (mismatches) nelle prime 3 basi? se si, elimina riga (o printa ID sequenza)
"""

header = """

+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		May 2013
 Contact:	andrea.calabria@hsr.it
 Revision:	0.2
+--------------------------------------------------+

 Description
  - Paired Ends ONLY!
  
 Note:
  - all BED reads must be contained into BAM alignments
  - BED format -> derived from BedTools converting BAM to BED -> names are HEADER/1 or /2
  - General MD reg exp: [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)* 
  - Min threshold for MATCHES of MD tag: 3 (<=3 -> discard read)
  - it uses PySAM and not HTSeq -> HTSeq has an error on PAIR TYPE evaluation (aln.pe_which if often None even if it is first read in pair)!!!!! Whereas PySAM does not have problems (aln.is_read1).
  
 TODO:
  - Single read implementation

 Steps
	My filter:: depending on strand, 
	1. check CIGAR: it must start/end with at least 4M else discard read
	2. check MD:
		if MD is digit -> ok: tutti buoni
		else, check the first/end number is <=4 -> discard read
""" 

description = "This application will filter ISs by CIGAR score and MD tag."

usage_example = """
Examples of usage:
 (1) Filter BED by CIGAR and MD tag using BED as driver [from path /home/andrea/Dropbox/TIGET/Workspace/Andrea/qLAMdNEF/_bam_p5e6_t1/]
    APP --bed LTR22.LC74.sorted.rel.pg.dupflag.tagfilter.bed --bam LTR22.LC74.sorted.rel.pg.dupflag.tagfilter.bam -o LTR22.LC74.sorted.rel.pg.dupflag.tagfilter.cigarfilter.v2.bed

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--bed', dest="bedfile", help="BED file to process. No default option.", action="store", required=True)
parser.add_argument('--bam', dest="bamfile", help="BAM file to process. No default option.", action="store", required=True)
parser.add_argument('--minStartingMatches', dest="minStartingMatches", help="Minimum starting matches to consider the read valid [Integer]. If a read has less than minStartingMatches, it will be discarded. Default = 3 (-> >=3 are valid, else discarded).", action="store", required=False, default=3)
#parser.add_argument('--headerformat', dest="headerformat", help="Header format of the filter file. Cases: type 1:: '@M00174:25:000000000-A0T21:1:1:15041:1491 1:N:0:0'; type 2:: '@M00571:5:000000000-A267F:1:1101:14475:1610/1'. Select type number {1, .., N}. No default value assigned because you must explicitely be aware of what you are doing.", action="store", required=True)
#parser.add_argument('--singleread', dest="singleread", action="store_true", help="In case of single read experiment, not paired-end. Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--sortBAM', dest="sortBAM", action="store_true", help="Sort BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--indexBAM', dest="indexBAM", action="store_true", help="Index BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('-o', '--outfilename', dest="outfilename", help="Output file name of the BED file. No default option.", action="store", required=True)

args = parser.parse_args()

#########################################################################################
####### GLOBAL VARS
#########################################################################################

# init db values


#########################################################################################
### MY FUNCTIONS
#########################################################################################

def checkArgs(args):
	"""
	Check file path and do required pre-operations.
	"""
	if not os.path.isfile(args.bedfile) or not os.path.isfile(args.bamfile):
		print "\n[AP]\tError while reading files: no valid paths.\n\tExit\n"
		sys.exit()
	if args.indexBAM: # index BAM file
		indexBAM(args.bamfile)
	if args.sortBAM: # sort BAM and assign the new file to ARGS and index it
		sorted_bamfile = sortBAM(args.bamfile)
		args.bamfile = sorted_bamfile
		indexBAM(args.bamfile)
	
	
def writeBEDfromString(targetfile, filecontent):
	"""
	From a large string, write BED file
	
	!!! TODO: not working fine, repair. Use writeBEDfromStringByHand(targetfile, filecontent) instead. !!!
	
	"""
	targetfile = pybedtools.BedTool(filecontent, from_string=True)

def writeBEDfromStringByHand(targetfile, filecontent):
	"""
	From a large string, write BED file
	"""
	ft = open(targetfile, 'w')
	ft.write(filecontent)
	ft.close()

def getBAMmapping_sliceBased(bed_row, bam_obj):
	"""
	----------------------------
	|	!!!PRECONDITION!!!!    |
	|	BAM MUST BE SORTED!    |
	----------------------------
	Given BAMaln (alignment from HTSeq), look for seqheader and get this row/s
	Slice BAM by position written into the BED element:: 
		get position and let HTSeq acquire them: p = HTSeq.GenomicPosition( "1", 145439814, "+" ) 
		Get window of the slice: window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
		filter BAM: for almnt in sortedbamfile[ window ]: DO...
	"""
	flanking_bases = 100 # interval of bases to extend search
	outarray_bam_aln = []
	# Get window of the slice: window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
	window = HTSeq.GenomicInterval( bed_row.chrom, bed_row.start - flanking_bases, bed_row.end + flanking_bases, bed_row.strand )
	# filter BAM: for almnt in sortedbamfile[ window ]: DO...
	for bam_aln in bam_obj[ window ]:
		if bed_row.name.strip().startswith( bam_aln.read.name ) and bed_row.strand == bam_aln.iv.strand: # if sequence header in BED contains seq header in BAM (because of /1 or /2 in BED header) and the strand is the same -> take this element
			outarray_bam_aln.append(bam_aln)
	# return data
	if len(outarray_bam_aln) == 0:
		print "[AP]\tERROR:\tSequence not found. Looking for sequences from BED interval (%s:%d-%d/%s), no reads found into BAM. This is NOT POSSIBLE!!! Because BED comes from BAM file.\n" %(bed_row.chrom, bed_row.start - flanking_bases, bed_row.end + flanking_bases, bed_row.strand)
		sys.exit()
	elif len(outarray_bam_aln) > 1:
		print "[AP]\tWARNING:\tMore than 1 sequences matched conditions. This is strange since we look for header and strand together! Please check it carefully"
	else:
		return outarray_bam_aln

def evaluateCigar_OBSOLETE(array_bam_aln):
	"""
	Input: array of HTSeq alignments
	Output: array of HTSeq alignments (len our <= len in)
	Based on strand, evaluate possible problems in CIGAR string. Given a as HTSeq BAM reader aln,
	a.cigar[0].type, a.cigar[0].size, a.cigar[0].ref_iv, a.cigar[0].query_from, a.cigar[0].query_to, a.cigar[0].check()
	"""
	valid = False
	out_array_bam_aln = []
	for e in array_bam_aln:
		if e.cigar[0].type is "M": # only if CIGAR is M, go ahead
			valid = True
			out_array_bam_aln.append(e)
	return valid, out_array_bam_aln

def evaluateCigar(bam_aln):
	"""
	Input: pysam alignment object
	Output: true or false
	Based on strand, evaluate possible problems in CIGAR string. 
	Cigar def: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/glossary.html#term-cigar
	Operations and len are referred to main manual: http://samtools.sourceforge.net/SAM1.pdf
	We could take into account number ops: 0, 7, 8 -> M, =, X. So far is only M [TODO]
	"""
	valid = False
	for tag in bam_aln.cigar:
		if tag[0] == 0 and int(tag[1]) > 0: # only if you have M -> ok
			valid = True
	return valid


def evaluateMDtag_OBSOLETE(array_bam_aln, min_starting_matches = 3):
	"""
	Input: array of HTSeq alignments
	Output: array of HTSeq alignments (len our <= len in)
	Based on strand, evaluate possible problems in CIGAR string. Given a as HTSeq BAM reader aln,
	a.optional_field('MD')
	"""
	fwd_basesmatch_re = '\A[0-9]+' # how many correct alignments (bases) do we have at the beginning of the read? -> returns array of max 1 element: {'', 0->}
	rev_basesmatch_re = '[0-9]+\Z' # how many correct alignments (bases) do we have at the end of the read? -> returns array of max 1 element: {'', 0->}
	re_to_use = None # re to be used, one of the previous 2
	valid = False
	out_array_bam_aln = []
	for e in array_bam_aln:
		mdstring = e.optional_field('MD')
		# TODO: evaluate MD using reg exps
		# first understand strand specific read orientation
		if e.iv.strand == '-':
			re_to_use = rev_basesmatch_re
		elif e.iv.strand == '+':
			re_to_use = fwd_basesmatch_re
		else:
			print "[AP]\tERROR:\tStrand NOT assigned to this read", e.read.name, " in BAM alignment file; check it first!!\n\n"
		# execute re finding all occurrences -> it creates an array -> only last/first exact matches. From this number you will classify the read as valid or not
		array_valid_matches = re.findall(re_to_use, mdstring)
		# evaluate returned matches of MD
		if len(array_valid_matches) > 0: # only if there is at least a match -> if not, this means that the read starts with a mismatch
			if int(array_valid_matches[0]) >= min_starting_matches:
				valid = True
				out_array_bam_aln.append(e)
			else:
				print "[AP]\t\t%s -> removed read by MD (Strand: %s, CIGAR: %s ; MD: %s)" %(e.read.name, e.iv.strand, e.cigar[0].type, e.optional_field('MD'))
		else:
			print "[AP]\t\t%s -> removed read by MD (Strand: %s, CIGAR: %s ; MD: %s)" %(e.read.name, e.iv.strand, e.cigar[0].type, e.optional_field('MD'))
	return valid, out_array_bam_aln

def evaluateMDtag(bam_aln, min_starting_matches = 3):
	"""
	Input: pysam alignment object
	Output: true or false
	Based on strand, evaluate possible problems in CIGAR string. Given a as HTSeq BAM reader aln,
	a.optional_field('MD')
	"""
	fwd_basesmatch_re = '\A[0-9]+' # how many correct alignments (bases) do we have at the beginning of the read? -> returns array of max 1 element: {'', 0->}
	rev_basesmatch_re = '[0-9]+\Z' # how many correct alignments (bases) do we have at the end of the read? -> returns array of max 1 element: {'', 0->}
	re_to_use = None # re to be used, one of the previous 2
	valid = False
	
	# find MD tag and string
	mdstring = None
	for tag in bam_aln.tags:
		if tag[0] == 'MD':
			mdstring = str(tag[1])
	# first understand strand specific read orientation
	if bam_aln.is_reverse:
		re_to_use = rev_basesmatch_re
	else:
		re_to_use = fwd_basesmatch_re
	# execute re finding all occurrences -> it creates an array -> only last/first exact matches. From this number you will classify the read as valid or not
	array_valid_matches = re.findall(re_to_use, mdstring)
	# evaluate returned matches of MD
	if len(array_valid_matches) > 0: # only if there is at least a match -> if not, this means that the read starts with a mismatch
		if int(array_valid_matches[0]) >= min_starting_matches:
			valid = True
		else:
			print "[AP]\t\t%s -> removed read by MD" %(bam_aln.qname,)
	else:
		print "[AP]\t\t%s -> removed read by MD" %(bam_aln.qname,)
	return valid


def convertBAMarrayToBEDstring(array_bam_aln, delimiter = "\t", newline = "\n"):
	"""
	Input: array of HTSeq alignments
	Output: string BEd like of the reads
	Notes:
		- BED like:  chr1	10083054	10083104	M00571:14:000000000-A346Y:1:1106:13538:20692/1	25	-
					 chr1	27111566	27111665	M00571:14:000000000-A346Y:1:1104:28699:19778/1	37	-
		- a.aQual, a.cigar[0].type, a.cigar[0].size, a.cigar[0].ref_iv, a.cigar[0].query_from, a.cigar[0].query_to, a.cigar[0].check(), a.iv.strand, a.iv.start, a.iv.start_d, a.read.name, a.read.seq, a.read.qual, a.pe_which
	"""
	outstring = ""
	for aln in array_bam_aln:
		# identify strand and obtain number
		pair_number = 0
		if aln.pe_which == 'second':
			pair_number = 2
		elif aln.pe_which == 'first':
			pair_number = 1
		else:
			print "[AP]\tERROR:\tNo valid pair in BAM read!! id: ", aln.read.name
		# now compose string
		content = delimiter.join(str(x) for x in [aln.iv.chrom, aln.iv.start, aln.iv.end, str(aln.read.name) + "/" + str(pair_number), aln.aQual, aln.iv.strand])
		outstring += content + newline
	
	# return final string
	return outstring

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
	
def getDictionaryFromBED(bed_object_array):
	"""
	Given a BED object, as BEDtools (-> array of bed elements composed by Interval(chrom, start, end, name=".", score=".", strand=".", otherfields=None)), create a dictionary as output.
	Input: BED object
	Output: dictionary of elements: k = read name, v = bed_object
	"""
	outdictionary = {}
	for k in bed_object_array:
		outdictionary[k.name] = k
	if len(outdictionary.keys())==0:
		print "[AP]\tERROR:\tBED file returned empty elements into the dictionary. Why?? Check it please!\n"
	return outdictionary
		

#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	# first check args and file paths
	checkArgs(args)
	minStartingMatches = int(args.minStartingMatches)
	# acquire BED and BAM files: INPUT data
	print "[AP]\tChecked inputs, now acquiring data, both BED and BAM"
	bed = pybedtools.BedTool(args.bedfile) # in bed[index].name you can find the sequence header + /1 or /2 in case of paired reads
	bam = pysam.Samfile(args.bamfile)
	# init output file
	filetarget = open(args.outfilename, 'w')
	
	print "[AP]\tAcquire BED reads (names) into a dictionary"
	bed_dictionary = getDictionaryFromBED(bed)
	print "[AP]\tNow looping over BAM alignments..."
	for bam_alingment in bam:
		#print bam_alingment, bam_alingment.pe_which
		# compose read name according to BED format: name/1 or name/2 based on read pair
		# first, get pair number: 1 or 2?
		pair_number = 0
		if bam_alingment.is_read1:
			pair_number = 1
		elif bam_alingment.is_read2:
			pair_number = 2
		else:
			print "[AP]\tWARNING:\tNo valid pair ID in BAM alignment!! Read id and flag: ", bam_alingment.qname, bam_alingment.flag
		
		# then compose string
		query_aln = bam_alingment.qname.strip() + "/" + str(pair_number)
		#print "query:", query_aln
		# now that you have the query read of the BAM alignment according to BED format, query BED for this read and do what you need
		if query_aln in bed_dictionary.keys():
			valid_cigar = evaluateCigar(bam_alingment) # if M -> ok, else discard it
			valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches) # if MD tag has correct match -> ok, else discard it ## TODO!
			if valid_cigar and valid_MDtag:
				# write data into BED file
				##bed_dictionary[query_aln] # from dictionary, get info and write data
				content = '\t'.join(str(x) for x in [bed_dictionary[query_aln].chrom, bed_dictionary[query_aln].start, bed_dictionary[query_aln].end, bed_dictionary[query_aln].name, bed_dictionary[query_aln].score, bed_dictionary[query_aln].strand])
				filetarget.write(content + "\n")
	
	# finalize file
	filetarget.close()
	
	print "\n[AP]\tTask Finished, closing.\n"

# sentinel
if __name__ == "__main__":
    main()
