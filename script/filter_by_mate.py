#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, os, re, argparse, csv, time
import HTSeq, pybedtools, pysam
from operator import itemgetter, attrgetter
from string import Template
#from Bio import SeqIO
from time import gmtime, strftime
from multiprocessing import Process, Lock, Queue
import multiprocessing, math

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
         ***   FILTER BY MATE   ***
+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		October 2013
 Contact:	andrea.calabria@hsr.it
 Revision:	0.1 (MP)
+--------------------------------------------------+
  
 Note:
  - it uses PySAM and not HTSeq -> HTSeq has an error on PAIR TYPE evaluation (aln.pe_which if often None even if it is first read in pair)!!!!! Whereas PySAM does not have problems (aln.is_read1).
  - BAM file will be indexed (otherwise PySAM will not work)
  
 TODO:
  - 
 
 Logics:
  - if R2 is mapped, is paired, is properly paired:
  	if overlaps with R1:
  		if is a chimera OR the overlap R2 starts before the end of R1 (this is an anomaly! -> insert size will be negative!!):
  		then it is reported in the output file
""" 

description = "This application will filter ISs by Mate pair overlap and SA alternatives."

usage_example = """
Examples of usage:
 (1) 

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--bam', dest="bamfile", help="BAM file to process. No default option.", action="store", required=True)
parser.add_argument('--mateoverlap', dest="mateoverlap", help="Maximum number of non-overlapping bases to consider 2 paired-end reads overlapping. Default = 5 (-> if R1 and R2 overlap leaving at most 5 bp of difference, they are considered fully overlapping).", action="store", required=False, default=5)
parser.add_argument('--maxsamapq', dest="maxsamapq", help="Maximum SA MapQ (mapping quality) to consider this read as chimera due to the alternative. Default = 7 (-> if in the SA a hit has got a MapQ > maxsamapq => this read is a chimera).", action="store", required=False, default=7)
parser.add_argument('--maxinsertsizethreshold', dest="maxinsertsizethreshold", help="Maximum insert size threshold to be considered as valid proper paired reads: if R2 end is within R1 mapping then this mapping is suspicious and this threshold set this limit. Default = 1 (-> if R2 start within R1 mapping +threshold => discard this product).", action="store", required=False, default=1)
parser.add_argument('--sortBAM', dest="sortBAM", action="store_true", help="Sort BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--indexBAM', dest="indexBAM", action="store_true", help="Index BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--pairedEnds', dest="single_end", action="store_true", help="Are your reads single-end and not paired-end? If your reads are single-end, activate this option. Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--SAalnAction', dest="SAalnAction", help="If the BAM file contains SA alignments (Secondary Alignments, that could be also chimera), you can choose which action to do among {ignore, remove, removeByComparison}. This parameter will be used only if the option compareSubOptimal is active/True. Default: ignore.", action="store", default="ignore")
parser.add_argument('--suboptimalThreshold', dest="suboptimalThreshold", help="Specify the suboptimal threshold over which to keep the read [Integer]. If a read has less the best hit alignment score that is different (higher) from the suboptimal alignment score for a percentage suboptimalThreshold then this read is kept, otherwise is discarded. Example: AS=100, XS=80 delta=(100 - XS/AS*100)=20, delta >= suboptimalThreshold? If yes, this read passes filter (is kept). Default: 60 (-> delta>=60 :: valid, else discarded).", action="store", required=False, default=60)
parser.add_argument('--deltaMaxSubOptimal', dest="deltaMaxSubOptimal", help="Threshold in XS/AS between R2 and R1 to accept the R2 Xs alignment as real or not in caso of fully overlapping sequences (R1 and R2). When R2 maps in the same position of R1, then you expect to have the same alignment score (or very similar). If R1 does not have any XS whereas R2 yes, this read is discarded. If delta=(XS/AS|R2)/(XS/AS|R1) >= thisThreshold, the read is discarded, since it means that R2 has got some XS values that R1 does not have even if they are fully overlapping. This threshold  is Integer. Default: 5 (-> delta>=5 -> read discarded).", action="store", required=False, default=5)
parser.add_argument('-o', '--outfilename', dest="outfilename", help="Output file name. The output will be a list of headers to be removed. No default option.", action="store", required=True)
parser.add_argument('-p', '--processes', dest="processes", help="Number of processes to use. Default 2.", action="store", default=2)
parser.add_argument('-r', '--regions', dest="regions", help="Regions to process as CSV string in the format: CHR:START-END. Default value is for HG19: '1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566'.", action="store", default="1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566")
parser.add_argument('-g', '--genome', dest="genome", help="Select the whole genome with relative regions (this is alterative to regions). If you want to load regions, write this as 'regions'. Default: regions. Values: {regions, hg19, mm9, mm10, mfa5, ce, lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada, lv743, lv3029, hiv}.", action="store", default="regions")

args = parser.parse_args()

#########################################################################################
####### GLOBAL VARS
#########################################################################################

# init db values -> from UCSC database ChromInfo
hg19="1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566;M:1-16571"
mm9="1:1-197195432;2:1-181748087;3:1-159599783;4:1-155630120;5:1-152537259;6:1-149517037;7:1-152524553;8:1-131738871;9:1-124076172;10:1-129993255;11:1-121843856;12:1-121257530;13:1-120284312;14:1-125194864;15:1-103494974;16:1-98319150;17:1-95272651;18:1-90772031;19:1-61342430;X:1-166650296;Y:1-15902555;M:1-16299"
mm10="1:1-195471971;2:1-182113224;3:1-160039680;4:1-156508116;5:1-151834684;6:1-149736546;7:1-145441459;8:1-129401213;9:1-124595110;10:1-130694993;11:1-122082543;12:1-120129022;13:1-120421639;14:1-124902244;15:1-104043685;16:1-98207768;17:1-94987271;18:1-90702639;19:1-61431566;X:1-171031299;Y:1-91744698"
# BGI macaca fascicularis genome
bgice="1:1-232298473;2:1-192287759;3:1-199730254;4:1-169826420;5:1-185160406;6:1-180941566;7:1-171133315;8:1-150570238;9:1-133536003;10:1-96524937;11:1-136666246;12:1-108081426;13:1-137686314;14:1-134145161;15:1-108996402;16:1-81259162;17:1-95628244;18:1-75918782;19:1-65364038;20:1-89811811;X:1-155426724;Ur:1-453853"
# NCBI macaca fascicularis genome
mfa5="1:1-227556264;2:1-192460366;3:1-192294377;4:1-170955103;5:1-189454096;6:1-181584905;7:1-171882078;8:1-146850525;9:1-133195287;10:1-96509753;11:1-137757926;12:1-132586672;13:1-111193037;14:1-130733371;15:1-112612857;16:1-80997621;17:1-96864807;18:1-75711847;19:1-59248254;20:1-78541002;X:1-152835861;M:1-16575"
# The following genomes come from TIGET sequences, coded as follows:
# - lv.backbone.fa: --> lv
# 	chr1 = LTR + vecore genome; 
# 	chr2 = dnef
lv="1:1-1747;2:1-82"
# - lv.plasmid.kana.fa --> lvkana
# 	chr1 = LTR
# 	chr2 = vector genome
# 	chr3 = dnef
# 	chr4 = plasmid backbone kana
lvkana="1:1-234;2:1-1513;3:1-82;4:1-3928"
# - lv.plasmid.amp.fa --> lvamp
# 	chr1 = LTR
# 	chr2 = vector genome
# 	chr3 = dnef
# 	chr4 = plasmid backbone amp
lvamp="1:1-234;2:1-1513;3:1-82;4:1-3915"
# - lv.backbone.hpgk.arsa.wprem.fa: --> lvarsa
# 	chr1 = LTR + vecore genome; 
# 	chr2 = dnef
# 	chr3 = hPGK ARSA WPREmut
lvarsa="1:1-1748;2:1-82;3:1-2807"
# - lv.backbone.wasp.was.wprem.fa: --> lvwas
# 	chr1 = LTR + vecore genome; 
# 	chr2 = dnef
# 	chr3 = hWASp WAS Wpre_mut
lvwas="1:1-1748;2:1-82;3:1-3818"
# - transposon.cd19.sb11.fa --> transposon
# 	chr1 = full plasmid sequence
transposon="1:1-6794;2:1-4366"
# - retro.ada.fa --> giada
# 	chr1 = full plasmid sequence
giada="1:1-6185"

## HIV genome
hiv="1:1-9718"

# Custom genome for RNA datasets
# - lv.743.fa
lv743="743.pCCLsin.PPT.hPGK.GFP.Wpre_mut_AMP:1-3968"
# - lv.743.fa
lv3029="3029.pCCL.223.142.3pLTRmir.sin.PPT.BpA.GFP.PGKnoSA.wpre:1-4513"

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
	if args.indexBAM: # index BAM file
		indexBAM(args.bamfile)
	if args.sortBAM: # sort BAM and assign the new file to ARGS and index it
		sorted_bamfile = sortBAM(args.bamfile)
		args.bamfile = sorted_bamfile
		indexBAM(args.bamfile)

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
	valid = False # checking variable to return
	
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
	# only if the MD string is not NULL -> a null string may happen ONLY if the aligner does not give you the MD tag
	if mdstring is not None:
		# execute re finding all occurrences -> it creates an array -> only last/first exact matches. From this number you will classify the read as valid or not
		array_valid_matches = re.findall(re_to_use, mdstring)
		# evaluate returned matches of MD
		if len(array_valid_matches) > 0: # only if there is at least a match -> if not, this means that the read starts with a mismatch
			if int(array_valid_matches[0]) >= min_starting_matches:
				valid = True # return this read as valid
			else:
				print "[AP]\t\t%s -> removed read by MD\t[chr%d:%d-%s, len %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen)
		else:
			print "[AP]\t\t%s -> removed read by MD\t[chr%d:%d-%s, len %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen)
	else: #
		print "[AP]\t\t%s -> WARNING no MD tag found, this read will be kept anyway because I cannot evaluate the MD score\t[chr%d:%d-%s, len %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen)
		valid = True # return this read as valid
	return valid

def evaluateChimeraSA(bam_aln, sa_max_mapq = 7):
	"""
	Input: only the SA tag string (e.g. chr6,142629962,+,66S24M109S,20,0;chr18,54898335,+,152S22M25S,0,0;chr1,20664235,+,141S21M37S,9,0;chr19,9056388,+,46S19M134S,0,0;) , in the format <chr,start,strand,cigar,mapq,nm>
	Output: bool, true if this read is not chimeric
	Logics:
		sa_max_mapq = 10, sa_min_delta = 10. -> if the input read has at least one SA that overcomes mapq (>=10) or delta mapq (<=10), then this read has got a chimera -> not valid
	"""
	ischimera = False
	chimera_counter = 0 # how many mapq overgo the input threshold
	# get SA string, if exists
	SA = None
	for tag, val in bam_aln.tags:
		if tag == 'SA':
			SA = val
	#print "SA:", SA
	# now process SA if exists
	if SA is not None: # is exists, then analyze it
		for x in SA.split(';'):
			if x != '':
				chr, start, strand, cigar, mapq, nm = x.split(',')
				if len(chr.split('_'))<2: # no chimera on random chrs
					if int(mapq) >= sa_max_mapq:
						chimera_counter += 1
						#print "\tchimera found"
		# now assess delta from counter
		if chimera_counter > 0:
			ischimera = True
	else: # is no SA exists, this read is valid by definition
		ischimera = False		
	return ischimera

def getChimeraSAOccurrences(SA_string, ):
	"""
	Input: only the SA tag string (e.g. chr6,142629962,+,66S24M109S,20,0;chr18,54898335,+,152S22M25S,0,0;chr1,20664235,+,141S21M37S,9,0;chr19,9056388,+,46S19M134S,0,0;) , in the format <chr,start,strand,cigar,mapq,nm>
	Output: returns the number of chimeric values from SA string
	Logics:
		sa_max_mapq = 10, sa_min_delta = 10. -> if the input read has at least one SA that overcomes mapq (>=10) or delta mapq (<=10), then this read has got a chimera -> not valid
	"""
	ischimera = False
	chimera_counter = 0 # how many mapq overgo the input threshold
	# now process SA if exists
	if SA_string is not None: # is exists, then analyze it
		for x in SA_string.split(';'):
			if x != '':
				chr, start, strand, cigar, mapq, nm = x.split(',')
				if len(chr.split('_'))<2: # no chimera on random chrs
					chimera_counter += 1
	return chimera_counter

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

def splitChromosomeListByProcess(chr_string):
	"""
	Given chromosome string CSV, split it by number of processes and return list of ToDo.
	"""
	# build an array of chromosomes, sorted by numeric index (X=23, Y=24)
	chrboundaries = []
	chrboundaries.append( int(x) for x in chr_string.strip().split(',') )
	# process todo list: dictionary of chr indexes to process: k = pid, v = (index start, index end) of chroms to compute
	pid_dict = {}
	step = len(chrboundaries) / nproc
	pid = 0
	for i in range(1, len(chrboundaries)+1, step):
		if i+step > len(chrboundaries):
			pid_dict[pid] = (i, len(chrboundaries))
		else:
			pid_dict[pid] = (i, i+step)
		pid += 1

def getRegions(regions_str):
	"""
	Given a string CSV of regions in the form CHR:START-END, return the list of tuples of regions
	"""
	# build an array of regions
	regions_list = [] # array of tuples
	for reg in regions_str.split(';'):
		chr = "chr" + reg.strip().split(':')[0].strip()
		start = int(reg.strip().split(':')[1].strip().split('-')[0])
		end = int(reg.strip().split(':')[1].strip().split('-')[1])
		regions_list.append( (chr, start, end) )
	return regions_list

def getOverlap(bam_aln):
	"""
	Extract overlap of this read with its pair. Values: [0-100]
	Overlap is defined as:
		abs(insert size) / abs(position start - position end)
	From BAM and PySAM, start and end are independent from the orientation, as well as the mate position. These are the fields: b.positions[0], b.positions[-1], b.pnext, b.rnext, b.mate_is_reverse, b.tlen,
	"""
	overlap = None # overlap percentage
	start = bam_aln.positions[0]
	end = bam_aln.positions[-1]
	mapped_len = abs(int(start)-int(end))
	insert_size = float(bam_aln.tlen)
	overlap = float(mapped_len)/abs(insert_size)
	if overlap <= 0.0 or overlap > 100.0:
		print "[AP]\t\tERROR!!\t%s Negative or >100%% overlap! Why?!?! Check it please!" %(bam_aln.qname, )
		sys.exit()
	return overlap

def getCigarTagValFirstOccurrenceByStrand(bam_aln, usertag_string):
	"""
	Extract the value for the specific input tag in the CIGAR string.
	Given the Pysam codes for CIGAR string conversion.
	"""
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
	usertag = cigar_flag[usertag_string]
	tagval = None
	# based on the orientation, get the first tag value
	if not bam_aln.is_reverse: # if 5' - 3', start from first occurrence of M
		seen = False
		for tag, val in bam_aln.cigar:
			if tag == usertag and not seen:
				tagval = val
				seen = True
	else: # reverse case, take the last occurrence -> revert list, 
		seen = False
		for tag, val in bam_aln.cigar[::-1]:
			if tag == usertag and not seen:
				tagval = val
				seen = True
	return tagval # can be none!

def evaluateMateOverlap(bam_aln, Mval, maxNonOverlappingBases = 5):
	"""
	Mval can be None!!
	"""
	overlapping = False # return var
	# vars from bam alignment
	read_start = bam_aln.positions[0] # genome start (5')
	mate_start = bam_aln.pnext # starting position of the mate
	# compute the formula to assess overlap
	if Mval is None: # if no M val available, then the formula is simplified
		D = abs(bam_aln.tlen)
	else:
		D = abs( abs(bam_aln.tlen) - abs(int(Mval)) )
	# mate overlap value
	MateOv = abs( read_start - mate_start ) - D
	if MateOv <= maxNonOverlappingBases:
		# the read is considered overlapping the mate
		overlapping = True
	return overlapping

def evaluateR2startsWithinR1(bam_aln, r1_dictionary, r1_end_threshold = 1):
	"""
	Given the BAM alignment of R2, check if the insert size len is compatible with the R2 start -> R2 start MUST not start within the R1 mapped read
	"""
	r2_starts_withinr1 = False # in default this read is assumed correct, thus it DOES not start within R1 mapping (before the end of R1)
	r1_end = None
	r2_start = None
	if r1_dictionary.has_key(bam_aln.qname): # if not, no IS will be reported as well -> go ahead
		if bam_aln.is_reverse: # rev/- strand
			r1_end = r1_dictionary[bam_aln.qname]['end']
			r2_start = bam_aln.positions[-1]
			if r2_start < (r1_end + r1_end_threshold):
				r2_starts_withinr1 = True
		else: # fwd/+ strand
			r1_end = r1_dictionary[bam_aln.qname]['end']
			r2_start = bam_aln.positions[0]
			if r2_start > (r1_end - r1_end_threshold):
				r2_starts_withinr1 = True
	return r2_starts_withinr1
	
def evaluateR2startsWithinR1_fromDictionary(reads_dictionary, ):
	"""
	Returns a set from an input that is a disciontary of pairs, r1 and r2, with 'start', 'end' and 'strand' keys.
	"""
	set_todiscard = set()
	for qname, pairs in reads_dictionary.iteritems():
		if pairs.has_key('r2') and pairs.has_key('r1'):
			# which strand?
			if (pairs['r1']['strand'] == '+' and pairs['r2']['start'] < pairs['r1']['end']) or (pairs['r1']['strand'] == '-' and pairs['r2']['start'] > pairs['r1']['end']):
				set_todiscard.add(qname)
				#print "Discard by R2 start: R1 strand, R2 start, R1 end", pairs['r1']['strand'], pairs['r2']['start'], pairs['r1']['end']
	return set_todiscard

# def evaluateFullOverlap(reads_dictionary):
# 	"""
# 	If R2 is fully overlapping R1 (start<->end) then put this read in the output set.
# 	This flag will be useful for further analysis, for example for XS/AS ratio.
# 	"""
# 	outset = set()
# 	for qname, pairs in reads_dictionary.iteritems():
# 		if pairs.has_key('r2') and pairs.has_key('r1'): # only paired end reads
# 			# which strand?
# 			if 	(pairs['r1']['strand'] == '+' and pairs['r2']['start'] >= pairs['r1']['end'] and pairs['r2']['end'] == pairs['r1']['start'] ) 
# 				or 
# 				(pairs['r1']['strand'] == '-' and pairs['r2']['start'] > pairs['r1']['end']):
# 				set_todiscard.add(qname)
# 				#print "Discard by R2 start: R1 strand, R2 start, R1 end", pairs['r1']['strand'], pairs['r2']['start'], pairs['r1']['end']
# 	return outset

def evaluateR2boundaries(reads_dictionary, suboptimalThreshold, deltaMaxSubOptimal):
	"""
	Output: set of reads that do not satisfy our requirements (strand dependent):
		0. [properly paired] if strands are identical, discard it
		1. [correct insert size, >=0, where 0 = full overlap, identical pairs len] R2 must not start within R1 mapping read -> else discard
		2. [correct insert size, >=0] R2 end must not go over the IS, that is R1 start -> else discard
		3. [full ovelap] if R2 start = R1 end => R2 end = R1 start, else discard the read
		4. [full ovelap] if R2 end == R1 start (-> IS) then R2 end must be at least >= R1 end => full overlap -> if ok, then:
			A. if R1 does not have XS but R2 yes => discard it (WARNING: if R2 start is out of the R1 end, thus longer than R1, this read could map elsewhere too)
			B. if R2 has got XS as well as R1, then compare XS/AS ratio between R2 and R1. This comparison must be in line, that is ratioR2 must not be over a threshold (e.g. 5) than ratioR1
			B1. if B is true (ration in the valid range and you do not need to discard the reed), then count the number of SA: if R2 has a number of SA that is higher than R1 => discard the product
			C. evaluate R2 AlnScore alone (XS/SA ratio) -> discard if ratio is bad (over threshold, i.e RatioR2 < 60)
			D. R2 must not start with cigar Soft Clipped
	"""
	bp_extension = 3
	set_todiscard = set()
	for qname, pairs in reads_dictionary.iteritems():
		if pairs.has_key('r2') and pairs.has_key('r1'):
			# 0. [properly paired] if strands are identical, discard it
			if (pairs['r1']['strand'] == pairs['r2']['strand']):
				set_todiscard.add(qname)
				print "[AP]\t\t%s -> removed read by in-proper pairing\t[R1 <chr%s'%s',s:%d;e:%d> R2 <chr%s,'%s',s:%d;e:%d>]" %(qname, pairs['r1']['chr'], pairs['r1']['strand'], pairs['r1']['start'], pairs['r1']['end'], pairs['r2']['chr'], pairs['r2']['strand'], pairs['r2']['start'], pairs['r2']['end'])
			
			# 1. [correct insert size, >=0] R2 must not start within R1 mapping read -> else discard
			if (pairs['r1']['strand'] == '+' and pairs['r2']['start'] < pairs['r1']['end']) or (pairs['r1']['strand'] == '-' and pairs['r2']['start'] > pairs['r1']['end']):
				set_todiscard.add(qname)
				print "[AP]\t\t%s -> removed read by invalid insert size (R2 starts within R1 span)\t[Prd Strand: %s; R1 <s:%d;e:%d> R2 <s:%d;e:%d>]" %(qname, pairs['r1']['strand'], pairs['r1']['start'], pairs['r1']['end'], pairs['r2']['start'], pairs['r2']['end'])
			
			# 2. [correct insert size, >=0] R2 end must not go over the IS, that is R1 start -> else discard
			if (pairs['r1']['strand'] == '+' and pairs['r2']['end'] < (pairs['r1']['start']-bp_extension)) or (pairs['r1']['strand'] == '-' and pairs['r2']['end'] > (pairs['r1']['start']+bp_extension)):
				set_todiscard.add(qname)
				print "[AP]\t\t%s -> removed read by invalid insert size (R2 ends over IS_R1 start)\t[Prd Strand: %s; R1 <chr%s,s:%d;e:%d> R2 <chr%s,s:%d;e:%d>]" %(qname, pairs['r1']['strand'], pairs['r1']['chr'], pairs['r1']['start'], pairs['r1']['end'], pairs['r2']['chr'], pairs['r2']['start'], pairs['r2']['end'])
			
			## 3. [full ovelap] if R2 start = R1 end => R2 end = R1 start, else discard the read. [ITA note: sembra molto stringente, ma se ci pensi ci sta: se con R1 vedo la read mappata, con R2 devo avere la stessa a meno di errori di trimming. questi errori lo posso considerare veri errori e quindi la elimino, qundi questo filtro cautela anche da sbavature di trimming LC/LTR]
			#if pairs['r2']['start'] == pairs['r1']['end']:
			#	if pairs['r2']['end'] <> pairs['r1']['start']:
			#		set_todiscard.add(qname)
			#		print "[AP]\t\t%s -> removed read by full overlap position (R2s=R1e & R1s<>R2e)\t[R1 <s:%d;e:%d> R2 <s:%d;e:%d>]" %(qname, pairs['r1']['start'], pairs['r1']['end'], pairs['r2']['start'], pairs['r2']['end'])
			
			## 3. [full ovelap] if R2 start = R1 end => R2 must not start with Soft Clipped cigar score and must not end with soft clipped cigar score
			if pairs['r2']['start'] == pairs['r1']['end']:
				if pairs['r2']['cigar'] is not None:
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
					# independently by the strand, do not start with cigar S
					tag, val = pairs['r2']['cigar'][0]
					if tag == 4 and val > bp_extension:
						set_todiscard.add(qname)
						print "[AP]\t\t%s -> removed read by CIGAR soft clip start (R2 start = R1 end, R2 CIGAR start S)\t[R2 <'%s',s:%d;e:%d>, cigar=%s; R1 chr%s, R2 chr%s]" %(qname, pairs['r2']['strand'], pairs['r2']['start'], pairs['r2']['end'], pairs['r2']['cigar'], pairs['r1']['chr'], pairs['r2']['chr'])
					tag, val = pairs['r2']['cigar'][-1]
					if tag == 4 and val > bp_extension:
						set_todiscard.add(qname)
						print "[AP]\t\t%s -> removed read by CIGAR soft clip end (R2 start = R1 end, R2 CIGAR end S)\t[R2 <'%s',s:%d;e:%d>, cigar=%s; R1 chr%s, R2 chr%s]" %(qname, pairs['r2']['strand'], pairs['r2']['start'], pairs['r2']['end'], pairs['r2']['cigar'], pairs['r1']['chr'], pairs['r2']['chr'])
			
			# 4. [full ovelap] if R2 end == R1 start (-> IS) then R2 end must be at least >= R1 end => full overlap -> this must return also the corresponding list that will be passed to the AlnScore (XS/SA ratio)
			# {'r1': {'start': 166650228L, 'AS': 105, 'end': 166650124L, 'strand': '-', 'XS': 69}, 'r2': {'start': 166650124L, 'AS': 62, 'end': 166650200L, 'strand': '+', 'XS': 65}}}
			if (pairs['r1']['strand'] == '+' and pairs['r2']['end'] == pairs['r1']['start'] and pairs['r2']['start'] == pairs['r1']['end']) or (pairs['r1']['strand'] == '-' and pairs['r2']['end'] == pairs['r1']['start'] and pairs['r2']['start'] == pairs['r1']['end']): # you can also avoid to do the strand specific analysis because it is redundant, but just in case of adding a threshold/tolerance everything is ready and easier to integrate
				# A. if R1 does not have XS but R2 yes => discard it (WARNING: if R2 start is out of the R1 end, thus longer than R1, this read could map elsewhere too)
				if pairs['r2']['XS'] is not None and pairs['r1']['XS'] is None:
					set_todiscard.add(qname)
					print "[AP]\t\t%s -> removed read by full overlap position by XS (XS|R2 not None, XS|R1 None)\t[R1 <chr%s,AS=%d;XS=%s> R2 <chr%s,AS=%d;XS=%s>]" %(qname, pairs['r1']['chr'], pairs['r1']['AS'], pairs['r1']['XS'], pairs['r2']['chr'], pairs['r2']['AS'], pairs['r2']['XS'])
				# B. if R2 has got XS as well as R1, then compare XS/AS ratio between R2 and R1. This comparison (it is a difference) must be in line, that is ratioR2 must not be over a threshold (e.g. 5) than ratioR1
				ratio_aln_R1 = None
				ratio_aln_R2 = None
				if pairs['r2']['XS'] is not None and pairs['r1']['XS'] is not None and pairs['r1']['AS'] > 0 and pairs['r2']['AS'] > 0:
					ratio_aln_R1 = pairs['r1']['XS']/float(pairs['r1']['AS'])*100
					ratio_aln_R2 = pairs['r2']['XS']/float(pairs['r2']['AS'])*100
					if (ratio_aln_R2-ratio_aln_R1) > deltaMaxSubOptimal:
						set_todiscard.add(qname)
						print "[AP]\t\t%s -> removed read by full overlap position by ratio XS/AS|R1/R2\t[R1 <chr%s,AS=%d;XS=%s;ratio=%.2f> R2 <chr%s,AS=%d;XS=%s;ratio=%.2f>]" %(qname, pairs['r1']['chr'], pairs['r1']['AS'], pairs['r1']['XS'], ratio_aln_R1, pairs['r2']['chr'], pairs['r2']['AS'], pairs['r2']['XS'], ratio_aln_R2)
					else:
						# B1. if B is true (ration in the valid range and you do not need to discard the reed), then count the number of SA: if R2 has a number of SA that is higher than R1 => discard the product
						R1_SA_number = getChimeraSAOccurrences(pairs['r1']['SA'])
						R2_SA_number = getChimeraSAOccurrences(pairs['r2']['SA'])
						if R2_SA_number > R1_SA_number:
							#print R1_SA_number, R2_SA_number
							set_todiscard.add(qname)
							print "[AP]\t\t%s -> removed read by full overlap position by number of SA (SA|R2>SA|R1)\t[R1 <chr%s,AS=%d;XS=%s;SAn=%d> R2 <chr%s,AS=%d;XS=%s;SAn=%d>]" %(qname, pairs['r1']['chr'], pairs['r1']['AS'], pairs['r1']['XS'], R1_SA_number, pairs['r2']['chr'], pairs['r2']['AS'], pairs['r2']['XS'], R2_SA_number)
				# C. evaluate R2 AlnScore alone (XS/SA ratio) -> discard if ratio is bad (over threshold, i.e RatioR2 < 60)
				# evaluate XS/AS
				AS = None
				XS = None
				if pairs['r2'].has_key('XS'):
					XS = pairs['r2']['XS']
				if pairs['r2'].has_key('AS'):
					AS = int(pairs['r2']['AS'])
				if XS is not None:
					if AS > XS: # only if the best hit alignment score is better than subotimal one, go ahead, else this read is not valid
						# find delta: as - (xs/as*100)
						delta = 100-(XS/float(AS)*100)
						# compare delta with threshold
						if delta <= float(suboptimalThreshold):
							set_todiscard.add(qname)
							print "[AP]\t\t%s -> removed read by full overlap XS/AS score (R2 delta>threshold)\t[R2 <s:%d;e:%d>, AS=%d, XS=%d, delta=%.2f, threshold=%s; R1 chr%s, R2 chr%s]" %(qname, pairs['r2']['start'], pairs['r2']['end'], pairs['r2']['AS'], pairs['r2']['XS'], delta, suboptimalThreshold, pairs['r1']['chr'], pairs['r2']['chr'])
					else: # in the case in which AS <= XS, exclude the read a priori
						set_todiscard.add(qname)
						print "[AP]\t\t%s -> removed read by full overlap (R2 AS<=XS)\t[R2 <s:%d;e:%d>, AS=%d, XS=%d]" %(qname, pairs['r2']['start'], pairs['r2']['end'], pairs['r2']['AS'], pairs['r2']['XS'])
				# D. R2 must not start with cigar Soft Clipped
				if pairs['r2']['cigar'] is not None:
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
						if pairs['r2']['strand'] == '+':
							tag, val = pairs['r2']['cigar'][0]
							if tag == 4 and val > bp_extension:
								set_todiscard.add(qname)
								print "[AP]\t\t%s -> removed read by full overlap CIGAR soft clip start\t[R2 <'%s',s:%d;e:%d>, cigar=%s; R1 chr%s, R2 chr%s]" %(qname, pairs['r2']['strand'], pairs['r2']['start'], pairs['r2']['end'], pairs['r2']['cigar'], pairs['r1']['chr'], pairs['r2']['chr'])
						if pairs['r2']['strand'] == '-':
							tag, val = pairs['r2']['cigar'][-1]
							if tag == 4 and val > bp_extension:
								set_todiscard.add(qname)
								print "[AP]\t\t%s -> removed read by full overlap CIGAR soft clip start\t[R2 <'%s',s:%d;e:%d>, cigar=%s; R1 chr%s, R2 chr%s]" %(qname, pairs['r2']['strand'], pairs['r2']['start'], pairs['r2']['end'], pairs['r2']['cigar'], pairs['r1']['chr'], pairs['r2']['chr'])
						
			#else: # if R2 end == R1 start then R2 start MUST be at least R1 end -> else discard it (this is the else branch) [same as the rule 3]
			#	set_todiscard.add(qname)
			#	print "[AP]\t\t%s -> removed read by full overlap position (R2s=R1e & R1s<>R2e)\t[R1 <s:%d;e:%d> R2 <s:%d;e:%d>]" %(qname, pairs['r1']['start'], pairs['r1']['end'], pairs['r2']['start'], pairs['r2']['end'])
			# debug
			if qname == 'M00571:39:000000000-A4EBU:1:2111:23888:14795' or qname == 'M00571:39:000000000-A4EBU:1:2107:16785:4759' or qname == 'M00571:39:000000000-A4EBU:1:2106:4448:15143' or qname == 'M00571:39:000000000-A4EBU:1:2107:10572:27095' or qname == 'M00571:39:000000000-A4EBU:1:2107:24180:19334':
				print "\n\n", qname, reads_dictionary[qname], "\n\n"
	return set_todiscard

def evaluateR2fullOverlapsAndR2endsBeforeR1start_fromDictionary():
	"""
	1. R2 end non puo andare oltre R1 start.
	2. se R2 end parte da R1 end, allora R2 end DEVE essere R1 start (ovvero full overlap)
	
	Se 1 o 2 sono false allora la read e' da scartare.
	"""
	pass

def worker(todo_regions, bam, maxbptocondidermateoverlap, maxsamapq, maxinsertsizethreshold, suboptimalThreshold, deltaMaxSubOptimal, out_q, SAalnAction): 
	""" 
	The worker function, invoked in a process. 
		region: array of tuples with coordinate to slice chromosome
		bam: the input bam to splice
		out_q: queue of the process
		mateoverlap_perc: percentage of mate overlap between R1 and R2
	"""
	outdata = []
	for region in todo_regions:
		chr, start, end = region
		r1 = {} # k=qname of R1 (you will have the same id for R2!!!), v='start', 'end' based on the orientation!!!
		product = {} # k = header, v = 'r1': {'start', 'end'}, 'r2': {'start', 'end'}
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t...processing region (chr, start, end)", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		
		print "[AP]\t\t\tAcquiring dictionary of reads in the region (chr, start, enc)\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		r1_nonvalid_byfirstif = set()
		for bam_alignment in mybam: # for each alignment in the BAM file (both R1 and R2!!!)
			# acquire data structure of mapped pp R1 of all reads. this data will be useful for understanding if the START of R2 is within the mapped R1. If yes, this will be removed!
			if bam_alignment.is_paired and bam_alignment.is_proper_pair and bam_alignment.is_read1 and not bam_alignment.is_secondary and not  bam_alignment.is_unmapped:
				AS = None # best hit alignment score
				XS = None # suboptimal alignment score
				SA = None
				for tag, val in bam_alignment.tags: # acquire data/values from tags
					if tag == 'AS': # best hit alignment AS score
						AS = val
					if tag == 'XS': # XS tag, suboptimal alignment score
						XS = val
					if tag == 'SA': # XS tag, suboptimal alignment score
						SA = val
				if bam_alignment.is_reverse: # bam_alignment.positions[0], bam_alignment.positions[-1]
					r1[bam_alignment.qname] = {	'start': bam_alignment.positions[-1], 'end': bam_alignment.positions[0]}
					if product.has_key(bam_alignment.qname):
						product[bam_alignment.qname]['r1'] = {		'start': bam_alignment.positions[-1], 
																	'end': bam_alignment.positions[0], 
																	'strand': '-' , 
																	'AS': AS,
																	'XS': XS,
																	'SA': SA,
																	'chr': bam_alignment.tid, 
																	'cigar': bam_alignment.cigar[::-1],
															}
					else:
						product[bam_alignment.qname] = {	'r1': {		'start': bam_alignment.positions[-1], 
																	'end': bam_alignment.positions[0], 
																	'strand': '-' , 
																	'AS': AS,
																	'XS': XS,
																	'SA': SA,
																	'chr': bam_alignment.tid, 
																	'cigar': bam_alignment.cigar[::-1],
															}
													}
				else:
					r1[bam_alignment.qname] = {'start': bam_alignment.positions[0], 'end': bam_alignment.positions[-1]}
					if product.has_key(bam_alignment.qname):
						product[bam_alignment.qname]['r1'] = {	'start': bam_alignment.positions[0], 
																'end': bam_alignment.positions[-1], 
																'strand': '+',
																'AS': AS,
																'XS': XS,
																'SA': SA,
																'chr': bam_alignment.tid,
																'cigar': bam_alignment.cigar,
																}
					else:
						product[bam_alignment.qname] = {	'r1': {	'start': bam_alignment.positions[0], 
																'end': bam_alignment.positions[-1], 
																'strand': '+',
																'AS': AS,
																'XS': XS,
																'SA': SA,
																'chr': bam_alignment.tid,
																'cigar': bam_alignment.cigar,
																}
													}
			
			# now R2
			if bam_alignment.is_paired and bam_alignment.is_proper_pair and bam_alignment.is_read2 and not bam_alignment.is_secondary and not  bam_alignment.is_unmapped:
				AS = None # best hit alignment score
				XS = None # suboptimal alignment score
				SA = None
				for tag, val in bam_alignment.tags: # acquire data/values from tags
					if tag == 'AS': # best hit alignment AS score
						AS = val
					if tag == 'XS': # XS tag, suboptimal alignment score
						XS = val
					if tag == 'SA': # XS tag, suboptimal alignment score
						SA = val
				#print "trovato R2"
				if bam_alignment.is_reverse: # bam_alignment.positions[0], bam_alignment.positions[-1]
					#print "reverse"
					if product.has_key(bam_alignment.qname):
						product[bam_alignment.qname]['r2'] = {	'start': bam_alignment.positions[-1], 
																'end': bam_alignment.positions[0], 
																'strand': '-',
																'AS': AS,
																'XS': XS,
																'SA': SA,
																'chr': bam_alignment.tid,
																'cigar': bam_alignment.cigar[::-1],
																} 
						#print bam_alignment.qname
					else:
						product[bam_alignment.qname] = {	'r2': {	'start': bam_alignment.positions[-1], 
																'end': bam_alignment.positions[0], 
																'strand': '-',
																'AS': AS,
																'XS': XS,
																'SA': SA,
																'chr': bam_alignment.tid,
																'cigar': bam_alignment.cigar[::-1],
																}
													}
				else:
					if product.has_key(bam_alignment.qname):
						product[bam_alignment.qname]['r2'] = {	'start': bam_alignment.positions[0], 
																'end': bam_alignment.positions[-1], 
																'strand': '+',
																'AS': AS,
																'XS': XS,
																'SA': SA,
																'chr': bam_alignment.tid,
																'cigar': bam_alignment.cigar,
																}
					else:
						product[bam_alignment.qname] = {	'r2': {	'start': bam_alignment.positions[0], 
																'end': bam_alignment.positions[-1], 
																'strand': '+',
																'AS': AS,
																'XS': XS,
																'SA': SA,
																'chr': bam_alignment.tid,
																'cigar': bam_alignment.cigar,
																}
														}
					#else:
						#print "[AP]\t\t\tNot paired read (given R2)", bam_alignment.qname
			## debug
			#if bam_alignment.qname == 'M00571:39:000000000-A4EBU:1:2107:24180:19334':
			#	if bam_alignment.is_paired and bam_alignment.is_proper_pair and bam_alignment.is_read2 and not bam_alignment.is_secondary:
			#		print "R2", bam_alignment

		
		# from the dictionary, get the reads to exclude a priori due to invalid insert size: R2 starts within the R1 product/mapping
		#print product
		#tobe_discarded_by_insertsize = evaluateR2startsWithinR1_fromDictionary(product)
		todiscard = evaluateR2boundaries(product, suboptimalThreshold, deltaMaxSubOptimal)
		todiscard = todiscard | r1_nonvalid_byfirstif
		# now loop over the alignments
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t\tProcessing R2 reads in the region (chr, start, enc)\t\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		for bam_alignment in mybam: # for each alignment in the BAM file (both R1 and R2!!!)
			# process only R2 now evaluating mate filtering
			if bam_alignment.is_paired and  bam_alignment.is_proper_pair and  bam_alignment.is_read2 and not  bam_alignment.is_secondary and not  bam_alignment.is_unmapped: # if read is paired, properly paired, is R2, not secondary and mapped -> process it, else DO NOT report it to be removed
				### obsolete::: starts_withinR1 = evaluateR2startsWithinR1( bam_alignment, r1, r1_end_threshold = maxinsertsizethreshold) # does R2 starts within R1 mapped region? true/false
				### now use dictionary instead: the reason is that we need to have this check on a global var, not in the for loop.
				if bam_alignment.qname in todiscard:
					outdata.append(bam_alignment.qname)
				else: # if not, then starts with the rest: evaluate chimera due to R2 SA
					# check if SA tag is available for R2 and if yes, process this read R2
					SA = None
					Mtag = 0
					for tag, val in  bam_alignment.tags:
						if tag == 'SA':
							SA = val
					if SA is not None: # if SA exists, go ahead
						Mval = getCigarTagValFirstOccurrenceByStrand( bam_alignment, 'M') # if exists, M value in the cigar is collected here, orientation-based. can be None!!
						is_mateoverlap = evaluateMateOverlap( bam_alignment, Mval, maxNonOverlappingBases = maxbptocondidermateoverlap) # is the mate (R1) overlapping this read (R2)?
						if is_mateoverlap: # if it is overlapping, then tell me if it is also a chimera. if yes, then THIS pair will be removed, else NO.
							is_chimera = evaluateChimeraSA( bam_alignment, sa_max_mapq = maxsamapq ) # is this read a chimera? (=> that is it has got other SA above the given MAPQ)
							if is_chimera: # then this pair will be removed
								outdata.append(bam_alignment.qname)
		print "[AP]\t\t\tDone, region completed (chr, start, end)\t\t\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" 
																
	# queue results at the end of this process
	out_q.put(outdata)


#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	start_time = time.time() # get read of starting time
	print "[AP]\tYour settings:"
	for k, v in vars(args).iteritems():
		print "\t\t%s::\t\t%s" %(k, v)
	# first check args and file paths
	print "[AP]\tChecking inputs."
	checkArgs(args)
	# acquire thresholds
	maxbptocondidermateoverlap = float(args.mateoverlap)
	maxsamapq = int(args.maxsamapq)
	maxinsertsizethreshold = int(args.maxinsertsizethreshold)
	suboptimalThreshold = int(args.suboptimalThreshold)
	deltaMaxSubOptimal = int(args.deltaMaxSubOptimal)
	print "[AP]\tIndexing BAM (this is required if you did not do it)."
	indexing = pysam.index(args.bamfile) # indexing is BLANK
	# acquire BED and BAM files: INPUT data
	print "[AP]\tAcquiring data, both BED and BAM."
	bam = pysam.Samfile(args.bamfile) # the bam file to process
	print "[AP]\tAcquiring region list to use while slicing BAM."
	regions = None
	if args.genome is not "regions":
		if args.genome == "hg19":
			regions = getRegions(hg19)
		elif args.genome == "mm9":
			regions = getRegions(mm9)
		elif args.genome == "mm10":
			regions = getRegions(mm10)
		elif args.genome == "mfa5":
			regions = getRegions(mfa5)
		elif args.genome == "ce":
			regions = getRegions(bgice)
		elif args.genome == "lv":
			regions = getRegions(lv)
		elif args.genome == "lvwas":
			regions = getRegions(lvwas)
		elif args.genome == "lvarsa":
			regions = getRegions(lvarsa)
		elif args.genome == "lvamp":
			regions = getRegions(lvamp)
		elif args.genome == "lvkana":
			regions = getRegions(lvkana)
		elif args.genome == "transposon":
			regions = getRegions(transposon)
		elif args.genome == "giada":
			regions = getRegions(giada)
		elif args.genome == "lv743":
			regions = getRegions(lv743)
		elif args.genome == "lv3029":
			regions = getRegions(lv3029)
		elif args.genome == "hiv":
			regions = getRegions(hiv)
		else:
			print "\t\tLoading regions error! Please read the help (-h)."
			sys.exit()
	else:
		regions = getRegions(args.regions)
	print "[AP]\tNow looping over BAM alignments Splitting run in processes..."
	out_q = Queue() # init process queue
	nprocs = int(args.processes)
	chunksize = int(math.ceil(len(regions) / float(nprocs))) # find the chunk wrt processes for running regions
	procs = [] # array of processes
	# run processes in a range of chunks
	for i in range(nprocs):
		p = multiprocessing.Process(target=worker, args=(regions[chunksize * i:chunksize * (i + 1)], bam, maxbptocondidermateoverlap, maxsamapq, maxinsertsizethreshold, suboptimalThreshold, deltaMaxSubOptimal, out_q, args.SAalnAction )) # prototype: todo_regions, bam, bed_dictionary, minStartingMatches, out_q
		procs.append(p)
		p.start()
	# Collect all results into a single result list.
	print "[AP]\tGetting results from different processes."
	resultdata = [] # this is what will be written into the file
	for i in range(nprocs):
		resultdata += out_q.get()
	# Wait for all worker processes to finish
	print "[AP]\tWaiting that all processes end."
	for p in procs:
		p.join()
	# write output results from joined array of strings
	print "[AP]\tWriting results into output file", args.outfilename
	filetarget = open(args.outfilename, 'w') # init output file
	for row in resultdata:
		filetarget.write(row + "\n") # each row already contains \n
	filetarget.close() # finalize file
	
	elapsed_time = time.time() - start_time # get read of finishing time
	print "\n[AP]\tTask Finished, closing.\n\tElapsed time: %.2f [seconds]\n" %(elapsed_time)
	

# sentinel
if __name__ == "__main__":
    main()
