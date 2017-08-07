#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, os, re, argparse, csv, time, sqlite3
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
        ***   FILTER BY CIGAR   ***
+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		May 2013 - October 2013
 Contact:	andrea.calabria@hsr.it
 Revision:	0.4 (MP, platform independent, indels, aln score, SA/chimeras)
 Note:		inverse logics of previous version!
 			now write in file the reads to prune!
+--------------------------------------------------+
  
 Note:
  - Min threshold for MATCHES of MD tag: 3 (<=3 -> discard read)
  - it uses PySAM and not HTSeq -> HTSeq has an error on PAIR TYPE evaluation (aln.pe_which if often None even if it is first read in pair)!!!!! Whereas PySAM does not have problems (aln.is_read1).
  - BAM file will be indexed (otherwise PySAM will not work)
  
 TODO:
  - 

 Steps
	My filter:: depending on strand, 
	1. check CIGAR: it must start/end with at least 4M else discard read
	2. check MD:
		if MD is digit -> ok: tutti buoni
		else, check the first/end number is <=4 -> discard read
	3. check INDELS in first bases
	4. check alignment scores (best vs suboptimal)
	5. check SA/Chimeras
""" 

description = "This application will filter ISs by CIGAR score and MD tag."

usage_example = """
Examples of usage:
 (1) Illumina paired-end reads
 (2) 454 reads

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--bam1', dest="bamfile1", help="BAM file 1 to process. No default option.", action="store", required=True)
parser.add_argument('--bam2', dest="bamfile2", help="BAM file 2 to process. No default option.", action="store", required=True)
parser.add_argument('-g1', '--genome1', dest="genome1", help="Select the whole genome 1 with relative regions (this is alterative to regions). If you want to load regions, write this as 'regions'. Default: regions. Values: {regions, hg19, mm9, mfa5, ce, lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada, lv743, lv3029}.", action="store", default="regions")
parser.add_argument('-g2', '--genome2', dest="genome2", help="Select the whole genome 2 with relative regions (this is alterative to regions). If you want to load regions, write this as 'regions'. Default: regions. Values: {regions, hg19, mm9, mfa5, ce, lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada, lv743, lv3029}.", action="store", default="regions")
parser.add_argument('--sortBAM', dest="sortBAM", action="store_true", help="Sort BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--indexBAM', dest="indexBAM", action="store_true", help="Index BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)

parser.add_argument('--deltaAlignment', dest="deltaAlignment", help="Max threshold for alignment overlap of the same read between the two input genomoes. This is useful when a read is mapped in both genomes but with some conflict bases. Example: given a read of 100 bp, it maps 80M20S on HG19 and 77S23M on LV. This means that 3 bp are in common between the two genoms. This option is aimed at setting up the number of maximum accepted bases in overla. Default = 5 (-> <=5 are valid reads, else discard the read from both genomes).", action="store", required=False, default=3)
parser.add_argument('--minStartingMatches', dest="minStartingMatches", help="Minimum starting matches to consider the read valid [Integer]. If a read has less than minStartingMatches, it will be discarded. Default = 3 (-> >=3 are valid, else discarded).", action="store", required=False, default=3)
parser.add_argument('--minStartingBasesNoIndels', dest="minStartingBasesNoIndels", help="Minimum number of starting bases in which no insertions nor deletions must occur (orientation dependent). This number contains minStartingMatches. Write an Integer value. If a read has an indel within minStartingBasesNoIndels bass, it will be discarded. Default = 7 (-> from 1 to 7 no indels. <=7 are valid, else discarded).", action="store", required=False, default=7)
parser.add_argument('--endClipThreshold', dest="endClipThreshold", help="Maximum number of bases for which R1 can end without M matches. This is typically useful for soft clip bases: if the end of R1 has a number of soft clipped bases greater than threshold, this read is removed. Default = 10 (-> >10 are filtered out).", action="store", required=False, default=10)

parser.add_argument('--singleEnd', dest="single_end", action="store_true", help="Are your reads single-end and not paired-end? If your reads are single-end, activate this option. Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--pruneByPair', dest="pruneByPair", action="store", help="If you have a specific pair that drives the filter (e.g. read1 or read2), specify it here {both, read1, read2}. This is useful if your experimental design attemps to analyze at the end only one of the two pairs, for example all integration sites come from read1 only. If your reads are single-end this option is not activated. Default: 'both'.", default='both')
parser.add_argument('--compareSubOptimal', dest="compareSubOptimal", action="store_true", help="Do you want to compare the alignment score versus the suboptimal alignment score? If yes, here you are going to compare the AS tag (that you can change in the ASlikeTag option) and the XS tag (set it in the XSlikeTag option), as output from BWA MEM. Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--suboptimalThreshold', dest="suboptimalThreshold", help="Specify the suboptimal threshold over which to keep the read [Integer]. If a read has less the best hit alignment score that is different (higher) from the suboptimal alignment score for a percentage suboptimalThreshold then this read is kept, otherwise is discarded. Example: AS=100, XS=80 delta=(AS - XS/AS*100)=20, delta >= suboptimalThreshold? If yes, this read passes filter (is kept). Default: 10 (-> delta>=10 :: valid, else discarded).", action="store", required=False, default=10)
parser.add_argument('--ASlikeTag', dest="ASlikeTag", help="If the BAM file contains an alignment score tag different from 'AS', change it here. This parameter will be used only if the option compareSubOptimal is active/True. Default: AS (AS is usually written by BWA-MEM)", action="store", default="AS")
parser.add_argument('--XSlikeTag', dest="XSlikeTag", help="If the BAM file contains a suboptimal alignment score tag different from 'XS', change it here. This parameter will be used only if the option compareSubOptimal is active/True. Default: XS (XS is usually written by BWA-MEM)", action="store", default="XS")
parser.add_argument('--SAalnAction', dest="SAalnAction", help="If the BAM file contains SA alignments (Secondary Alignments, that could be also chimera), you can choose which action to do among {ignore, remove, removeByComparison}. This parameter will be used only if the option compareSubOptimal is active/True. Default: ignore.", action="store", default="ignore")

parser.add_argument('-o', '--outfilename', dest="outfilename", help="Output file name. The output will be a list of headers to be removed. No default option.", action="store", required=True)
parser.add_argument('-p', '--processes', dest="processes", help="Number of processes to use. Default 2.", action="store", default=2)


args = parser.parse_args()

#########################################################################################
####### GLOBAL VARS
#########################################################################################

# init db values -> from UCSC database ChromInfo
hg19="1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566;M:1-16571"
mm9="1:1-197195432;2:1-181748087;3:1-159599783;4:1-155630120;5:1-152537259;6:1-149517037;7:1-152524553;8:1-131738871;9:1-124076172;10:1-129993255;11:1-121843856;12:1-121257530;13:1-120284312;14:1-125194864;15:1-103494974;16:1-98319150;17:1-95272651;18:1-90772031;19:1-61342430;X:1-166650296;Y:1-15902555;M:1-16299"
# BGI macaca fascicularis genome
bgice="1:1-232298473;2:1-192287759;3:1-199730254;4:1-169826420;5:1-185160406;6:1-180941566;7:1-171133315;8:1-150570238;9:1-133536003;10:1-96524937;11:1-136666246;12:1-108081426;13:1-137686314;14:1-134145161;15:1-108996402;16:1-81259162;17:1-95628244;18:1-75918782;19:1-65364038;20:1-89811811;X:1-155426724;Ur:1-453853"
# NCBI macaca fascicularis genome
mfa5="1:1-227556264;2:1-192460366;3:1-192294377;4:1-170955103;5:1-189454096;6:1-181584905;7:1-171882078;8:1-146850525;9:1-133195287;10:1-96509753;11:1-137757926;12:1-132586672;13:1-111193037;14:1-130733371;15:1-112612857;16:1-80997621;17:1-96864807;18:1-75711847;19:1-59248254;20:1-78541002;X:1-152835861;M:1-16575"
# The following genomes come from TIGET sequences, coded as follows:
# - lv.backbone.fa: --> lv
# 	chr1 = LTR + vecore genome; 
# 	chr2 = dnef
lv="1:1-1748;2:1-82"
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
# - transposon.cd123.fa --> transposon
# 	chr1 = full plasmid sequence
transposon="1:1-4887"
# - retro.ada.fa --> giada
# 	chr1 = full plasmid sequence
giada="1:1-6185"

# Custom genome for RNA datasets
# - lv.743.fa
lv743="743.pCCLsin.PPT.hPGK.GFP.Wpre_mut_AMP:1-3968"
# - lv.743.fa
lv3029="3029.pCCL.223.142.3pLTRmir.sin.PPT.BpA.GFP.PGKnoSA.wpre:1-4513"


readorder = ['chr', 'start', 'end', 'strand', 'RG', 'quality', 'NM', 'flag', 'cigar', 'MD', 'insert_size', 'AS', 'XS', 'SA', 'nasequence'] # this is the order of the content to written in the DB for each read (and pair)
locus_prefix = 'isread_' # prefix of the integration site read for DB columns
mate_prefix = 'mate_' # prefix of the mate read for DB columns
readdbschema = {	'chr': "varchar(50) DEFAULT NULL",
					'start': "int(30) NOT NULL",
					'end': "int(30) NOT NULL",
					'strand': "varchar(50) DEFAULT NULL",
					'quality': "int(10) NOT NULL",
					'NM': "int(10) DEFAULT NULL",
					'flag': "varchar(100) DEFAULT NULL",
					'cigar': "varchar(255) DEFAULT NULL",
					'MD': "varchar(255) DEFAULT NULL",
					'insert_size': "int(30) DEFAULT NULL",
					'AS': "int(30) DEFAULT NULL",
					'XS': "int(30) DEFAULT NULL",
					'SA': "varchar(255) DEFAULT NULL",
					'quality': "int(30) DEFAULT NULL",
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
	if args.indexBAM: # index BAM file
		indexBAM(args.bamfile)
	if args.sortBAM: # sort BAM and assign the new file to ARGS and index it
		sorted_bamfile = sortBAM(args.bamfile)
		args.bamfile = sorted_bamfile
		indexBAM(args.bamfile)
	if args.locusRead not in ['read1', 'read2', 'both', 'alnstart']:
		print "\n[AP]\tError while reading locusRead option: no valid alternative selected. See help for instructions. You set %s\n" %(args.locusRead)
		sys.exit()
	if (args.locusRead in ['alnstart'] and not args.singleEnd) or (args.locusRead not in ['alnstart'] and args.singleEnd):
		print "\n[AP]\tError in input request, it seems that you are asking to import single end data but that the locus reads is set as paired-ends read. With --singleEnd you must use --locusRead alnstart. See help for instructions.\n"
		sys.exit()
# 	if args.poolID in ['', None] or args.associationID in ['', None]:
# 		print "\n[AP]\tError in IDs options: Association and Pool ID must be not NULL/None nor empty.\n"
# 		sys.exit()
	if args.tag is None or args.tag == '':
		print "\n[AP]\tError in barcode/tag IDs options: it must be not NULL/None nor empty.\n"
		sys.exit()

def parseAssociationfile(infile):
	"""
	IN: association file
	OUT: dictionary of TAGs (k=tag, v=dict: k={tissue, sample, timepoint}, v=values)
	FILE FORMAT: <barcode_id> <barcode_sequence> <tissue> <sample> <timepoint> <lam_id> <lam_fullname> <marker> <enzyme> <vector> <associationID> <poolID>
	"""
	assodict = {}
	with open(infile, 'rb') as inf:
		reader = csv.reader(inf, delimiter = "\t")
		for row in reader:
			if len(row) > 0:
				barcode_id = row[0] 
				barcode_sequence = row[1]
				tissue = row[2]
				sample = row[3]
				timepoint = row[4]
				assodict[barcode_sequence] = {
					"tissue": tissue,
					"sample": sample,
					"timepoint": timepoint,
					"barcode_id": barcode_id,
					"lam_id": row[5],
					"lam_fullname": row[6],
					"marker": row[7],
					"enzyme": row[8],
					"vector": row[9],
					"associationID": row[-2],
					"poolID": row[-1],
					}
			else:
				print "[AP]\tWarning: this file contais a blank row!!!! Check it please."
	return assodict

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

def complementString(instr, sequencetype = 'DNA'): 
	"""Return the complementary sequence string.""" 
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'} 
	if sequencetype == 'RNA':
		basecomplement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'N': 'N', '-': '-'} 
	letters = list(instr) 
	letters = [basecomplement[base] for base in letters] 
	return ''.join(letters) 

def reversecomplement(instr, sequencetype = 'DNA'):
	"""Return the reverse complement of the dna string.""" 
	reverse_sequence = revertString(instr)
	reversecomplement_sequence = complementString(reverse_sequence)
	return reversecomplement_sequence

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


def worker(todo_regions, bam, singleEnd, locusRead, out_q, ): 
	""" 
	The worker function, invoked in a process. 
		region: array of tuples with coordinate to slice chromosome
		bam: the input bam to splice
		out_q: queue of the process
	"""
	outdata = []
	for region in todo_regions:
		chr, start, end = region
		product = {} # k = header, v = 'r1': {'start', 'end'}, 'r2': {'start', 'end'}
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t...processing region (chr, start, end)", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		
		print "[AP]\t\t\tAcquiring dictionary of reads in the region (chr, start, enc)\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		for bam_alignment in mybam: # for each alignment in the BAM file (both R1 and R2!!!)
			# get orientation and check for correct position (get end from orientation and start)
			orientation = '+' # default, if fwd
			aln_start = bam_alignment.pos # default, if fwd
			aln_end = bam_alignment.pos + bam_alignment.alen # default, if fwd
			nasequence = bam_alignment.seq # set up read sequence, default fwd strand
			if bam_alignment.is_reverse:
				orientation = '-'
				## aln_end = bam_alignment.pos
				## aln_start = bam_alignment.pos + bam_alignment.alen
				## nasequence = reversecomplement(bam_alignment.seq) 

			# acquire tags
			AS = "NULL" # best hit alignment score
			XS = "NULL" # suboptimal alignment score
			SA = "NULL" # secondary alignments / chimera
			NM = "NULL" # number of mismatches
			MD = "NULL" # md score from cigar
			RG = "NULL" # read group
			for tag, val in bam_alignment.tags: # acquire data/values from tags
				if tag == 'AS': # best hit alignment AS score
					AS = val
				if tag == 'XS': # XS tag, suboptimal alignment score
					XS = val
				if tag == 'SA': # XS tag, suboptimal alignment score
					SA = val
				if tag == 'MD': # XS tag, suboptimal alignment score
					MD = val
				if tag == 'NM': # XS tag, suboptimal alignment score
					NM = val
				if tag == 'RG': # XS tag, suboptimal alignment score
					RG = val
			# paired ends experiment case
			if not singleEnd:
				if not bam_alignment.is_read2: # all but not read 2: this means read 1 and potentially all singletones from r1
					if product.has_key(bam_alignment.qname):
						product[bam_alignment.qname]['r1'] = {		'start': aln_start, 
																	'end': aln_end, 
																	'strand': orientation , 
																	'AS': AS,
																	'XS': XS,
																	'SA': SA,
																	'NM': NM,
																	'MD': MD,
																	'RG': RG,
																	'chr': bam.getrname(bam_alignment.tid), 
																	'cigar': bam_alignment.cigarstring,
																	'flag': bam_alignment.flag,
																	'insert_size': bam_alignment.tlen ,
																	'quality': bam_alignment.mapq,
																	'nasequence': nasequence,
															}
					else:
						product[bam_alignment.qname] = {	'r1': {	'start': aln_start, 
																	'end': aln_end, 
																	'strand': orientation , 
																	'AS': AS,
																	'XS': XS,
																	'SA': SA,
																	'NM': NM,
																	'MD': MD,
																	'RG': RG,
																	'chr': bam.getrname(bam_alignment.tid), 
																	'cigar': bam_alignment.cigarstring,
																	'flag': bam_alignment.flag,
																	'insert_size': bam_alignment.tlen ,
																	'quality': bam_alignment.mapq,
																	'nasequence': nasequence,
															}
													}
				else: # the r2 reads
					if product.has_key(bam_alignment.qname):
						product[bam_alignment.qname]['r2'] = {	'start': aln_start, 
																'end': aln_end, 
																'strand': orientation , 
																'AS': AS,
																'XS': XS,
																'SA': SA,
																'NM': NM,
																'MD': MD,
																'RG': RG,
																'chr': bam.getrname(bam_alignment.tid), 
																'cigar': bam_alignment.cigarstring,
																'flag': bam_alignment.flag,
																'insert_size': bam_alignment.tlen ,
																'quality': bam_alignment.mapq,
																'nasequence': nasequence,
																} 
						#print bam_alignment.qname
					else:
						product[bam_alignment.qname] = {	'r2': {	'start': aln_start, 
																	'end': aln_end, 
																	'strand': orientation , 
																	'AS': AS,
																	'XS': XS,
																	'SA': SA,
																	'NM': NM,
																	'MD': MD,
																	'RG': RG,
																	'chr': bam.getrname(bam_alignment.tid), 
																	'cigar': bam_alignment.cigarstring,
																	'flag': bam_alignment.flag,
																	'insert_size': bam_alignment.tlen ,
																	'quality': bam_alignment.mapq,
																	'nasequence': nasequence,
																}
													}
			# single end experiment
			else:
				product[bam_alignment.qname] = { 'r': {	'start': aln_start, 
														'end': aln_end, 
														'strand': orientation , 
														'AS': AS,
														'XS': XS,
														'SA': SA,
														'NM': NM,
														'MD': MD,
														'RG': RG,
														'chr': bam.getrname(bam_alignment.tid), 
														'cigar': bam_alignment.cigarstring,
														'flag': bam_alignment.flag,
														'insert_size': bam_alignment.tlen ,
														'quality': bam_alignment.mapq,
														'nasequence': nasequence,
														}
												}
		# now loop over the dictionary keys to create the output array based on the driver read
		# forma dell'output: 'prod_header', 'prod_chr', 'prod_locus', 'prod_end', 'prod_strand', 'ref_associationid', 'ref_matrixid', 'ref_poolid', 'isread_start', 'isread_end', 'isread_strand', 'isread_chr', 'isread_cigar', 'isread_flag', 'isread_insert_size', 'isread_AS', 'isread_XS', 'isread_SA', 'isread_NM', 'isread_MD', 'mate_start', 'mate_end', 'mate_strand', 'mate_chr', 'mate_cigar', 'mate_flag', 'mate_insert_size', 'mate_AS', 'mate_XS', 'mate_SA', 'mate_NM', 'mate_MD',
		
		if locusRead == 'read1': # locus is in r1
			locus = 'r1'
			mate = 'r2'
			for k, v in product.iteritems():
				# set prod end 
				prod_end = None
				if v.has_key(mate): # if it has the pair -> this is the end
					prod_end = v[mate]['start']
				if v.has_key(locus):
					if prod_end is None:
						prod_end = v[locus]['end'] # if no pair -> this prod makes the end
					is_array = [k, v[locus]['chr'], v[locus]['start'], prod_end, v[locus]['strand'], 'associationID', 'NULL', dict_association[barcodeID]['poolID'], ]
					for x in readorder: # first the locus read
						is_array.append( str(v[locus][x]) )
					if v.has_key(mate):# then the mate read, if it exists
						for x in readorder: 
							is_array.append( str(v[mate][x]) )
					else:
						for x in readorder: 
							is_array.append( 'NULL' )
						####### -> ???? is_array.append( str('NULL') for x in range(0,len(readorder)))
					outdata.append(is_array) # append results to output
		elif locusRead == 'read2': # locus in r2
			locus = 'r2'
			mate = 'r1'
			for k, v in product.iteritems():
				# set prod end 
				prod_end = None
				if v.has_key(mate): # if it has the pair -> this is the end
					prod_end = v[mate]['start']
				if v.has_key(locus):
					if prod_end is None:
						prod_end = v[locus]['end'] # if no pair -> this prod makes the end
					is_array = [k, v[locus]['chr'], v[locus]['start'], prod_end, v[locus]['strand'], dict_association[barcodeID]['associationID'], 'NULL', dict_association[barcodeID]['poolID'], ]
					for x in readorder: # first the locus read
						is_array.append( str(v[locus][x]) )
					if v.has_key(mate):# then the mate read, if it exists
						for x in readorder: 
							is_array.append( str(v[mate][x]) )
					else:
						for x in readorder: 
							is_array.append( 'NULL' )
						####### -> ???? is_array.append( str('NULL') for x in range(0,len(readorder)))
					outdata.append(is_array) # append results to output
		elif locusRead == 'alnstart': # locus by alignment start
			locus = 'r'
			for k, v in product.iteritems():
				is_array = [k, v[locus]['chr'], v[locus]['start'], v[locus]['end'], v[locus]['strand'], dict_association[barcodeID]['associationID'], 'NULL', dict_association[barcodeID]['poolID'], ]
				for x in readorder: # first the locus read
					is_array.append( str(v[locus][x]) )
				for x in readorder: 
					is_array.append( 'NULL' )
				outdata.append(is_array) # append results to output
		else:
			print "[AP]\tERROR:\tNo valid locus read option.\n"
			sys.exit()
			
		print "[AP]\t\t\tDone, region completed (chr, start, end)\t\t\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" 
	
	#print len(outdata), outdata[0]
	# queue results at the end of this process
	out_q.put(outdata)

def importData(dbfile, dbtable, dbschema, localfile, resultdata):
	"""
	Import data into DB MySQL using system call.
	Now use the package MySQLdb because there are some problems with the system call: not all rows are written even if they exist. Running the script I can import 98 rows, from the shell the same call imports all 100 onces.
	"""
	conn = sqlite3.connect(dbfile)

	#cursor = conn.cursor ()
	cursor = conn.cursor()
	# create table
	## read order :: 'chr', 'start', 'end', 'strand', 'quality', 'NM', 'flag', 'cigar', 'MD', 'insert_size', 'AS', 'XS', 'SA'
	createtable = """
		CREATE TABLE IF NOT EXISTS %s (  
			read_header varchar(255) NOT NULL, 
			read_r1 int(1) DEFAULT 1,
			read_chr varchar(255) DEFAULT NULL, 
			read_start int(30) NOT NULL, 
			read_end int(30) NOT NULL, 
			read_strand varchar(255) DEFAULT NULL, 
			read_RG varchar(250) DEFAULT NULL, 
			read_quality int(10) NOT NULL, 
			read_NM int(30) DEFAULT NULL, 
			read_flag varchar(255) DEFAULT NULL, 
			read_cigar varchar(255) DEFAULT NULL, 
			read_MD varchar(255) DEFAULT NULL, 
			read_insert_size int(30) DEFAULT NULL, 
			read_AS int(30) DEFAULT NULL, 
			read_XS int(30) DEFAULT NULL, 
			read_SA varchar(255) DEFAULT NULL, 
			-- read_nasequence varchar(1000) DEFAULT NULL, 
			PRIMARY KEY  ( read_header , read_r1 , read_chr , read_start , read_end )
		) ENGINE=MyISAM DEFAULT CHARSET=latin1; 
		""" %(dbtable)
	#print createtable
	##os.system("""mysql -u %s --password=%s -h %s --port %s %s -e "%s" """ %(dbuser, dbpassword, dbhost, dbport, dbschema, createtable) )
	cursor.execute(createtable)

	# import data into table
	##os.system("mysqlimport -u %s --password=%s -h %s --local --port %s %s %s" %(dbuser, dbpassword, dbhost, dbport, dbschema, localfile))

	### now format data to import all and insert them all
	step = 5000
	bedstart = 0 # starting point of bed data array for each step
	bedend = 0 # starting point of bed data array for each step
	# check max len of beddata to set bedend, no more than beddata len!
	if len(resultdata) > step:
		bedend = step
	else:
		bedend = len(resultdata)
	# now process data into steps of a maximum number of blocks (step)
	for elems in range(0, len(resultdata), step):
		query_import = """
			INSERT INTO `%s`.`%s` (
				`read_header` ,
				`read_chr` ,
				`read_start` ,
				`read_end` ,
				`read_strand` ,
				`read_RG` ,
				`read_quality` ,
				`read_NM` ,
				`read_flag` ,
				`read_cigar` ,
				`read_MD` ,
				`read_insert_size` ,
				`read_AS` ,
				`read_XS` ,
				`read_SA`
				)
			VALUES
			""" %(dbschema, dbtable)

		# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
		for row in resultdata[bedstart:bedend]:
			query_import += " ('" + "', '".join( str(x).replace("'","") for x in row ) + "'), "
		#print query_import.rstrip(', ')
		importdata = query_import.rstrip(', ').replace("'NULL'", "NULL")
		#print importdata
		cursor.execute( importdata )
		# increment bed indexes
		bedstart += step
		bedend += step

	cursor.close()
	conn.close()
	# return value
	return 1


#########################################################################################
### MAIN
#########################################################################################

def main():
	"""
	Main part of the program.
	"""
	start_time = time.time() # get read of starting time
	# first check args and file paths
	print "[AP]\tChecking inputs."
	checkArgs(args)
	# acquire thresholds
	minStartingMatches = int(args.minStartingMatches)
	minStartingBasesNoIndels = int(args.minStartingBasesNoIndels)
	endClipThreshold = int(args.endClipThreshold)
	deltaAlignment = int(args.deltaAlignment)
	# if you do not find the index, do it first
	if not (os.path.isfile(args.bamfile1 + ".bai") or os.path.isfile(os.path.splitext(args.bamfile1)[0] + ".bai")):
		print "[AP]\tIndexing BAM 1 (this is required if you did not do it)."
		indexing = pysam.index(args.bamfile1) # indexing is BLANK
	if not (os.path.isfile(args.bamfile2 + ".bai") or os.path.isfile(os.path.splitext(args.bamfile2)[0] + ".bai")):
		print "[AP]\tIndexing BAM 2 (this is required if you did not do it)."
		indexing = pysam.index(args.bamfile2) # indexing is BLANK
	# acquire BAM files: INPUT data
	print "[AP]\tAcquiring data BAM from both files."
	bam1 = pysam.Samfile(args.bamfile1) # the bam file to process
	bam2 = pysam.Samfile(args.bamfile2) # the bam file to process
	print "[AP]\tAcquiring region list to use while slicing BAM."
	regions1 = None
	if args.genome1 is not "regions":
		if args.genome1 == "hg19":
			regions1 = getRegions(hg19)
		elif args.genome1 == "mm9":
			regions1 = getRegions(mm9)
		elif args.genome1 == "mfa5":
			regions1 = getRegions(mfa5)
		elif args.genome1 == "ce":
			regions1 = getRegions(bgice)
		elif args.genome1 == "lv":
			regions1 = getRegions(lv)
		elif args.genome1 == "lvwas":
			regions1 = getRegions(lvwas)
		elif args.genome1 == "lvarsa":
			regions1 = getRegions(lvarsa)
		elif args.genome1 == "lvamp":
			regions1 = getRegions(lvamp)
		elif args.genome1 == "lvkana":
			regions1 = getRegions(lvkana)
		elif args.genome1 == "transposon":
			regions1 = getRegions(transposon)
		elif args.genome1 == "giada":
			regions1 = getRegions(giada)
		elif args.genome1 == "lv743":
			regions1 = getRegions(lv743)
		elif args.genome1 == "lv3029":
			regions1 = getRegions(lv3029)
		else:
			print "\t\tLoading regions1 error! Please read the help (-h)."
			sys.exit()
	else:
		print "[AP]\tSorry, the option regions1 is not used here, please do not use it.\n\n"
		sys.exit()

	regions2 = None
	if args.genome2 is not "regions":
		if args.genome2 == "hg19":
			regions2 = getRegions(hg19)
		elif args.genome2 == "mm9":
			regions2 = getRegions(mm9)
		elif args.genome2 == "mfa5":
			regions2 = getRegions(mfa5)
		elif args.genome2 == "ce":
			regions2 = getRegions(bgice)
		elif args.genome2 == "lv":
			regions2 = getRegions(lv)
		elif args.genome2 == "lvwas":
			regions2 = getRegions(lvwas)
		elif args.genome2 == "lvarsa":
			regions2 = getRegions(lvarsa)
		elif args.genome2 == "lvamp":
			regions2 = getRegions(lvamp)
		elif args.genome2 == "lvkana":
			regions2 = getRegions(lvkana)
		elif args.genome2 == "transposon":
			regions2 = getRegions(transposon)
		elif args.genome2 == "giada":
			regions2 = getRegions(giada)
		elif args.genome2 == "lv743":
			regions2 = getRegions(lv743)
		elif args.genome2 == "lv3029":
			regions2 = getRegions(lv3029)
		else:
			print "\t\tLoading regions2 error! Please read the help (-h)."
			sys.exit()
	else:
		print "[AP]\tSorry, the option regions2 is not used here, please do not use it.\n\n"
		sys.exit()

	print "[AP]\tNow looping over BAM alignments Splitting run in processes..."
	out_q = Queue() # init process queue
	nprocs = int(args.processes)
	chunksize = int(math.ceil(len(regions) / float(nprocs))) # find the chunk wrt processes for running regions
	procs = [] # array of processes
	# run processes in a range of chunks
	for i in range(nprocs):
		p = multiprocessing.Process(target=worker, args=(regions[chunksize * i:chunksize * (i + 1)], bam, minStartingMatches, minStartingBasesNoIndels, endClipThreshold, args.pruneByPair, out_q, args.single_end, args.compareSubOptimal, args.suboptimalThreshold, args.ASlikeTag, args.XSlikeTag, args.SAalnAction )) # prototype: todo_regions, bam, bed_dictionary, minStartingMatches, out_q
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
		filetarget.write(row + "\n")
	filetarget.close() # finalize file
	
	elapsed_time = time.time() - start_time # get read of finishing time
	print "\n[AP]\tTask Finished, closing.\n\tElapsed time: %.2f [seconds]\n" %(elapsed_time)
	

# sentinel
if __name__ == "__main__":
    main()
