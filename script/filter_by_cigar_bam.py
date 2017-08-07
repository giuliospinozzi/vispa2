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
        ***   FILTER BY CIGAR   ***
+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		May 2013 - October 2013
 Contact:	andrea.calabria@hsr.it
 Revision:	0.4 (MP, platform independent, indels, aln score, SA/chimeras)
 			0.5 18 June 2014 (bug fixed for single end evaluation and CIGAR rules)
 			0.6 19 June 2014 (new rule for NM evaluation)
 			0.7 19 August 2014 (change of mismatch evaluation: added fromBP toBP)
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
parser.add_argument('--bam', dest="bamfile", help="BAM file to process. No default option.", action="store", required=True)
parser.add_argument('--minStartingMatches', dest="minStartingMatches", help="Minimum starting matches to consider the read valid [Integer]. If a read has less than minStartingMatches, it will be discarded. Default = 3 (-> >=3 are valid, else discarded).", action="store", required=False, default=3)
parser.add_argument('--minStartingMatches_fromNbp', dest="minStartingMatches_fromNbp", help="Initial tolerance threshold before considering the minimum starting matches to consider the read valid [Integer]. A read must have the first N bases (minStartingMatches_fromNbp) of tolerance before the N bases of minStartingMatches. Default = 0.", action="store", required=False, default=0)
parser.add_argument('--minStartingBasesNoIndels', dest="minStartingBasesNoIndels", help="Minimum number of starting bases in which no insertions nor deletions must occur (orientation dependent). This number contains minStartingMatches and minStartingMatches_fromNbp, indeed the program will parse for indels in the interval (minStartingMatches_fromNbp + minStartingMatches, minStartingBasesNoIndels]. Write an Integer value. If a read has an indel within minStartingBasesNoIndels bass, it will be discarded. Default = 7 (-> from 1 to 7 no indels. <=7 are valid, else discarded).", action="store", required=False, default=7)
parser.add_argument('--endClipThreshold', dest="endClipThreshold", help="Maximum number of bases for which R1 can end without M matches. This is typically useful for soft clip bases: if the end of R1 has a number of soft clipped bases greater than threshold, this read is removed. Default = 10 (-> >10 are filtered out).", action="store", required=False, default=10)
#parser.add_argument('--headerformat', dest="headerformat", help="Header format of the filter file. Cases: type 1:: '@M00174:25:000000000-A0T21:1:1:15041:1491 1:N:0:0'; type 2:: '@M00571:5:000000000-A267F:1:1101:14475:1610/1'. Select type number {1, .., N}. No default value assigned because you must explicitely be aware of what you are doing.", action="store", required=True)
#parser.add_argument('--singleread', dest="singleread", action="store_true", help="In case of single read experiment, not paired-end. Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--sortBAM', dest="sortBAM", action="store_true", help="Sort BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--indexBAM', dest="indexBAM", action="store_true", help="Index BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--singleEnd', dest="single_end", action="store_true", help="Are your reads single-end and not paired-end? If your reads are single-end, activate this option. Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--pruneByPair', dest="pruneByPair", action="store", help="If you have a specific pair that drives the filter (e.g. read1 or read2), specify it here {both, read1, read2}. This is useful if your experimental design attemps to analyze at the end only one of the two pairs, for example all integration sites come from read1 only. If your reads are single-end this option is not activated. Default: 'both'.", default='both')

parser.add_argument('--compareSubOptimal', dest="compareSubOptimal", action="store_true", help="Do you want to compare the alignment score versus the suboptimal alignment score? If yes, here you are going to compare the AS tag (that you can change in the ASlikeTag option) and the XS tag (set it in the XSlikeTag option), as output from BWA MEM. Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--suboptimalThreshold', dest="suboptimalThreshold", help="Specify the suboptimal threshold over which to keep the read [Integer]. If a read has less the best hit alignment score that is different (higher) from the suboptimal alignment score for a percentage suboptimalThreshold then this read is kept, otherwise is discarded. Example: AS=100, XS=80 delta=(AS - XS/AS*100)=20, delta >= suboptimalThreshold? If yes, this read passes filter (is kept). Default: 10 (-> delta>=10 :: valid, else discarded).", action="store", required=False, default=10)
parser.add_argument('--ASlikeTag', dest="ASlikeTag", help="If the BAM file contains an alignment score tag different from 'AS', change it here. This parameter will be used only if the option compareSubOptimal is active/True. Default: AS (AS is usually written by BWA-MEM)", action="store", default="AS")
parser.add_argument('--XSlikeTag', dest="XSlikeTag", help="If the BAM file contains a suboptimal alignment score tag different from 'XS', change it here. This parameter will be used only if the option compareSubOptimal is active/True. Default: XS (XS is usually written by BWA-MEM)", action="store", default="XS")
parser.add_argument('--SAalnAction', dest="SAalnAction", help="If the BAM file contains SA alignments (Secondary Alignments, that could be also chimera), you can choose which action to do among {ignore, remove, removeByComparison}. This parameter will be used only if the option compareSubOptimal is active/True. Default: ignore.", action="store", default="ignore")

parser.add_argument('-o', '--outfilename', dest="outfilename", help="Output file name. The output will be a list of headers to be removed. No default option.", action="store", required=True)
parser.add_argument('-p', '--processes', dest="processes", help="Number of processes to use. Default 2.", action="store", default=2)

parser.add_argument('-r', '--regions', dest="regions", help="Regions to process as CSV string in the format: CHR:START-END. Default value is for HG19: '1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566'.", action="store", default="1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566")
parser.add_argument('-g', '--genome', dest="genome", help="Select the whole genome with relative regions (this is alterative to regions). If you want to load regions, write this as 'regions'. Default: regions. Values: {regions, hg19, mm9, mm10, mfa5, ce, lv, lvarsa, lvwas, lvkana, lvamp, transposon, giada, lv743, lv3029, hiv}.", action="store", default="regions")

parser.add_argument('--compareNM', dest="compareNM", action="store_true", help="Do you want to compare the alignment NM versus the suboptimal alignment? This is valid only if compareSubOptimal is True. The option is boolean and will activate a contrain relaxation of the suboptimalThreshold parameter by a 20 percent and will evaluate the number of mismatches (NM): if aNM/100bp and a=1, then keep it, else discard the read. Default: 'False'; alternative: 'True' or just activate the option.", default=False)

args = parser.parse_args()

#########################################################################################
####### GLOBAL VARS
#########################################################################################

# init db values -> from UCSC database ChromInfo
### source backup genomes
# hg19="1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566;M:1-16571"
# mm9="1:1-197195432;2:1-181748087;3:1-159599783;4:1-155630120;5:1-152537259;6:1-149517037;7:1-152524553;8:1-131738871;9:1-124076172;10:1-129993255;11:1-121843856;12:1-121257530;13:1-120284312;14:1-125194864;15:1-103494974;16:1-98319150;17:1-95272651;18:1-90772031;19:1-61342430;X:1-166650296;Y:1-15902555;M:1-16299"
# no chrM
hg19="1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566"
mm9="1:1-197195432;2:1-181748087;3:1-159599783;4:1-155630120;5:1-152537259;6:1-149517037;7:1-152524553;8:1-131738871;9:1-124076172;10:1-129993255;11:1-121843856;12:1-121257530;13:1-120284312;14:1-125194864;15:1-103494974;16:1-98319150;17:1-95272651;18:1-90772031;19:1-61342430;X:1-166650296;Y:1-15902555"
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
	if args.indexBAM: # index BAM file
		indexBAM(args.bamfile)
	if args.sortBAM: # sort BAM and assign the new file to ARGS and index it
		sorted_bamfile = sortBAM(args.bamfile)
		args.bamfile = sorted_bamfile
		indexBAM(args.bamfile)
	if args.compareSubOptimal and args.SAalnAction is not None:
		if args.SAalnAction not in ['ignore', 'remove', 'removeByComparison']:
			print "\n[AP]\tError in the option SAalnAction: no valid argument passed, you specified [%s]\n\n\tExiting...\n" %(args.SAalnAction)
			sys.exit()
	if args.pruneByPair not in ['both', 'read1', 'read2']:
		print "\n[AP]\tError in the option --pruneByPair: no valid argument passed, you specified [%s]\n\n\tExiting...\n" %(args.SAalnAction)
		sys.exit()
	
	
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

def evaluateCigar_OBSOLETE(bam_aln):
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

def evaluateCigar_v06(bam_aln, min_starting_matches = 3, minStartingBasesNoIndels = 3, endClipThreshold = 10):
	"""
	Input: pysam alignment object
	Output: true or false
	Logics: deletions or indels not in the first N bases (N = minStartingBasesNoIndels) with respect to the orientation
		- find mapping orientation
		- based on the orientation, if you have at least M in the first positions and NOT D or I in the second position after at least min_bases_aftermatches_forgaps bases, and no H within the first 20 bp, consider this read as valid (-> it means: if minStartingBasesNoIndels=15, 20M2D is valid, 10M2D is not valid)
	Based on strand, evaluate possible problems in CIGAR string. 
	Cigar def: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/glossary.html#term-cigar
	Operations and len are referred to main manual: http://samtools.sourceforge.net/SAM1.pdf
	We could take into account number ops: 0, 7, 8 -> M, =, X. So far is only M [TODO]
	If the R1 ends with tag not M of val > endClipThreshold => discard it.

	Cases (if reverse strand, the CIGAR must be reverted, in the stdout reverse strand are labeld with /rev:
	(0) nM -> ok | n!M -> discarded
	(1) nM -> ok iif n >= min_starting_matches | nM with n < min_starting_matches OR n!M(iI|iD|iH)+ [that is: aln does not start with M] -> discarded
	(2) nM(aI|aD)+: indels in the interval minStartingBasesNoIndels (n+a) -> remove
	(3) nM(nI|nD)*(aX|aH)+: ends with hard clipping or padding -> remove
	(4) nM(nI|nD)*iS: ends with soft clipping where i > endClipThreshold -> remove, if i <= endClipThreshold -> keep it
	"""
	valid = True
	# the following values of the cigar flag come from the PYSAM documentation 
	#  http://sqt.googlecode.com/git-history/hashing/sqt/cigar.py
	#  http://www.cgat.org/~andreas/documentation/pysam/api.html

	if len(bam_aln.cigar) > 1: # if you have an alignment with more than 1 cigar tag
		if bam_aln.is_reverse: # if the read is in reverse, start counting from the end of the CIGAR list	
			if (bam_aln.cigar[-1][0] == cigar_flag['M'] and int(bam_aln.cigar[-1][1]) >= min_starting_matches): # if first is M: nM where min(n) = min_starting_matches
				if int(bam_aln.cigar[-1][1]) <= minStartingBasesNoIndels and bam_aln.cigar[-2][0] in [cigar_flag['I'], cigar_flag['D']] and int(bam_aln.cigar[-2][1]) > 0:
					print "[AP]\t\t%s -> removed read by CIGAR: case (2) indels within (<=) %d bp\t[chr%d:%d-%s/rev, len %s, CIGAR %s]" %(bam_aln.qname, minStartingBasesNoIndels, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
					valid = False
				elif bam_aln.cigar[0][0] in [cigar_flag['H'], cigar_flag['P']] and int(bam_aln.cigar[0][1]) > 0:
					print "[AP]\t\t%s -> removed read by CIGAR: case (3) indels within (<=) %d bp\t[chr%d:%d-%s/rev, len %s, CIGAR %s]" %(bam_aln.qname, minStartingBasesNoIndels, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
					valid = False
				elif bam_aln.cigar[0][0] in [cigar_flag['S']] and int(bam_aln.cigar[0][1]) > endClipThreshold:
					print "[AP]\t\t%s -> removed read by CIGAR: case (4) ends with S soft clipping > %d bp\t[chr%d:%d-%s/rev, len %s, CIGAR %s]" %(bam_aln.qname, endClipThreshold, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
					valid = False
				# if all conditions above are passed, then this read is valid (default)
			else: # nM with n<=min_starting_matches
				print "[AP]\t\t%s -> removed read by CIGAR: case (1) nM has n <= %d\t[chr%d:%d-%s/rev, len %s, CIGAR %s]" %(bam_aln.qname, min_starting_matches, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
				valid = False
		else:
			if (bam_aln.cigar[0][0] == cigar_flag['M'] and int(bam_aln.cigar[0][1]) >= min_starting_matches): # if first is M: nM where min(n) = min_starting_matches
				if int(bam_aln.cigar[0][1]) <= minStartingBasesNoIndels and bam_aln.cigar[1][0] in [cigar_flag['I'], cigar_flag['D']] and int(bam_aln.cigar[1][1]) > 0:
					print "[AP]\t\t%s -> removed read by CIGAR: case (2) indels within (<=) %d bp\t[chr%d:%d-%s, len %s, CIGAR %s]" %(bam_aln.qname, minStartingBasesNoIndels, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
					valid = False
				elif bam_aln.cigar[-1][0] in [cigar_flag['H'], cigar_flag['P']] and int(bam_aln.cigar[-1][1]) > 0:
					print "[AP]\t\t%s -> removed read by CIGAR: case (3) indels within (<=) %d bp\t[chr%d:%d-%s, len %s, CIGAR %s]" %(bam_aln.qname, minStartingBasesNoIndels, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
					valid = False
				elif bam_aln.cigar[-1][0] in [cigar_flag['S']] and int(bam_aln.cigar[-1][1]) > endClipThreshold:
					print "[AP]\t\t%s -> removed read by CIGAR: case (4) ends with S soft clipping > %d bp\t[chr%d:%d-%s, len %s, CIGAR %s]" %(bam_aln.qname, endClipThreshold, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
					valid = False
				# if all conditions above are passed, then this read is valid (default)
			else: # nM with n<=min_starting_matches
				print "[AP]\t\t%s -> removed read by CIGAR: case (1) nM has n <= %d\t[chr%d:%d-%s, len %s, CIGAR %s]" %(bam_aln.qname, min_starting_matches, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
				valid = False
	else: # since you have only 1 tag for the cigar string, orientation (in this case) does not matter
		if (bam_aln.cigar[0][0] == cigar_flag['M'] and int(bam_aln.cigar[-1][1]) > 0): # only if you have M -> ok
			valid = True
		else:
			print "[AP]\t\t%s -> removed read by CIGAR: case (0) n!M \t[chr%d:%d-%s, len %s, CIGAR %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring)
			valid = False
			#print bam_aln.qname, bam_aln.cigar, bam_aln.tags
	
	if bam_aln.cigar is None:
		print "[AP]\tWarning: NO CIGAR (=None) for read %s" %(bam_aln.qname)
		valid = False
		
	return valid

def getCigarSequenceFromCigarList(bam_aln_cigar_list, ):
	"""
	from the cigar list like [(4,2), (0,5), (4,4)] tp the corresponding sequence of chars SSMMMMMSSSS
	"""
	cigar_sequence = ''.join( str(cigar_number[tag])*n_bp for tag, n_bp in bam_aln_cigar_list)
	return cigar_sequence


def evaluateCigar(bam_aln, min_starting_matches = 3, minStartingBasesNoIndels = 3, endClipThreshold = 10, min_starting_matches_afterNbp = 0):
	"""
	Input: pysam alignment object, min_starting_matches_afterNbp:: number of min bp of tolerance (0=no tolerance, 2=the first 2 bp of tolerance)
	Output: true or false
	Cigar def: http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/glossary.html#term-cigar
		List of cases for which the sequence is removed (all remaining cases and conditions are accepted):
		Given:
			|min_starting_matches_afterNbp|min_starting_matches|rest_of_the_sequence|end_of_the_sequence|
			|---|-----|---|----| sequence:: min_starting_matches_afterNbp:3, min_starting_matches:5, *:{all tags}, O:not M tag, {M,X,H,S,P}:tags, K:{X,P,H}
		Cases to discard the sequence:
			(1)	|***|MOMMM|***|****| -> in the safe zone [min_starting_matches after a number of min_starting_matches_afterNbp bp] the sequence cointains not only M
			(2) |KKK|*****|KKK|****| -> the sequence must not start nor contain X, P or H
			(3) |***|*ID**|ID*|****| -> indels in the no-indel zone (accepted only in the first bases)
			(4) |***|MMMMM|***|SSSS| -> end clipping size > endClipThreshold 

	"""
	valid = True
	# the following values of the cigar flag come from the PYSAM documentation 
	#  http://sqt.googlecode.com/git-history/hashing/sqt/cigar.py
	#  http://www.cgat.org/~andreas/documentation/pysam/api.html

	# revert cigar list if it is in reverse
	strand = '+'
	bam_aln_cigar_list = None
	if bam_aln.is_reverse: # if the read is in reverse, start counting from the end of the CIGAR list
		bam_aln_cigar_list = bam_aln.cigar[::-1]
		strand = '-'
	else:
		bam_aln_cigar_list = bam_aln.cigar
	# get cigar sequence and its list
	bam_aln_cigar_sequence = getCigarSequenceFromCigarList(bam_aln_cigar_list)
	bam_aln_cigar_sequence_list = list(bam_aln_cigar_sequence)
	# now evaluate cigar
	if len(bam_aln_cigar_list) > 1: # if you have an alignment with more than 1 cigar tag
		if len( set(bam_aln_cigar_sequence_list[min_starting_matches_afterNbp:min_starting_matches_afterNbp+min_starting_matches]) - set(['M']) ) > 0:
			print "[AP]\t\t%s -> removed read by CIGAR: case (1) |***|MOMMM|***|\t[chr%d:%d-%s, len %s, CIGAR %s, strand %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring, strand)
			valid = False
		else:
			if bam_aln_cigar_list[-1][0] in [cigar_flag['S']] and int(bam_aln_cigar_list[-1][1]) > endClipThreshold:
				print "[AP]\t\t%s -> removed read by CIGAR: case (4) ends with S soft clipping > %d bp\t[chr%d:%d-%s, len %s, CIGAR %s, strand %s]" %(bam_aln.qname, endClipThreshold, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring, strand)
				valid = False
		if 'X' in bam_aln_cigar_sequence_list[:min_starting_matches_afterNbp] or 'P' in bam_aln_cigar_sequence_list[:min_starting_matches_afterNbp] or 'H' in bam_aln_cigar_sequence_list[:min_starting_matches_afterNbp]:
			print "[AP]\t\t%s -> removed read by CIGAR: case (2) |KKK|*****|***|\t[chr%d:%d-%s, len %s, CIGAR %s, strand %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring, strand)
			valid = False
		if 'X' in bam_aln_cigar_sequence_list[min_starting_matches_afterNbp+min_starting_matches:] or 'P' in bam_aln_cigar_sequence_list[min_starting_matches_afterNbp+min_starting_matches:] or 'H' in bam_aln_cigar_sequence_list[min_starting_matches_afterNbp+min_starting_matches:]:
			print "[AP]\t\t%s -> removed read by CIGAR: case (2) |***|*****|KKK|\t[chr%d:%d-%s, len %s, CIGAR %s, strand %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring, strand)
			valid = False
		if minStartingBasesNoIndels > min_starting_matches_afterNbp+min_starting_matches_afterNbp: # enter here only if the value of noindels is > min_starting_matches_afterNbp+min_starting_matches_afterNbp, otherwise this case is related to the number 1 because indels will be searched in the safe zone of no mismatches. in other words, this setting mantains the tolerance of the min_starting_matches_afterNbp and looks at the indels only in the other bases beyond min_starting_matches_afterNbp+min_starting_matches_afterNbp
			if 'I' in bam_aln_cigar_sequence_list[min_starting_matches_afterNbp+min_starting_matches_afterNbp:minStartingBasesNoIndels] or 'D' in bam_aln_cigar_sequence_list[min_starting_matches_afterNbp+min_starting_matches_afterNbp:minStartingBasesNoIndels]:
				print "[AP]\t\t%s -> removed read by CIGAR: case (3) indels within (<=) %d bp\t[chr%d:%d-%s, len %s, CIGAR %s, strand %s]" %(bam_aln.qname, minStartingBasesNoIndels, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring, strand)
				valid = False
	else: # since you have only 1 tag for the cigar string, orientation (in this case) does not matter
		if (bam_aln_cigar_list[0][0] == cigar_flag['M'] and int(bam_aln_cigar_list[0][1]) >= min_starting_matches_afterNbp + min_starting_matches): # only if you have M -> ok
			valid = True
		else:
			print "[AP]\t\t%s -> removed read by CIGAR: case (1) n!M \t[chr%d:%d-%s, len %s, CIGAR %s, strand %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, bam_aln.cigarstring, strand)
			valid = False
			#print bam_aln.qname, bam_aln_cigar_list, bam_aln.tags
	
	if bam_aln_cigar_list is None:
		print "[AP]\tWarning: NO CIGAR (=None) for read %s" %(bam_aln.qname)
		valid = False
		
	return valid

def evaluateMDtag_v06(bam_aln, min_starting_matches = 3):
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
	strand = '+'
	# find MD tag and string
	mdstring = None
	for tag in bam_aln.tags:
		if tag[0] == 'MD':
			mdstring = str(tag[1])
	# first understand strand specific read orientation
	if bam_aln.is_reverse:
		re_to_use = rev_basesmatch_re
		strand = '-'
	else:
		re_to_use = fwd_basesmatch_re
		strand = '+'
	# only if the MD string is not NULL -> a null string may happen ONLY if the aligner does not give you the MD tag
	if mdstring is not None:
		# execute re finding all occurrences -> it creates an array -> only last/first exact matches. From this number you will classify the read as valid or not
		array_valid_matches = re.findall(re_to_use, mdstring)
		# evaluate returned matches of MD
		if len(array_valid_matches) > 0: # only if there is at least a match -> if not, this means that the read starts with a mismatch
			if int(array_valid_matches[0]) >= min_starting_matches:
				valid = True # return this read as valid
			else:
				print "[AP]\t\t%s -> removed read by MD\t[chr%d:%d-%s, len %s, MD %s, strand %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, mdstring, strand)
		else:
			print "[AP]\t\t%s -> removed read by MD\t[chr%d:%d-%s, len %s, MD %s, strand %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, mdstring, strand)
	else: #
		print "[AP]\t\t%s -> WARNING no MD tag found, this read will be kept anyway because I cannot evaluate the MD score\t[chr%d:%d-%s, len %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen)
		valid = True # return this read as valid
	return valid

def evaluateMDtag(bam_aln, min_starting_matches = 3, min_starting_matches_afterNbp = 0):
	"""
	Input: pysam alignment object
	Output: true or false
	Based on strand, evaluate possible problems in CIGAR string. Given a as HTSeq BAM reader aln,
	a.optional_field('MD')

	Note: rimane il dubbio sulla INSERZIONE che nella conversione in stringa NON viene considerata...! TODO!
	"""
	re_to_use = '([0-9]+)([A-Z]|[a-z]|\^[a-zA-Z]+)*' # find all occurrences of a number paired with a char or a null string and return a list of tuples
	valid = False # checking variable to return
	strand = '+'
	# find MD tag and string
	mdstring = None
	for tag in bam_aln.tags:
		if tag[0] == 'MD':
			mdstring = str(tag[1])
	# only if the MD string is not NULL -> a null string may happen ONLY if the aligner does not give you the MD tag
	if mdstring is not None:
		# execute re finding all occurrences -> it creates an array -> only last/first exact matches. From this number you will classify the read as valid or not
		array_valid_matches = re.findall(re_to_use, mdstring)
		# get the string sequence from the MD tag
		string_sequence = ""
		for n_bp, char_bp in array_valid_matches:
			string_sequence += 'M'*int(n_bp) + 'X'*len(char_bp.replace('^',''))
		# now sort by strand / orientation
		if bam_aln.is_reverse: # revert the array
			strand = '-'
			string_sequence = string_sequence[::-1]
		# evaluate returned matches of MD: 
		if 'X' not in string_sequence[min_starting_matches_afterNbp:min_starting_matches_afterNbp+min_starting_matches]:
			valid = True # return this read as valid
		else:
			print "[AP]\t\t%s -> removed read by MD\t[chr%d:%d-%s, len %s, MD %s, strand %s, seq slice %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, mdstring, strand, string_sequence[min_starting_matches_afterNbp:min_starting_matches_afterNbp+min_starting_matches] )
	else: #
		print "[AP]\t\t%s -> WARNING no MD tag found, this read will be kept anyway because I cannot evaluate the MD score\t[chr%d:%d-%s, len %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen)
		valid = True # return this read as valid
	return valid

def evaluateAlignmentScores(bam_aln, compareSubOptimal, suboptimalThreshold, ASlikeTag, XSlikeTag, SAalnAction):
	"""
	Input: pysam alignment object
	Output: true or false
	Given the alignment scores, if compareSubOptimal=True, from the tags ASlikeTag and XSlikeTag (usually AS and XS tags), compare them creating the delta and verify that the threshold is respected as follows:
		delta=(100 - XS/AS*100)
		if delta >= suboptimalThreshold :: keep the read (read is valid)
		else discard it
	a.optional_field('MD')
	"""
	valid = False
	AS = None # best hit alignment score
	XS = None # suboptimal alignment score
	chimera = False # chimera flag depending by USER option 
	SA = None # SA tag flag (also related to chimera) coming from data
	for tag, val in bam_aln.tags: # acquire data/values from tags
		if tag == ASlikeTag: # best hit alignment AS score
			AS = val
		if tag == XSlikeTag: # XS tag, suboptimal alignment score
			XS = val
		if tag == 'SA' and SAalnAction == 'ignore': # in case of chimeric values, if you want ignore them
			SA = val # assign value to tag only, do not change flag
		elif tag == 'SA' and SAalnAction == 'remove': # in case of chimeric values, if you want to remove them
			chimera = True
			SA = val # assign value to tag only and change flag chimera
		elif tag == 'SA' and SAalnAction == 'removeByComparison': # in case of chimeric values, if you want to remove them by comparing secondary alignments
			chimera = evaluateSAmapq(bam_aln, val, sa_max_mapq = 7, sa_min_delta = 10) # returns True or False based on the evaluation
			SA = val # assign value to tag only and change flag chimera
	if not chimera: # only in case of read that has not chimera (in SAM format is the SA tag)
		if AS > XS: # only if the best hit alignment score is better than subotimal one, go ahead, else this read is not valid
			# find delta: as - (xs/as*100)
			delta = 100.0 - (XS/float(AS)*100)
			# compare delta with threshold
			if delta >= float(suboptimalThreshold):
				valid = True
			else:
				print "[AP]\t\t%s -> removed read by Aln Score (delta<threshold)\t[chr%d:%d-%s, len %s, AS %s, XS %s, SA %s, delta %.2f]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, AS, XS, SA, delta)
		else: # in the case in which AS <= XS, exclude the read a priori
			print "[AP]\t\t%s -> removed read by Aln Score (AS<=XS)\t[chr%d:%d-%s, len %s, AS %s, XS %s, SA %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, AS, XS, SA)	
	else:
		print "[AP]\t\t%s -> removed read by Chimera evaluation\t[chr%d:%d-%s, len %s, SA %s]" %(bam_aln.qname, int(bam_aln.tid)+1, bam_aln.pos, bam_aln.aend, bam_aln.alen, SA)	
	return valid

def evaluateSAmapq(bam_aln, sa_string, sa_max_mapq = 7, sa_min_delta = 10):
	"""
	Input: only the SA tag string (e.g. chr6,142629962,+,66S24M109S,20,0;chr18,54898335,+,152S22M25S,0,0;chr1,20664235,+,141S21M37S,9,0;chr19,9056388,+,46S19M134S,0,0;) , in the format <chr,start,strand,cigar,mapq,nm>
	Output: bool, true if this read is not chimeric
	Logics:
		sa_max_mapq = 10, sa_min_delta = 10. -> if the input read has at least one SA that overcomes mapq (>=10) or delta mapq (<=10), then this read has got a chimera -> not valid
	"""
	valid = False
	chimera_counter = 0
	for x in sa_string.split(';'):
		if x != '':
			chr,start,strand,cigar,mapq,nm = x.split(',')
			if len(chr.split('_'))<2: # no chimera on random chrs
				if int(mapq) >= sa_max_mapq or (abs(int(bam_aln.mapq) - int(mapq)) <= sa_min_delta):
					chimera_counter += 1
	# now assess delta from counter
	if chimera_counter == 0:
		valid = True
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

def worker(todo_regions, bam, minStartingMatches, minStartingMatches_fromNbp, minStartingBasesNoIndels, endClipThreshold, pruneByPair, out_q, single_end, compareSubOptimal, suboptimalThreshold, ASlikeTag, XSlikeTag, SAalnAction):
	""" 
	The worker function, invoked in a process. 
		region: array of tuples with coordinate to slice chromosome
		bam: the input bam to splice
		out_q: queue of the process
	"""
	outdata = []
	for region in todo_regions:
		chr, start, end = region
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t...processing region (chr, start, end):", chr, start, end
		for bam_alingment in mybam:
			if not single_end: # ILLUMINA paired-end reads
				query_aln = bam_alingment.qname.strip()
				# which pair is the driver?
				if pruneByPair == 'both':
					#print "query:", query_aln
					if compareSubOptimal: # if you WANT to evaluate alignment scores
						valid_alnscore = evaluateAlignmentScores(bam_alingment, compareSubOptimal, suboptimalThreshold, ASlikeTag, XSlikeTag, SAalnAction) # evaluate alignment score AS and XS (best vs suboptimal score)
						if valid_alnscore: # if the read passes the alignment check, go ahead
							valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
							valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
							if not (valid_cigar and valid_MDtag):
								# report data into an array
								outdata.append(query_aln)
						else: # do not keep it -> append results to output
							outdata.append(query_aln)
					else: # if you DO NOT want to evaluate alignment scores
						valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
						valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
						if not (valid_cigar and valid_MDtag):
							# report data into an array
							outdata.append(query_aln)
				elif pruneByPair == 'read1':
					if bam_alingment.is_read1:
						if compareSubOptimal: # if you WANT to evaluate alignment scores
							valid_alnscore = evaluateAlignmentScores(bam_alingment, compareSubOptimal, suboptimalThreshold, ASlikeTag, XSlikeTag, SAalnAction) # evaluate alignment score AS and XS (best vs suboptimal score)
							if valid_alnscore: # if the read passes the alignment check, go ahead
								valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
								valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
								if not (valid_cigar and valid_MDtag):
									# report data into an array
									outdata.append(query_aln)
							else: # do not keep it -> append results to output
								outdata.append(query_aln)
						else: # if you DO NOT want to evaluate alignment scores
							valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
							valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
							if not (valid_cigar and valid_MDtag):
								# report data into an array
								outdata.append(query_aln)
				elif pruneByPair == 'read2':
					if bam_alingment.is_read2:
						if compareSubOptimal: # if you WANT to evaluate alignment scores
							valid_alnscore = evaluateAlignmentScores(bam_alingment, compareSubOptimal, suboptimalThreshold, ASlikeTag, XSlikeTag, SAalnAction) # evaluate alignment score AS and XS (best vs suboptimal score)
							if valid_alnscore: # if the read passes the alignment check, go ahead
								valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
								valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
								if not (valid_cigar and valid_MDtag):
									# report data into an array
									outdata.append(query_aln)
							else: # do not keep it -> append results to output
								outdata.append(query_aln)
						else: # if you DO NOT want to evaluate alignment scores
							valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
							valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
							if not (valid_cigar and valid_MDtag):
								# report data into an array
								outdata.append(query_aln)
			else: # 454/single end reads
				#print chr, aln.qname
				#print bam_alingment.qname, bam_alingment.is_reverse
				query_aln = bam_alingment.qname.strip()
				if compareSubOptimal: # if you WANT to evaluate alignment scores
					valid_alnscore = evaluateAlignmentScores(bam_alingment, compareSubOptimal, suboptimalThreshold, ASlikeTag, XSlikeTag, SAalnAction) # evaluate alignment score AS and XS (best vs suboptimal score)
					if valid_alnscore: # if the read passes the alignment check, go ahead
						valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
						valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
						#print "valid_cigar =", valid_cigar, "; valid_MDtag =", valid_MDtag
						if not (valid_cigar and valid_MDtag):
							outdata.append(query_aln)
					else: # do not keep it -> append results to output
						#print "remove it, not valid_alnscore"
						outdata.append(query_aln)
				else:
					valid_cigar = evaluateCigar(bam_alingment, min_starting_matches = minStartingMatches, minStartingBasesNoIndels = minStartingBasesNoIndels, endClipThreshold = endClipThreshold, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if M -> ok, else discard it
					valid_MDtag = evaluateMDtag(bam_alingment, min_starting_matches = minStartingMatches, min_starting_matches_afterNbp = minStartingMatches_fromNbp) # if MD tag has correct match -> ok, else discard it ## TODO!
					if not (valid_cigar and valid_MDtag):
						outdata.append(query_aln)

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
	# first check args and file paths
	print "[AP]\tYour settings:"
	for k, v in vars(args).iteritems():
		print "\t\t%s::\t\t%s" %(k, v)

	print "[AP]\tChecking inputs."
	checkArgs(args)
	# acquire thresholds
	minStartingMatches = int(args.minStartingMatches)
	minStartingMatches_fromNbp = int(args.minStartingMatches_fromNbp)
	minStartingBasesNoIndels = int(args.minStartingBasesNoIndels)
	endClipThreshold = int(args.endClipThreshold)
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
		p = multiprocessing.Process(target=worker, args=(regions[chunksize * i:chunksize * (i + 1)], bam, minStartingMatches, minStartingMatches_fromNbp, minStartingBasesNoIndels, endClipThreshold, args.pruneByPair, out_q, args.single_end, args.compareSubOptimal, args.suboptimalThreshold, args.ASlikeTag, args.XSlikeTag, args.SAalnAction )) # prototype: todo_regions, bam, bed_dictionary, minStartingMatches, out_q
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
