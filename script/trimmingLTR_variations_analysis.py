#!/usr/bin/python
# -*- coding: utf-8 -*-
import MySQLdb, sys, os, re, argparse, csv, time
import pysam
from operator import itemgetter, attrgetter
from string import Template
#from Bio import SeqIO
from time import gmtime, strftime, sleep
from multiprocessing import Process, Lock, Queue
import multiprocessing, math

"""
This programs must be run after the following script that aligns reads to target regions.

#!/bin/bash

CHANGEDIR="/home/andrea/test/LR40"
VECTORLTR="/opt/genome/vector/ltr/LTR.32bp.fa"
PREFIX="LR40"
R1_FASTQ="/storage/dx/backup/nas/LabDevelopment/HIV_patients/data/raw/20141003_OSR_MiSeq_LR40/LR40_S1_L001_R1_001.fastq.gz"
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;

### plot raw data stats
cd ${CHANGEDIR}
fastx_quality_stats -i <(zcat ${R1_FASTQ}) -o ${BNAME_R1}.statslog -Q 33
fastq_quality_boxplot_graph.sh -i ${BNAME_R1}.statslog -o ${BNAME_R1}.perbasequality.png -t "${BNAME_R1}"
fastx_nucleotide_distribution_graph.sh -i ${BNAME_R1}.statslog -o ${BNAME_R1}.nucleocomposition.png -t "${BNAME_R1}"


### plot LTR stats after alignment
# exmaple from my home, testing LR40
cd ${CHANGEDIR}

# align sequences to LTR genome VECTORLTR
bwa-7.5 mem -r 1 -T 15 -R "@RG\tID:${PREFIX}.RawLTR\tSM:LTR32bp\tCN:${PREFIX}.RawDataR1" -t 16 ${VECTORLTR} <(zcat ${R1_FASTQ} ) > ${PREFIX}.LTRgenome.rawR1.alm.sam

# remove reads that are:
#- not primary
#- not mapped
#- in reverse orientation
#- under quality alignment of 3 (phred)
samtools view -F 276 -q 3 -uS ${PREFIX}.LTRgenome.rawR1.alm.sam | samtools sort - ${PREFIX}.LTRgenome.rawR1.alm.sorted;
# index bam and do the MD filling
samtools index ${PREFIX}.LTRgenome.rawR1.alm.sorted.bam ;
samtools fillmd -b ${PREFIX}.LTRgenome.rawR1.alm.sorted.bam ${VECTORLTR} > ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.bam

# count per base (nucleotide) chenge and frequency (from a piled-up alignment)
# sembra funzionare ma NON benissimo....
samtools mpileup -d 1000000000 -f /opt/genome/vector/ltr/LTR.32bp.fa ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.bam > ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.pileup
/home/andrea/test/LR40/./sequenza-utils.py pileup2acgt ${PREFIX}.LTRgenome.rawR1.alm.sorted.md.pileup



-----------
Poi lanci questo dal DB

SELECT `r1_RG` , `r1_chr` , `r1_offtarget_type` , `r1_offtarget_interval` , `r1_offtarget_size` , `r2_offtarget_type` , `r2_offtarget_interval` , count( * ) AS groupreads, `r2_offtarget_size`
FROM `zinc_whole_t6`
WHERE `r1_NM` <=10 
GROUP BY `r1_RG` , `r1_chr` , `r1_offtarget_type` , `r1_offtarget_interval` , `r2_offtarget_type` , `r2_offtarget_interval`
ORDER BY `r1_RG` , `r1_chr` , `r1_offtarget_type` , `r2_offtarget_type` , groupreads

"""

header = """
+--------------------------------------------------+
  ***   ANALYZE LTR VARIATIONS FROM BAM    ***
+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		November 2014
 Contact:	andrea.calabria@hsr.it
 Revision:	0.1 (MP)
+--------------------------------------------------+
  
 Note:
  - it uses PySAM
  - BAM file will be indexed (otherwise PySAM will not work)
  - it creates a TSV file that will be imported in MySQL
  - target file format (TSV):
  <targetId> <targetNumber_0based> <start> <end>
  - from the associationfile it takes the last two columns (as -1 and -2)
  
 TODO:
  - SQLite access -> avoid MySQL
 
""" 

description = "This application will import ZF targets from BAM files."

usage_example = """
Examples of usage (in folder ALIEN ~/Dropbox/TIGET/Workspace/Andrea/ZincFTargets [http://172.25.39.2/workbench/ZincFTarget/test/bam/]):
 (1) APP --bam whole.experiment.sorted.md.rel.bam -p 10 -r targets.intervals --dbhost locahost --dbschema test --dbtable zinc.whole.t6 --dbuser andrea --dbpassword andrea

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--bam', dest="bamfile", help="BAM file to process. No default option.", action="store", required=True)
parser.add_argument('--sortBAM', dest="sortBAM", action="store_true", help="Sort BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--indexBAM', dest="indexBAM", action="store_true", help="Index BAM file (calling SAMtools). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('-o', '--outfilename', dest="outfilename", help="Output file name of the quantified targets. No default option.", action="store", required=True)
parser.add_argument('-p', '--processes', dest="processes", help="Number of processes to use. Default 2.", action="store", default=2)
parser.add_argument('-r', '--regions', dest="regions", help="Acquire the Regions to process from a CSV file in the standard (short) BED format: <targetId> <start> <end> <reference region lenght>, where targetId is the reference target ID (ideally the chromosome ID). No Default value; required.", action="store", required=True)
parser.add_argument('--bpoffset', dest="bpoffset", help="Number of bases to add before and after regions for each target. Default = 0.", action="store", default=0)
parser.add_argument('--maxNM', dest="maxNM", help="Number of maximum allowed mismatches. NB: deletions are accounted as mismatches! Default = 10.", action="store", default=10)
parser.add_argument('--dbhost', dest="dbhost", help="Target Database host in which importing data. E.g.: localhost. default=localhost.", action="store", default="localhost")
parser.add_argument('--dbport', dest="dbport", help="Target Database port in which importing data. E.g.: 3306. default=3306.", action="store", default="3306")
parser.add_argument('--dbschema', dest="dbschema", help="Target Database schema in which importing data. E.g.: MLD01. Mandatory.", action="store", required=True)
parser.add_argument('--dbtable', dest="dbtable", help="Target Database table in which importing data. E.g.: redundant_MLD01_ALL_NEW_NAMES. Mandatory.", action="store", required=True)
parser.add_argument('--dbuser', dest="dbuser", help="Target Database user. Mandatory.", action="store", required=True)
parser.add_argument('--dbpassword', dest="dbpassword", help="Target Database password. Mandatory.", action="store", required=True)
parser.add_argument('--controlSample', dest="controlSample", help="Control sample ID (from the RG groups) that will be the lowest abundance level: this sample ID will discriminate background noise from true positive values. Default is 'NOT SET', otherwise write the corresponding string of the RG from your BAM file group that will be treated as control (case sensitive).", action="store", default='NOT SET')
parser.add_argument('--targets', dest="targetSequenceFile", help="The target file of tab separated values: region ID (for example gene name or chr), reference sequences.", action="store", default="NOT SET")
parser.add_argument('--fastaFolder', dest="fastaFolder", help="The target FASTA files folder that will contain fasta files of splitted regions with MSA results. Default: current directory", action="store", default=os.getcwd())
parser.add_argument('--doClustering', dest="doClustering", help="To run the MSA with output plots and sequences, activate this option.", action="store_true", default=False)
args = parser.parse_args()

#########################################################################################
####### GLOBAL VARS
#########################################################################################

# init db values -> from UCSC database ChromInfo
hg19="1:1-249250621;2:1-243199373;3:1-198022430;4:1-191154276;5:1-180915260;6:1-171115067;7:1-159138663;8:1-146364022;9:1-141213431;10:1-135534747;11:1-135006516;12:1-133851895;13:1-115169878;14:1-107349540;15:1-102531392;16:1-90354753;17:1-81195210;18:1-78077248;19:1-59128983;20:1-63025520;21:1-48129895;22:1-51304566;X:1-155270560;Y:1-59373566;M:1-16571"
mm9="1:1-197195432;2:1-181748087;3:1-159599783;4:1-155630120;5:1-152537259;6:1-149517037;7:1-152524553;8:1-131738871;9:1-124076172;10:1-129993255;11:1-121843856;12:1-121257530;13:1-120284312;14:1-125194864;15:1-103494974;16:1-98319150;17:1-95272651;18:1-90772031;19:1-61342430;X:1-166650296;Y:1-15902555;M:1-16299"

readorder = ['RG', 'chr', 'start', 'end', 'strand', 'quality', 'NM', 'flag', 'cigar', 'cigar_components', 'offtarget', 'offtarget_type', 'offtarget_interval', 'offtarget_size', 'intarget_seq', 'MD', 'insert_size', 'AS', 'XS', 'SA'] # this is the order of the content to written in the DB for each read (and pair)
fields = ["sample", "target", "header", "var_type", "var_start", "var_end"]
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
	'C': 9, # Andrea Defined as BASE CHANGE
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
	9: 'C', # Andrea Defined as BASE CHANGE
	} 

field_translation = {
	'r1_chr': 'Region', 
	'r1_RG': 'Sample',
	'r1_offtarget_type': 'Variation Type', 
	'r1_offtarget_interval': 'Variation Interval', 
	'r1_offtarget_size': 'Variation Size', 
	'groupreads': 'Reads Count',
	'ratio': 'Variations vs Others',
	'percentage': 'Abundance Variation percentage',
	'ref_ratio': 'Other sites reads count',
	'ref_all': 'Overall Reads',
	}

string_where_r1_NM = "r1_NM <= %s" %(args.maxNM)

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
	if args.targetSequenceFile not in ['NOT SET']:
		if not os.path.isfile(args.targetSequenceFile):
			print "\n[AP]\tError while reading file of targets: no valid paths.\n\tInput files: %s\n\n\tExiting...\n" %(args.targets)
	if not os.path.isdir(args.fastaFolder):
		print "\n[AP]\tError while reading options: no valid paths.\n\tFolder FASTA: %s\n\n\tExiting...\n" %(args.fastaFolder)

def parseTargetFile(infile, bpoffset):
	"""
	IN: association file, BED file format
	OUT: array of regions as list of tuples as target:start-end
	FILE FORMAT:   <targetId> <start> <end> <ref len>
	"""
	targetlist = []
	with open(infile, 'rb') as inf:
		reader = csv.reader(inf, delimiter = "\t")
		for row in reader:
			if len(row) > 0 and not row[0].startswith('#'):
				target_start = int(row[1].strip()) - bpoffset
				target_end = int(row[2].strip()) + bpoffset
				ref_size = int(row[3].strip())
				if target_start < 0 and target_end > ref_size:
					print "[AP]\tERROR: wrong offset for target region, too large -> remaining starting position <0 and end is over region limit.\n"
					sys.exit()
				targetlist.append( (row[0].strip(),  target_start, target_end) )
			else:
				print "[AP]\tWarning: this file contais a blank row!!!! Check it please."
	targetdictionary = {} # k = region, v = [start, end]
	for (region, start, end) in targetlist:
		targetdictionary[region] = [start, end]
	return targetlist, targetdictionary

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
	os.system( "samtools sort %s %s" %(bamfile, sorted_bamfile) )
	return sorted_bamfile + ".bam"

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

def getCigarByOrientation(bam_aln, ):
	"""
	Given the Pysam codes for CIGAR string conversion, return the CIGAR list by orientation. The output will be the CIGAR list in 5'-3' always, so that all counts in the interval will be evaluated in the same orientation.
	"""
	tagval = None
	if not bam_aln.is_reverse: # if 5' - 3'
		tagval = bam_aln.cigar
	else: # reverse case, take the last occurrence -> revert list, 
		tagval = bam_aln.cigar[::-1]
	return tagval 

def worker(todo_regions, bam, out_q, ): 
	""" 
	The worker function, invoked in a process. 
		region: array of tuples with coordinate to slice chromosome
		bam: the input bam to splice
		out_q: queue of the process
	"""
	outdata = []
	for region in todo_regions:
		chr, start, end = region # now start and end are the boundaries of putative targets!
		product = {} # k = header, v = dict:: "sample", "target", "header", "var_type", "var_start", "var_end"
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t...processing region (chr, start, end)", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		
		print "[AP]\t\t\tAcquiring dictionary of reads in the region (chr, start, enc)\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		for bam_alignment in mybam: # for each alignment in the BAM file (both R1 and R2!!!)
			# get orientation and check for correct position (get end from orientation and start)
			orientation = '+' # default, if fwd
			aln_start = bam_alignment.pos +1 # default, if fwd
			aln_end = bam_alignment.pos + bam_alignment.alen # default, if fwd
			if bam_alignment.is_reverse:
				orientation = '-'
				# aln_end = bam_alignment.pos
				# aln_start = bam_alignment.pos + bam_alignment.alen

			# get the CIGAR list
			## cigarlist = getCigarByOrientation(bam_alignment) # this is true if you revert the start!!! otherwise do not use it!
			cigarlist = bam_alignment.cigar

			# acquire tags
			AS = "NULL" # best hit alignment score
			XS = "NULL" # suboptimal alignment score
			SA = "NULL" # secondary alignments / chimera
			NM = "NULL" # number of mismatches
			MD = "NULL" # md score from cigar
			RG = "NULL" # sample (here read group)
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

			# get the offtarget definition/evaluation
			component_start = aln_start
			#component_end = 0
			# scorri in ordine di mappatura il cigar per frammentare la regione nelle componenti
			# note: in caso di INSERZIONI non spostare in avanti il contatore delle basi in quanto non le vedrai; mentre nelle delezioni lo devi spostare (la logica e' inversa a quella comune di calcolo, ovvero se inserisco aggiungo mentre se deleto tolgo)
			components = [] # list of tuples
			for (component, component_size) in cigarlist:
				component_string = cigar_number[component]
				if component_start == aln_start: # if you are in the first component
					if component_string not in ['I', 'S']: # if the read is aligned and DOES NOT start with S, if the component is NOT an insertion, then increment start and end value
						components.append((component_string, int(component_start), component_start + int(component_size)-1, int(component_size)))
						component_start += int(component_size)
					else: # if the component is an insertion, then DO NOT increment start and end value (do not move over the genomic interval)
						components.append((component_string, int(component_start), int(component_start), int(component_size)))
				else: # for all the other components after the first one
					if component_string is not 'I': # if the component is NOT an insertion, then increment start and end value
						components.append((component_string, int(component_start), component_start + int(component_size)-1, int(component_size)))
						component_start += int(component_size)
					else: # if the component is an insertion, then DO NOT increment start and end value (do not move over the genomic interval)
						components.append((component_string, int(component_start), int(component_start), int(component_size)))
			
			# valuta se e' off target valido!! occhio: off target validi sono tutti quelli che sia sono completamente inclusi e se hanno delezioni che iniziano dentro il confine ma terminano oltre
			offtarget = False # initialize offtarget var
			offtarget_type_string = None # it will be the I or D (or in case of multiple types the concatenation of them)
			offtarget_type_list = []
			offtarget_interval_string = None # the interval of the offtarget region
			offtarget_interval_list = []
			offtarget_size_string = None
			offtarget_size_list = []
			for (component_string, component_start, component_end, component_size) in components:
				if component_string in ['D', 'I']: # it means in D or I (deletion or insertions)
					# casi: (1) eventi completamente contenuti nella regione target, (2) eventi che iniziano nella regione ma finiscono oltre l'end, (3) iniziano prima dello start e finiscono dentro la regione
					if (component_start >= start and component_end <= end) or (component_start >= start and component_start <= end and component_end >= end) or (component_end >= start and component_end <= end and component_start <= start) or (component_start <= start and component_end >= end) : # where start and and end are already cominig from the target region!!!
						offtarget = True
						offtarget_type_list.append(component_string)
						offtarget_interval_list.append(str(component_start) +"-"+ str(component_end))
						offtarget_size_list.append(str(component_size))
			offtarget_type_string = '|'.join(x for x in offtarget_type_list)
			offtarget_interval_string = '|'.join(x for x in offtarget_interval_list)
			offtarget_size_string = '|'.join(x for x in offtarget_size_list)
			if offtarget_type_string == "":
				offtarget_type_string = "NULL"
			if offtarget_interval_string == "":
				offtarget_interval_string = "NULL"
			if offtarget_size_string == "":
				offtarget_size_string = "NULL"

			# get the in-target dna sequence from the alignment -> only if the sequence starts with a S we MUST shift the aln start position! otherwise you will miss the target region!
			aln_sequence = None
			component_index = 0
			offtarget_sequence_start = 0
			offtarget_sequence_end = 0
			for (component_string, component_start, component_end, component_size) in components: # the first component
				if component_index == 0:
					if component_string not in ['S']: # if not Soft clipped
						#aln_sequence = bam_alignment.seq[start-bam_alignment.pos:end-bam_alignment.pos]
						offtarget_sequence_start = start-bam_alignment.pos
						offtarget_sequence_end = end-bam_alignment.pos
					else:
						#aln_sequence = bam_alignment.seq[start-bam_alignment.pos+component_size:end-bam_alignment.pos+component_size]
						offtarget_sequence_start = start-bam_alignment.pos+component_size
						offtarget_sequence_end = end-bam_alignment.pos+component_size
				else: # component index
					if component_end <= start:# if something deleted/inserted elsewhere, add/remove it. se la fine di questa componente e' entro l'inizio della target region, allora devi guardarla (se non entra in questo IF allora e' un target component vero, quindi non devo prenderlo!!!)
						if component_string in ['D']: # with a deletion before the target region you need to add the same interval len (=> component_size)
							offtarget_sequence_start -= component_size
							offtarget_sequence_end -= component_size
						if component_string in ['I']: # with an insertion before the target region you need to remove the same interval len (=> component_size)
							offtarget_sequence_start += component_size
							offtarget_sequence_end += component_size
						if component_string in ['M']: # with an insertion before the target region you need to remove the same interval len (=> component_size)
							offtarget_sequence_start += component_size
							offtarget_sequence_end += component_size
				component_index += 1
			aln_sequence = bam_alignment.seq[offtarget_sequence_start:offtarget_sequence_end]

			# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
			# 	print bam_alignment.qname, "is R2?", bam_alignment.is_read2

			# now acquire the DB of reads
			if product.has_key(bam_alignment.qname):
				## -- debug -- start --
				# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
				# 	print bam_alignment.qname, "in IF key exists, is R2?", bam_alignment.is_read2
				## -- debug -- end --
				if not bam_alignment.is_read2: # all but not read 2: this means read 1 and potentially all singletones from r1
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in IF key exists and IF NOT R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
					product[bam_alignment.qname]['r1'] = {	'start': aln_start, 
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
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
													}
				else:
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in IF key exists and ELSE R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
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
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
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
															'AS': AS,
															'XS': XS,
															'SA': SA,
															'NM': NM,
															'MD': MD,
															'RG': RG,
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
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
															'AS': AS,
															'XS': XS,
															'SA': SA,
															'NM': NM,
															'MD': MD,
															'RG': RG,
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
														}
													}

		# now loop over the dictionary keys to create the output array based on the driver read
		# forma dell'output: 'prod_header', 'prod_chr', 'prod_locus', 'prod_end', 'prod_strand', 'ref_associationid', 'ref_matrixid', 'ref_poolid', 'isread_start', 'isread_end', 'isread_strand', 'isread_chr', 'isread_cigar', 'isread_flag', 'isread_insert_size', 'isread_AS', 'isread_XS', 'isread_SA', 'isread_NM', 'isread_MD', 'mate_start', 'mate_end', 'mate_strand', 'mate_chr', 'mate_cigar', 'mate_flag', 'mate_insert_size', 'mate_AS', 'mate_XS', 'mate_SA', 'mate_NM', 'mate_MD',

		locus = 'r1'
		mate = 'r2'
		for k, v in product.iteritems():

			# if k in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
			# 	print k, v

			is_array = [k, ] # data array
			# R1
			if v.has_key(locus):
				for x in readorder: # first the locus read
					is_array.append( str(v[locus][x]) )
			else:
				for x in readorder: 
					is_array.append( 'NULL' )
			# R2
			if v.has_key(mate):# then the mate read, if it exists
				for x in readorder: 
					is_array.append( str(v[mate][x]) )
			else:
				for x in readorder: 
					is_array.append( 'NULL' )
					####### -> ???? is_array.append( str('NULL') for x in range(0,len(readorder)))
			outdata.append(is_array) # append results to output

		print "[AP]\t\t\tDone, region completed (chr, start, end)\t\t\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" 
	
	#print len(outdata), outdata[0]
	# queue results at the end of this process
	out_q.put(outdata)

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

def singleWorker(todo_regions, bam, out_q, ): 
	""" 
	The worker function, invoked in a process. 
		region: array of tuples with coordinate to slice chromosome
		bam: the input bam to splice
		out_q: queue of the process
	"""
	outdata = []
	for region in todo_regions:
		chr, start, end = region # now start and end are the boundaries of putative targets!
		product = {} # k = header, v = dict:: "sample", "target", "header", "var_type", "var_start", "var_end"
		mybam = bam.fetch(chr, start, end)
		print "[AP]\t\t...processing region (chr, start, end)", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		
		print "[AP]\t\t\tAcquiring dictionary of reads in the region (chr, start, enc)\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" %(mybam.mapped, mybam.unmapped, )
		for bam_alignment in mybam: # for each alignment in the BAM file (both R1 and R2!!!)
			# get orientation and check for correct position (get end from orientation and start)
			orientation = '+' # default, if fwd
			aln_start = bam_alignment.pos +1 # default, if fwd
			aln_end = bam_alignment.pos + bam_alignment.alen # default, if fwd
			if bam_alignment.is_reverse:
				orientation = '-'
				# aln_end = bam_alignment.pos
				# aln_start = bam_alignment.pos + bam_alignment.alen

			# get the CIGAR list
			## cigarlist = getCigarByOrientation(bam_alignment) # this is true if you revert the start!!! otherwise do not use it!
			cigarlist = bam_alignment.cigar

			# acquire tags
			AS = "NULL" # best hit alignment score
			XS = "NULL" # suboptimal alignment score
			SA = "NULL" # secondary alignments / chimera
			NM = "NULL" # number of mismatches
			MD = "NULL" # md score from cigar
			RG = "NULL" # sample (here read group)
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

			# get the offtarget definition/evaluation
			component_start = aln_start
			#component_end = 0
			# scorri in ordine di mappatura il cigar per frammentare la regione nelle componenti
			# note: in caso di INSERZIONI non spostare in avanti il contatore delle basi in quanto non le vedrai; mentre nelle delezioni lo devi spostare (la logica e' inversa a quella comune di calcolo, ovvero se inserisco aggiungo mentre se deleto tolgo)
			components = [] # list of tuples
			for (component, component_size) in cigarlist:
				component_string = cigar_number[component]
				if component_start == aln_start: # if you are in the first component
					if component_string not in ['I', 'S']: # if the read is aligned and DOES NOT start with S, if the component is NOT an insertion, then increment start and end value
						components.append((component_string, int(component_start), component_start + int(component_size)-1, int(component_size)))
						component_start += int(component_size)
					else: # if the component is an insertion, then DO NOT increment start and end value (do not move over the genomic interval)
						components.append((component_string, int(component_start), int(component_start), int(component_size)))
				else: # for all the other components after the first one
					if component_string is not 'I': # if the component is NOT an insertion, then increment start and end value
						# se hai M, valuta se hai cambi di base specifici, quindi se M ha MD tag che non e' un numero intero, allora splitta l'M tag nelle componenti dell'MD
						components.append((component_string, int(component_start), component_start + int(component_size)-1, int(component_size)))
						component_start += int(component_size)
					else: # if the component is an insertion, then DO NOT increment start and end value (do not move over the genomic interval)
						components.append((component_string, int(component_start), int(component_start), int(component_size)))
			
			# valuta se e' off target valido!! occhio: off target validi sono tutti quelli che sia sono completamente inclusi e se hanno delezioni che iniziano dentro il confine ma terminano oltre
			offtarget = False # initialize offtarget var
			offtarget_type_string = None # it will be the I or D (or in case of multiple types the concatenation of them)
			offtarget_type_list = []
			offtarget_interval_string = None # the interval of the offtarget region
			offtarget_interval_list = []
			offtarget_size_string = None
			offtarget_size_list = []
			for (component_string, component_start, component_end, component_size) in components:
				if component_string in ['D', 'I'] or (component_string in ['M'] and component_end-component_start < end-start): # it means in D or I (deletion or insertions)
					# casi: (1) eventi completamente contenuti nella regione target, (2) eventi che iniziano nella regione ma finiscono oltre l'end, (3) iniziano prima dello start e finiscono dentro la regione
					if (component_start >= start and component_end <= end) or (component_start >= start and component_start <= end and component_end >= end) or (component_end >= start and component_end <= end and component_start <= start) or (component_start <= start and component_end >= end) : # where start and and end are already cominig from the target region!!!
						offtarget = True
						offtarget_type_list.append(component_string)
						offtarget_interval_list.append(str(component_start) +"-"+ str(component_end))
						offtarget_size_list.append(str(component_size))
			offtarget_type_string = '|'.join(x for x in offtarget_type_list)
			offtarget_interval_string = '|'.join(x for x in offtarget_interval_list)
			offtarget_size_string = '|'.join(x for x in offtarget_size_list)
			if offtarget_type_string == "":
				offtarget_type_string = "NULL"
			if offtarget_interval_string == "":
				offtarget_interval_string = "NULL"
			if offtarget_size_string == "":
				offtarget_size_string = "NULL"

			# get the in-target dna sequence from the alignment -> only if the sequence starts with a S we MUST shift the aln start position! otherwise you will miss the target region!
			aln_sequence = None
			component_index = 0
			offtarget_sequence_start = 0
			offtarget_sequence_end = 0
			for (component_string, component_start, component_end, component_size) in components: # the first component
				if component_index == 0:
					if component_string not in ['S']: # if not Soft clipped
						#aln_sequence = bam_alignment.seq[start-bam_alignment.pos:end-bam_alignment.pos]
						offtarget_sequence_start = start-bam_alignment.pos
						offtarget_sequence_end = end-bam_alignment.pos
					else:
						#aln_sequence = bam_alignment.seq[start-bam_alignment.pos+component_size:end-bam_alignment.pos+component_size]
						offtarget_sequence_start = start-bam_alignment.pos+component_size
						offtarget_sequence_end = end-bam_alignment.pos+component_size
				else: # component index
					if component_end <= start:# if something deleted/inserted elsewhere, add/remove it. se la fine di questa componente e' entro l'inizio della target region, allora devi guardarla (se non entra in questo IF allora e' un target component vero, quindi non devo prenderlo!!!)
						if component_string in ['D']: # with a deletion before the target region you need to add the same interval len (=> component_size)
							offtarget_sequence_start -= component_size
							offtarget_sequence_end -= component_size
						if component_string in ['I']: # with an insertion before the target region you need to remove the same interval len (=> component_size)
							offtarget_sequence_start += component_size
							offtarget_sequence_end += component_size
						if component_string in ['M']: # with an insertion before the target region you need to remove the same interval len (=> component_size)
							offtarget_sequence_start += component_size
							offtarget_sequence_end += component_size
				component_index += 1
			aln_sequence = bam_alignment.seq[offtarget_sequence_start:offtarget_sequence_end]

			# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
			# 	print bam_alignment.qname, "is R2?", bam_alignment.is_read2

			# now acquire the DB of reads
			if product.has_key(bam_alignment.qname):
				## -- debug -- start --
				# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
				# 	print bam_alignment.qname, "in IF key exists, is R2?", bam_alignment.is_read2
				## -- debug -- end --
				if not bam_alignment.is_read2: # all but not read 2: this means read 1 and potentially all singletones from r1
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in IF key exists and IF NOT R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
					product[bam_alignment.qname]['r1'] = {	'start': aln_start, 
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
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
													}
				else:
					## -- debug -- start --
					# if bam_alignment.qname in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
					# 	print bam_alignment.qname, "in IF key exists and ELSE R2, is R2?", bam_alignment.is_read2
					## -- debug -- end --
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
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
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
															'AS': AS,
															'XS': XS,
															'SA': SA,
															'NM': NM,
															'MD': MD,
															'RG': RG,
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
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
															'AS': AS,
															'XS': XS,
															'SA': SA,
															'NM': NM,
															'MD': MD,
															'RG': RG,
															'chr': bam.getrname(bam_alignment.tid), 
															'cigar': bam_alignment.cigarstring,
															'cigar_components': components,
															'offtarget': offtarget,
															'offtarget_type': offtarget_type_string,
															'offtarget_interval': offtarget_interval_string,
															'offtarget_size': offtarget_size_string,
															'flag': bam_alignment.flag,
															'insert_size': bam_alignment.tlen ,
															'quality': bam_alignment.mapq,
															'intarget_seq': aln_sequence,
														}
													}

		# now loop over the dictionary keys to create the output array based on the driver read
		# forma dell'output: 'prod_header', 'prod_chr', 'prod_locus', 'prod_end', 'prod_strand', 'ref_associationid', 'ref_matrixid', 'ref_poolid', 'isread_start', 'isread_end', 'isread_strand', 'isread_chr', 'isread_cigar', 'isread_flag', 'isread_insert_size', 'isread_AS', 'isread_XS', 'isread_SA', 'isread_NM', 'isread_MD', 'mate_start', 'mate_end', 'mate_strand', 'mate_chr', 'mate_cigar', 'mate_flag', 'mate_insert_size', 'mate_AS', 'mate_XS', 'mate_SA', 'mate_NM', 'mate_MD',

		locus = 'r1'
		mate = 'r2'
		for k, v in product.iteritems():

			# if k in ['M00571:54:000000000-A6H4T:1:2105:8010:15099', 'M00571:54:000000000-A6H4T:1:2105:13972:17573']:
			# 	print k, v

			is_array = [k, ] # data array
			# R1
			if v.has_key(locus):
				for x in readorder: # first the locus read
					is_array.append( str(v[locus][x]) )
			else:
				for x in readorder: 
					is_array.append( 'NULL' )
			# R2
			if v.has_key(mate):# then the mate read, if it exists
				for x in readorder: 
					is_array.append( str(v[mate][x]) )
			else:
				for x in readorder: 
					is_array.append( 'NULL' )
					####### -> ???? is_array.append( str('NULL') for x in range(0,len(readorder)))
			outdata.append(is_array) # append results to output

		print "[AP]\t\t\tDone, region completed (chr, start, end)\t\t\t", chr, start, end #, " that contains:\n\t\t\t%d mapped reads\n\t\t\t%d unmapped reads" 
	
	#print len(outdata), outdata[0]
	# queue results at the end of this process
	return outdata

def importData(dbhost, dbport, dbuser, dbtable, dbschema, dbpassword, localfile, resultdata):
	"""
	Import data into DB MySQL using system call.
	Now use the package MySQLdb because there are some problems with the system call: not all rows are written even if they exist. Running the script I can import 98 rows, from the shell the same call imports all 100 onces.
	"""
	# set up a tmp MySQL connection
	conn = MySQLdb.connect( 
		host = dbhost,
		port = int(dbport),
		user = dbuser,
		passwd = dbpassword,
		db = "information_schema",
		)
	#cursor = conn.cursor ()
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	
	# create database
	createdatabase = "CREATE DATABASE IF NOT EXISTS %s;" %(dbschema)
	##os.system("""mysql -u %s --password=%s -h %s --port %s -e "%s" """ %(dbuser, dbpassword, dbhost, dbport, createdatabase) )
	cursor.execute( createdatabase )

	# closing operations for tmp conn
	cursor.close()
	conn.close()
	
	# now that the schema exists, set up a MySQL connection to the target schema
	conn = MySQLdb.connect( 
		host = dbhost,
		port = int(dbport),
		user = dbuser,
		passwd = dbpassword,
		db = dbschema,
		)
	#cursor = conn.cursor ()
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	# create table
	## read order :: 'chr', 'start', 'end', 'strand', 'quality', 'NM', 'flag', 'cigar', 'MD', 'insert_size', 'AS', 'XS', 'SA'
	createtable = """
		CREATE TABLE IF NOT EXISTS %s (  
			prod_header varchar(255) NOT NULL, 
			r1_RG varchar(250) DEFAULT NULL,
			r1_chr  varchar(50) DEFAULT NULL,
			r1_start int(30) DEFAULT NULL,
			r1_end int(30) DEFAULT NULL,
			r1_strand varchar(50) DEFAULT NULL,
			r1_quality int(30) DEFAULT NULL,
			r1_NM int(30) DEFAULT NULL,
			r1_flag int(30) DEFAULT NULL,
			r1_cigar varchar(50) DEFAULT NULL,
			r1_cigar_components text,
			r1_offtarget varchar(50) DEFAULT NULL,
			r1_offtarget_type varchar(50) DEFAULT NULL,
			r1_offtarget_interval varchar(250) DEFAULT NULL,
			r1_offtarget_size varchar(250) DEFAULT NULL,
			r1_intarget_seq varchar(600) DEFAULT NULL,
			r1_MD varchar(250) DEFAULT NULL,
			r1_insert_size int(30) DEFAULT NULL,
			r1_AS int(30) DEFAULT NULL,
			r1_XS int(30) DEFAULT NULL,
			r1_SA varchar(255) DEFAULT NULL,
			r2_RG varchar(250) DEFAULT NULL,
			r2_chr varchar(50) DEFAULT NULL,
			r2_start int(30) DEFAULT NULL,
			r2_end int(30) DEFAULT NULL,
			r2_strand varchar(50) DEFAULT NULL,
			r2_quality int(30) DEFAULT NULL,
			r2_NM int(30) DEFAULT NULL,
			r2_flag int(30) DEFAULT NULL,
			r2_cigar varchar(50) DEFAULT NULL,
			r2_cigar_components text,
			r2_offtarget varchar(50) DEFAULT NULL,
			r2_offtarget_type varchar(50) DEFAULT NULL,
			r2_offtarget_interval varchar(250) DEFAULT NULL,
			r2_offtarget_size varchar(250) DEFAULT NULL,
			r2_intarget_seq varchar(600) DEFAULT NULL,
			r2_MD varchar(250) DEFAULT NULL,
			r2_insert_size int(30) DEFAULT NULL,
			r2_AS int(30) DEFAULT NULL,
			r2_XS int(30) DEFAULT NULL,
			r2_SA varchar(255) DEFAULT NULL,
			PRIMARY KEY (prod_header),
  			KEY r1_intarget_seq (r1_intarget_seq)
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
			INSERT INTO `%s`.`%s` (`prod_header`, `r1_RG`, `r1_chr`, `r1_start`, `r1_end`, `r1_strand`, `r1_quality`, `r1_NM`, `r1_flag`, `r1_cigar`, `r1_cigar_components`, `r1_offtarget`, `r1_offtarget_type`, `r1_offtarget_interval`, r1_offtarget_size, r1_intarget_seq, `r1_MD`, `r1_insert_size`, `r1_AS`, `r1_XS`, `r1_SA`, `r2_RG`, `r2_chr`, `r2_start`, `r2_end`, `r2_strand`, `r2_quality`, `r2_NM`, `r2_flag`, `r2_cigar`, `r2_cigar_components`, `r2_offtarget`, `r2_offtarget_type`, `r2_offtarget_interval`, r2_offtarget_size, r2_intarget_seq, `r2_MD`, `r2_insert_size`, `r2_AS`, `r2_XS`, `r2_SA`)
			VALUES 
			""" %(dbschema, dbtable)

		# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
		for row in resultdata[bedstart:bedend]:
			query_import += ' ("' + '", "'.join( str(x) for x in row ) + '"), '
		#print query_import.rstrip(', ')
		importdata = query_import.rstrip(', ').replace('"NULL"', "NULL")
		#print importdata
		cursor.execute( importdata )
		# increment bed indexes
		bedstart += step
		bedend += step

	cursor.close()
	conn.close()
	# return value
	return 1

def getDataFromDB(dbhost, dbport, dbuser, dbtable, dbschema, dbpassword, field_translation, ):
	"""
	reference_reads is the denominator count reference number. if it is 'ontarget' it means that reference reads count are the reads that are NOT target, thus the percentage will NOT represent the abundance but only a proportion with respect to ON vs OFF target occurrences. using 'all' the base reference reads will be the sum of all reads of the experiment, representing the overall count of reads (100%), thus the off target occurrences will be a subset of these reads.
	"""
	# set up a tmp MySQL connection
	conn = MySQLdb.connect( 
		host = dbhost,
		port = int(dbport),
		user = dbuser,
		passwd = dbpassword,
		db = dbschema,
		)
	
	# with this query you will get the reference number of reads (denominator of percentage to compute)
	reference_dictionary_all = {} # k = (chr, rg); v = reads count
	reference_dictionary_ontarget = {} # k = (chr, rg); v = reads count
	query_reference_ontarget = """
			SELECT r1_chr, r1_RG, count( * ) AS groupreads
			FROM `%s`.`%s`
			WHERE %s
				AND r1_offtarget LIKE 'False'
			GROUP BY r1_chr, r1_RG
			ORDER BY r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, 'groupreads'
		""" %(dbschema, dbtable, string_where_r1_NM)
	query_reference_all = """
			SELECT r1_chr, r1_RG, count( * ) AS groupreads
			FROM `%s`.`%s`
			WHERE %s
			GROUP BY r1_chr, r1_RG
			ORDER BY r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, 'groupreads'
		""" %(dbschema, dbtable, string_where_r1_NM)
	# get the reference counts
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	#print query_reference_all
	cursor.execute(query_reference_all)
	results = cursor.fetchall()
	# now fill the reference dictionary
	for result in results:
		reference_dictionary_all[(result["r1_chr"], result["r1_RG"])] = float(result['groupreads'])
	cursor.close()
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	# get the reference counts
	cursor.execute(query_reference_ontarget)
	results = cursor.fetchall()
	# now fill the reference dictionary
	for result in results:
		reference_dictionary_ontarget[(result["r1_chr"], result["r1_RG"])] = float(result['groupreads'])
	cursor.close()
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)

	# now query target data. r1_NM is max allowed mismatches, set to 10, is a contraint for read quality in terms of mapping and base quality. r1_offtarget_type null means to count only off-target sites and not control reads (without any variation off-target). LENGTH(r1_offtarget_type) <= 1 because here we avoid potential double variation classification (D|I, ecc).
	query_target = """
			SELECT r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, r1_offtarget_size, count( * ) AS groupreads
			FROM `%s`.`%s`
			WHERE %s
				AND r1_offtarget NOT LIKE 'False'
				AND LENGTH(r1_offtarget_type) <= 1
			GROUP BY r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval
			ORDER BY r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, 'groupreads'
			""" %(dbschema, dbtable, string_where_r1_NM)
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	cursor.execute(query_target)
	results = cursor.fetchall()
	offtarget_dictionary = {} # k = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval); v = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, r1_offtarget_size, groupreads) + (abundance percentage, ratio on vs off)
	for result in results:
		ratio_on_off = result["groupreads"]/float(reference_dictionary_ontarget[(result["r1_chr"], result["r1_RG"])])
		abundance_perc_offtarget = 100*(result["groupreads"]/float(reference_dictionary_all[(result["r1_chr"], result["r1_RG"])]))
		res_data = ()
		for col in ['r1_chr', 'r1_RG', 'r1_offtarget_type', 'r1_offtarget_interval', 'r1_offtarget_size', 'groupreads']:
			res_data += (result[col], )
		offtarget_dictionary[(result["r1_chr"], result["r1_RG"], result["r1_offtarget_type"], result["r1_offtarget_interval"])] = res_data + (ratio_on_off, abundance_perc_offtarget, reference_dictionary_ontarget[(result["r1_chr"], result["r1_RG"])], reference_dictionary_all[(result["r1_chr"], result["r1_RG"])] )
	cursor.close()
	conn.close()

	outdata = [] # list of data as array of array, to be easily written
	for k, v in offtarget_dictionary.iteritems():
		outdata.append(v)
	outdatasorted = sorted(outdata, key=itemgetter(0, 1, 2, 3), reverse = False)
	outdatalist = [list(x) for x in outdatasorted] # list of lists

	return outdatalist, offtarget_dictionary

def checkSampleGroupID(controlSample, dbhost, dbport, dbuser, dbtable, dbschema, dbpassword, field_translation, ):
	"""
	From the DB check if the input RG group exists.
	Output: true or false
	"""
	controlSample_inDB = False
	# set up a tmp MySQL connection
	conn = MySQLdb.connect( 
		host = dbhost,
		port = int(dbport),
		user = dbuser,
		passwd = dbpassword,
		db = dbschema,
		)
	query = """
		SELECT DISTINCT r1_RG
		FROM `%s`.`%s`
		WHERE 1
		"""
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	cursor.execute(query)
	results = cursor.fetchall()
	for result in results:
		if controlSample == result['r1_RG']:
			controlSample_inDB = True
	cursor.close()
	conn.close()
	return controlSample_inDB

def acquireMinimumValueFromControl(controlSample, dbhost, dbport, dbuser, dbtable, dbschema, dbpassword, field_translation, ):
	"""
	From the DB get the minimum value from the control sample. 
	This function returns a dictionary:: k = (region, sample), v = min value.
	"""
	sample_minval_dictionary = {} # k = (region, sample), v = (region, sample, all reads, off target reads, percentage)
	# set up a tmp MySQL connection
	conn = MySQLdb.connect( 
		host = dbhost,
		port = int(dbport),
		user = dbuser,
		passwd = dbpassword,
		db = dbschema,
		)
	
	query_reference_control = """
			SELECT r1_chr, r1_RG, count( * ) AS groupreads
			FROM `%s`.`%s`
			WHERE %s
				AND r1_RG = '%s'
			GROUP BY r1_chr, r1_RG
			ORDER BY r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, 'groupreads'
		""" %(dbschema, dbtable, string_where_r1_NM, controlSample)
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	cursor.execute(query_reference_control)
	results = cursor.fetchall()
	for result in results:
		sample_minval_dictionary[(result["r1_chr"], result["r1_RG"])] = (result["r1_chr"], result["r1_RG"], result["groupreads"])
	cursor.close()

	query_target_control = """
			SELECT r1_chr, r1_RG, count( * ) AS groupreads
			FROM `%s`.`%s`
			WHERE %s
				AND r1_offtarget NOT LIKE 'False'
				AND LENGTH(r1_offtarget_type) <= 1
				AND r1_RG = '%s'
			GROUP BY r1_chr, r1_RG
			ORDER BY r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, 'groupreads'
		""" %(dbschema, dbtable, string_where_r1_NM, controlSample)
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	cursor.execute(query_target_control)
	results = cursor.fetchall()
	for result in results:
		sample_minval_dictionary[(result["r1_chr"], result["r1_RG"])] += (result["groupreads"],)
	cursor.close()

	# now re-parse the dictionary to add 0 values for not found off-target reads (in the previous query not all samples return a value greater than 0)
	for k, v in sample_minval_dictionary.iteritems():
		if len(v) < 4:
			v += (0,)

	conn.close()
	return sample_minval_dictionary

def getHeaderReadsFromOffTargets(dbhost, dbport, dbuser, dbtable, dbschema, dbpassword, ):
	"""
	output: dictionary
	"""
	dictionary = {} # k = (region, sample), v = header
	read_sample_dictionary = {} # k = header, v = 
	# set up a tmp MySQL connection
	conn = MySQLdb.connect( 
		host = dbhost,
		port = int(dbport),
		user = dbuser,
		passwd = dbpassword,
		db = dbschema,
		)
	
	query_reference_control = """
			SELECT prod_header, r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, r1_offtarget_size, count( * ) AS groupreads
			FROM `%s`.`%s`
			WHERE %s
				AND r1_offtarget NOT LIKE 'False'
				AND LENGTH(r1_offtarget_type) <= 1
			GROUP BY r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval
			ORDER BY r1_NM asc, r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval, 'groupreads'
		""" %(dbschema, dbtable, string_where_r1_NM)
	cursor = conn.cursor (MySQLdb.cursors.DictCursor)
	cursor.execute(query_reference_control)
	results = cursor.fetchall()
	for result in results:
		# data for region-sample dictionary
		if dictionary.has_key((result["r1_chr"], result["r1_RG"])):
			dictionary[(result["r1_chr"], result["r1_RG"])].append(result["prod_header"])
		else:
			dictionary[(result["r1_chr"], result["r1_RG"])] = [ result["prod_header"] ]
		# data for reads dictionary
		res_data = ()
		for col in ['r1_chr', 'r1_RG', 'r1_offtarget_type', 'r1_offtarget_interval']:
			res_data += (result[col], )
		read_sample_dictionary[result["prod_header"]] = (result["r1_chr"], result["r1_RG"], result["r1_offtarget_type"], result["r1_offtarget_interval"])
	cursor.close()

	conn.close()
	return dictionary, read_sample_dictionary

def getTargetSequences(targetFile, ):
	"""
	"""
	target_sequence = {}
	try:
		reader = csv.reader(open(targetFile, 'rb'), delimiter = "\t")
		for row in reader:
			if len(row)>1:
				target_sequence[row[0]] = row[1]
	except csv.Error, e:
		sys.exit("[AP]\tError while writing CSV file %s. Check paths and names." % (targetFile))
	return target_sequence

def findReadSequencesFromBAM_CLI(readsfile, bamfile, readslist):
	"""
	Command line run pf the filtering process. This option will create the BAM tmp file.
	"""
	readseqs_dictionary = {}
	tmpbamfile = bamfile + ".tmp"
	tmpcsvfile = readsfile + ".csv"
	# PICARD
	print "[AP]\tCreating BAM from reads list."
	#os.system("cat %s" %(readsfile))
	picard_cmd = "FilterSamReads INPUT=\"%s\" FILTER=includeReadList RLF=\"%s\" SO=coordinate O=\"%s\"" %(bamfile, readsfile, tmpbamfile)
	print picard_cmd
	os.system(picard_cmd)
	# SAMTOOLS
	print "[AP]\tExtracting sequences and creating CSV."
	os.system("samtools view %s | cut -f 1,10 > %s" %(tmpbamfile, tmpcsvfile))
	# now parse the output to get sequence dictionary
	try:
		reader = csv.reader(open(tmpcsvfile, 'rb'), delimiter = "\t")
		for row in reader:
			readseqs_dictionary[row[0]] = row[1]
	except csv.Error, e:
		sys.exit("[AP]\tError while writing CSV file %s. Check paths and names." % (targetFile))
	return readseqs_dictionary # k = read header , v = sequence

def findReadSequencesFromBAM(readsfile, bamfile, readslist):
	"""
	"""
	readseqs_dictionary = {}
	readaln_dictionary = {} # k = header, v = bam pysam object
	bamtmp = pysam.Samfile(bamfile) # the bam file to process
	for aln in bamtmp:
		if aln.qname in readslist:
			readseqs_dictionary[aln.qname] = aln.seq
			readaln_dictionary[aln.qname] = aln
	bamtmp.close()
	return readseqs_dictionary, readaln_dictionary # k = read header , v = sequence

def createFastaFiles_runMSA(fastaFolder, targetseqs_dictionary, readseqs_dictionary, groupheader_dictionary, regions_dictionary, offtarget_dictionary, header_samplerg_dictionary, split_files_by_offtarget_type = True):
	"""
	create the FASTA files from dictionaries: 1 file for sample/group with first row as the reference region
	## createFastaFiles(args.fastaFolder, targetseqs_dictionary, readseqs_dictionary, groupheader_dictionary) # 
	"""
	if split_files_by_offtarget_type: # if it is required to split fasta files in D and I (off target types), then create 2 files each sample-region, else only 1
		for (region, sample), readheader_list in groupheader_dictionary.iteritems():
			print "[AP]\tCreating Fasta:", region, sample
			# init fasta file
			filefasta_I = os.path.join(fastaFolder, region + "." + sample + ".I.fa")
			ffo_I = open(filefasta_I, 'w')
			ffo_I.write(">%s\n%s\n" %(region, targetseqs_dictionary[region]) )
			
			filefasta_D = os.path.join(fastaFolder, region + "." + sample + ".D.fa")
			ffo_D = open(filefasta_D, 'w')
			ffo_D.write(">%s\n%s\n" %(region, targetseqs_dictionary[region]) )
			
			for header in readheader_list:
				if header_samplerg_dictionary[header][2] == 'D':
					# ffo_D.write(">%s|%s|%s|%.4f|%s\n%s\n" %(sample, header_samplerg_dictionary[header][2], offtarget_dictionary[header_samplerg_dictionary[header]][5], offtarget_dictionary[header_samplerg_dictionary[header]][7], header, readseqs_dictionary[header]) ) # RG,D/I,interval,#reads,abundance% | header_samplerg_dictionary:: k = header, v = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval)
					ffo_D.write(">%s|%s|%s|%s|%.4f\n%s\n" %(region, sample, header_samplerg_dictionary[header][2], header_samplerg_dictionary[header][3], offtarget_dictionary[header_samplerg_dictionary[header]][7], readseqs_dictionary[header]) ) # RG,D/I,interval,#reads,abundance% | header_samplerg_dictionary:: k = header, v = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval)
				if header_samplerg_dictionary[header][2] == 'I':
					ffo_I.write(">%s|%s|%s|%s|%.4f\n%s\n" %(region, sample, header_samplerg_dictionary[header][2], header_samplerg_dictionary[header][3], offtarget_dictionary[header_samplerg_dictionary[header]][7], readseqs_dictionary[header]) ) # RG,D/I,interval,#reads,abundance% | header_samplerg_dictionary:: k = header, v = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval)
			
			ffo_I.close()
			ffo_D.close()

			# run MSA
			print "[AP]\tRunning MSA (clustalw2):", region, sample
			filefasta_postmsa_I = os.path.join(fastaFolder, region + "." + sample + ".I.msa.fa" )
			filefasta_postmsa_D = os.path.join(fastaFolder, region + "." + sample + ".D.msa.fa" )
			# CLUSTALW2
			clustal_cmd_I = "clustalw2 -ALIGN -TYPE=DNA -OUTFILE=%s -OUTPUT=FASTA -OUTORDER=INPUT -INFILE=%s" %(filefasta_postmsa_I, filefasta_I, )
			clustal_cmd_D = "clustalw2 -ALIGN -TYPE=DNA -OUTFILE=%s -OUTPUT=FASTA -OUTORDER=INPUT -INFILE=%s" %(filefasta_postmsa_D, filefasta_D, )
			print clustal_cmd_I, clustal_cmd_D
			os.system(clustal_cmd_I)
			os.system(clustal_cmd_D)
			# MVIEW
			print "[AP]\tCreating HTML files, Running MVIEW:", region, sample
			# no colors in output
			outhtml_I = os.path.join(fastaFolder, region + "." + sample + ".I.msa.html" )
			outhtml_D = os.path.join(fastaFolder, region + "." + sample + ".D.msa.html" )
			os.system("mview -in pearson -ruler on -html head -css on -range %d:%d %s > %s " %(int(regions_dictionary[region][0])-10, int(regions_dictionary[region][1])+10, filefasta_postmsa_I, outhtml_I) )
			os.system("mview -in pearson -ruler on -html head -css on -range %d:%d %s > %s " %(int(regions_dictionary[region][0])-10, int(regions_dictionary[region][1])+10, filefasta_postmsa_D, outhtml_D) )
			# adding colors
			outhtml_I_col = os.path.join(fastaFolder, region + "." + sample + ".I.msa.colors.html" )
			outhtml_D_col = os.path.join(fastaFolder, region + "." + sample + ".D.msa.colors.html" )
			os.system("mview -in pearson -ruler on -html head -css on -coloring any -colormap CLUSTAL -range %d:%d %s > %s " %(int(regions_dictionary[region][0])-10, int(regions_dictionary[region][1])+10, filefasta_postmsa_I, outhtml_I_col) )
			os.system("mview -in pearson -ruler on -html head -css on -coloring any -colormap CLUSTAL -range %d:%d %s > %s " %(int(regions_dictionary[region][0])-10, int(regions_dictionary[region][1])+10, filefasta_postmsa_D, outhtml_D_col) )
	else: # case of not split_files_by_offtarget_type (split_files_by_offtarget_type=False)
		for (region, sample), readheader_list in groupheader_dictionary.iteritems():
			print "[AP]\tCreating Fasta:", region, sample
			# init fasta file
			filefasta = os.path.join(fastaFolder, region + "." + sample + ".fa")
			ffo = open(filefasta, 'w')
			ffo.write(">%s\n%s\n" %(region, targetseqs_dictionary[region]) )
			for header in readheader_list:
				ffo.write(">%s|%s|%s|%s|%.4f\n%s\n" %(region, sample, header_samplerg_dictionary[header][2], header_samplerg_dictionary[header][3], offtarget_dictionary[header_samplerg_dictionary[header]][7], readseqs_dictionary[header]) ) # RG,D/I,interval,#reads,abundance% | header_samplerg_dictionary:: k = header, v = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval)
			ffo.close()
			# run MSA
			print "[AP]\tRunning MSA (clustalw2):", region, sample
			filefasta_postmsa = os.path.join(fastaFolder, region + "." + sample + ".msa.fa" )
			# CLUSTALW2
			clustal_cmd = "clustalw2 -ALIGN -TYPE=DNA -OUTFILE=%s -OUTPUT=FASTA -OUTORDER=INPUT -INFILE=%s" %(filefasta_postmsa, filefasta, )
			print clustal_cmd
			os.system(clustal_cmd)
			# MVIEW
			print "[AP]\tCreating HTML files, Running MVIEW:", region, sample
			outhtml = os.path.join(fastaFolder, region + "." + sample + ".msa.html" )
			os.system("mview -in pearson -ruler on -html head -css on -coloring any -colormap CLUSTAL -range %d:%d %s > %s " %(int(regions_dictionary[region][0])-10, int(regions_dictionary[region][1])+10, filefasta_postmsa, outhtml) )
	return True

def createBAMFiles(fastaFolder, targetseqs_dictionary, readseqs_dictionary, groupheader_dictionary, regions_dictionary, offtarget_dictionary, header_samplerg_dictionary, readbam_dictionary, input_bamfile):
	"""
	
	"""
	##### versione corretta se funzionasse l'header info...
	# inbam = pysam.Samfile(input_bamfile) # the bam file to process
	# inbam_header = inbam.header # get the header from inbam
	# inbam.close()
	# for (region, sample), readheader_list in groupheader_dictionary.iteritems():
	# 	outbam_filename = os.path.join(fastaFolder, region + "." + sample + ".bam")
	# 	print "[AP]\tCreating BAM file:", region, sample, outbam_filename
	# 	outbam = pysam.Samfile( outbam_filename, "wh", header = inbam_header )
	# 	for header in readheader_list:
	# 		outbam.write(readbam_dictionary[header])
	# 	outbam.close()
	# 	#print "[AP]\t\tindexing this BAM file:", region, sample
	# 	#pysam.index(outbam_filename)
	for (region, sample), readheader_list in groupheader_dictionary.iteritems():
		print "[AP]\tCreating BAM files:", region, sample
		filereads = os.path.join(fastaFolder, region + "." + sample + ".readlist")
		ffo = open(filereads, 'w')
		for header in readheader_list:
			ffo.write(header+"\n")
		ffo.close()
		outbam_filename = os.path.join(fastaFolder, region + "." + sample + ".bam")
		# PICARD
		print "[AP]\tCreating BAM from reads list."
		#os.system("cat %s" %(readsfile))
		picard_cmd = "FilterSamReads INPUT=\"%s\" FILTER=includeReadList RLF=\"%s\" SO=coordinate O=\"%s\"" %(input_bamfile, filereads, outbam_filename)
		print picard_cmd
		os.system(picard_cmd)
		os.system("samtools index %s" %(outbam_filename))

	return True

def createBEDFiles(fastaFolder, targetseqs_dictionary, readseqs_dictionary, groupheader_dictionary, regions_dictionary, offtarget_dictionary, header_samplerg_dictionary):
	"""
	"""
	for (region, sample), readheader_list in groupheader_dictionary.iteritems():
		print "[AP]\tCreating BED file:", region, sample
		# init fasta file
		filebed = os.path.join(fastaFolder, region + "." + sample + ".bed")
		ffo = open(filebed, 'w')
		for header in readheader_list:
			# draw a different range for D or I (with I DO NOT remove a bp from start)
			if header_samplerg_dictionary[header][2] == 'D':
				ffo.write("%s\t%d\t%d\t%s|%s|%.4f|%d\n" %(	region, 
													int(header_samplerg_dictionary[header][3].split('-')[0])-1, 
													int(header_samplerg_dictionary[header][3].split('-')[1]), 
													sample, 
													header_samplerg_dictionary[header][2], 
													offtarget_dictionary[header_samplerg_dictionary[header]][7], 
													offtarget_dictionary[header_samplerg_dictionary[header]][5], ) ) # RG,D/I,interval,#reads,abundance% | header_samplerg_dictionary:: k = header, v = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval) 
			else:
				ffo.write("%s\t%d\t%d\t%s|%s|%.4f|%d\n" %(	region, 
													int(header_samplerg_dictionary[header][3].split('-')[0]), 
													int(header_samplerg_dictionary[header][3].split('-')[1]), 
													sample, 
													header_samplerg_dictionary[header][2], 
													offtarget_dictionary[header_samplerg_dictionary[header]][7], 
													offtarget_dictionary[header_samplerg_dictionary[header]][5], ) ) # RG,D/I,interval,#reads,abundance% | header_samplerg_dictionary:: k = header, v = (r1_chr, r1_RG, r1_offtarget_type, r1_offtarget_interval) 
		ffo.close()
	return True

	

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
	print "[AP]\tIndexing BAM (this is required if you did not do it)."
	indexing = pysam.index(args.bamfile) # indexing is BLANK
	print "[AP]\tAcquiring data BAM."
	bam = pysam.Samfile(args.bamfile) # the bam file to process
	print "[AP]\tAcquiring region list to use while slicing BAM."
	regions, regions_dictionary = parseTargetFile(args.regions, int(args.bpoffset))
	print "[AP]\tNow looping over BAM alignments in a SINGLE process..."
	resultdata = singleWorker(regions, bam, out_q=None)

	# print "[AP]\tNow looping over BAM alignments Splitting run in processes..."
	# out_q = Queue() # init process queue
	# nprocs = int(args.processes)
	# chunksize = int(math.ceil(len(regions) / float(nprocs))) # find the chunk wrt processes for running regions
	# procs = [] # array of processes
	# # run processes in a range of chunks
	# for i in range(nprocs):
	# 	p = multiprocessing.Process(target=worker, args=(regions[chunksize * i:chunksize * (i + 1)], bam, out_q, )) 
	# 	procs.append(p)
	# 	p.start()
	# # Collect all results into a single result list.
	# print "[AP]\tGetting results from different processes."
	# resultdata = [] # this is what will be written into the file
	# for i in range(nprocs):
	# 	resultdata += out_q.get()
	# # Wait for all worker processes to finish
	# print "[AP]\tWaiting that all processes end."
	# for p in procs:
	# 	p.join()
	
	# # write output results from joined array of strings
	header = ['header'] + ['r1_'+str(x) for x in readorder] + ['r2_'+str(x) for x in readorder]
	tmpfile = args.dbtable + ".csv"
	print "[AP]\tWriting results into output file", tmpfile
	try:
		writer = csv.writer(open(tmpfile, 'wb'), delimiter = "\t")
		writer.writerow(header)
		for row in resultdata:
			writer.writerow(row)
	except csv.Error, e:
		sys.exit("[AP]\tError while writing CSV BED file %s. Check paths and names." % (tmpfile))

	# import data into DB
	print "[AP]\tImporting data into DB [dbhost <%s>, dbport <%s>, dbuser <%s>, dbschema <%s>, dbtable <%s>]\n\t Number of reads to be imported: %d" %(args.dbhost, args.dbport, args.dbuser, args.dbschema, args.dbtable, len(resultdata))
	dbimportexit = importData(args.dbhost, args.dbport, args.dbuser, args.dbtable, args.dbschema, args.dbpassword, tmpfile, resultdata)

	# create the output percentage file of the off-targets
	##sleep(3)
	sampleGroup_valid = False
	if args.controlSample not in ['NOT SET']:
		print "[AP]\tYou set the control sample (background noise base):", args.controlSample, ", now on check in the DB."
		sampleGroup_valid = checkSampleGroupID(args.controlSample, args.dbhost, args.dbport, args.dbuser, args.dbtable, args.dbschema, args.dbpassword, )

	offtargetlist, offtargetdict = getDataFromDB(args.dbhost, args.dbport, args.dbuser, args.dbtable, args.dbschema, args.dbpassword, field_translation, ) # list of lists
	# # write output results from joined array of strings
	header = ['r1_chr', 'r1_RG', 'r1_offtarget_type', 'r1_offtarget_interval', 'r1_offtarget_size', 'groupreads', 'ratio', 'percentage', 'ref_ratio', 'ref_all']
	# print "[AP]\tWriting results into output file", args.outfilename
	# try:
	# 	writer = csv.writer(open(args.outfilename, 'wb'), delimiter = "\t")
	# 	writer.writerow([field_translation[x] for x in header])
	# 	for row in offtargetlist:
	# 		writer.writerow(row)
	# except csv.Error, e:
	# 	sys.exit("[AP]\tError while writing CSV file %s. Check paths and names." % (args.outfilename))
	
	print "[AP]\tWriting results into output file", args.outfilename
	fo = open(args.outfilename, 'w')
	fo.write('\t'.join( str(field_translation[x]) for x in header) + "\n")
	for row in offtargetlist:
		fo.write('\t'.join(str(x) for x in row) + "\n")
	fo.close()
	
	# section for MSA and reads extraction/formatting (fasta)
	if args.doClustering:
		# 1. get headers of sample reads from DB groups
		groupheader_dictionary, header_samplerg_dictionary = getHeaderReadsFromOffTargets(args.dbhost, args.dbport, args.dbuser, args.dbtable, args.dbschema, args.dbpassword, ) # k = groups (tuple), v = list of reads headers
		#print groupheader_dictionary
		# 2. get the reference sequences
		targetseqs_dictionary = getTargetSequences(args.targetSequenceFile, ) # k = target id, v = sequence
		# 3. get reads sequences from input BAM file into a dictionary
		reads_headers_list = groupheader_dictionary.values() # list of reads headers
		reads_headers_simple_list = []
		# tmp file of reads 
		tmpfile = os.path.join( os.getcwd(), strftime("%Y%m%d%H%M%S") + ".reads.list")
		print "[AP]\tWriting tmp file of reads (sample reads from BAM)", tmpfile
		tmpfo = open(tmpfile, 'w')
		for row in reads_headers_list:
			for h in row:
				tmpfo.write(h+"\n")
				reads_headers_simple_list.append(h)
		tmpfo.close()
		readseqs_dictionary, readbam_dictionary = findReadSequencesFromBAM(tmpfile, args.bamfile, reads_headers_simple_list) # k = read header, v = sequence
		##print readseqs_dictionary
		# 4. create the FASTA files from dictionaries: 1 file for sample/group with first row as the reference region
		createFastaFiles_runMSA(args.fastaFolder, targetseqs_dictionary, readseqs_dictionary, groupheader_dictionary, regions_dictionary, offtargetdict, header_samplerg_dictionary) # 
		# create bam of reads
		createBAMFiles(args.fastaFolder, targetseqs_dictionary, readseqs_dictionary, groupheader_dictionary, regions_dictionary, offtargetdict, header_samplerg_dictionary, readbam_dictionary, args.bamfile) # bam files with indexes
		createBEDFiles(args.fastaFolder, targetseqs_dictionary, readseqs_dictionary, groupheader_dictionary, regions_dictionary, offtargetdict, header_samplerg_dictionary)

	elapsed_time = time.time() - start_time # get read of finishing time
	print "\n[AP]\tTask Finished, closing.\n\tElapsed time: %.2f [seconds]\n" %(elapsed_time)
	

# sentinel
if __name__ == "__main__":
    main()
