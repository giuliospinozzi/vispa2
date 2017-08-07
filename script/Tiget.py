#!/usr/bin/python
# -*- coding: utf-8 -*-
import MySQLdb
from time import gmtime, strftime
import sys, os, re, csv, email, math
import subprocess
from string import Template
from operator import itemgetter, attrgetter
from Bio import Entrez, SeqIO
import random
from threading import Thread
import threading, multiprocessing, subprocess, signal
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

class Linkage:
	"""
	Linkage Analysis class.
	
	Notes:
	- this class is useful to format files for Merlin from Illumina Genome Studio exported data using Merlin plugin (that has been designed for associacion studies)
	- the main data structure is a dictionary, designed in this way (so far, this works fine for a single family, but is also extended to multiple families only if subjects among families have different IDs:
		subj_geno = {} -> this is the reference structure! composition:
			k = subject id
			v = dictionary:
				k = geno: genotype data, all markers, all chr -> v = list of markers
					fam: familiar data -> v = list of data
					chr 1..N: slice of genotypes for the chr key -> v = list of markers
					genopos: relative positions of slices for each chr -> v = dictionary:
						k = chr -> v = list of positions relative to geno markers
					markerpos: list of positions of only filtered markers for this chr
						k = chr -> v = list of positions of markers refered to the overall map file previously acquired as tuple array
					markerlist: list of ids
						k = chr -> v = tuple array
	- the idea is to acquire data from PED/MAP/DAT of genotyped subjects and build the basic data structure, and also add the "non-genotyped" subjects to the main data structure, so that I can write a real PED with also 0s/missing data
	- since just formatting PED, MAP and DAT of all markers does not work, I created methods to overcome different problems:
		1. each chromosome is better to be separated from the others
		2. only 1 marker for each position in cM -> if I have 2 markers in the same position, make a filter
		3. due to FATAL MEMORY FAULT, I also added parameters to filter markers by at least a static step (default is 0.2 cM by distance)
	
	Example of usage:
		#########################################################################################
		## input data
		pedg = "/home/andrea/Dropbox/HSR/Martinelli/MS/Linkage/input/MerlinMS/MerlinGO/InfiniumPP_NEUFAM05_23-4-12.ped"
		pedf = "/home/andrea/Dropbox/HSR/Martinelli/MS/Linkage/input/ms.pedonly.linux.ped"
		mapf = "/home/andrea/Dropbox/HSR/Martinelli/MS/Linkage/input/MerlinMS/MerlinGO/InfiniumPP_NEUFAM05_23-4-12.map"
		path = "/Users/andrea/Dropbox/HSR/Martinelli/MS/Linkage/input"
		pednew = os.path.join(path, "ms.full.subf.ped")
		prefix = "ms5.full.chr"
		#########################################################################################
		## first acquire data
		dsg = {}
		acquirePed(pedg, dsg) # -> from genotypes
		acquirePed(pedf, dsg) # -> from ped only file
		## now sort map by chr
		snpl = sortByChr(mapf)
		
		## extract/filter chr that you want, i.e. 21
		### -> alternative n1
		# one by one
		extractChr(snpl, dsg, 20) # -> exclude NaN markers as much as you can!
		## write data: PED with only chr 21, MAP/DAT with sorted/filtered markers
		writePed(ped20, dsg, what = 20) # write only 21
		writeMapDat(map20, dat20, dsg, 20, diseasecol = 0 )
		## now is filter distance default 0.2 ;)
		extractChr(snpl, dsg, 21, filter = True) # -> exclude NaN markers as much as you can!
		
		### -> alternative n2
		# bulk writing
		## extract all available chr from data
		available_chr = set()
		for m in snpl:
			available_chr.add( m[0] )
		## without any filter
		for k in available_chr:
			extractChr(snpl, dsg, k) # -> exclude NaN markers as much as you can!
			writePed(os.path.join(path, prefix + str(k) +".ped"), dsg, what = k) # write only 21
			writeMapDat(os.path.join(path, prefix + str(k) +".map"), os.path.join(path, prefix + str(k) +".dat"), dsg, k, diseasecol = 0 )
		## with filter ON
		for k in available_chr:
			extractChr(snpl, dsg, k, filter = True) # -> exclude NaN markers as much as you can!
			writePed(os.path.join(path, prefix + str(k) +".ped"), dsg, what = k) # write only 21
			writeMapDat(os.path.join(path, prefix + str(k) +".map"), os.path.join(path, prefix + str(k) +".dat"), dsg, k, diseasecol = 0 )
		#########################################################################################
		## then, add annotation to tabular files
		dbschema = 'ucsc_hg19'
		conn = MySQLdb.connect(
				  host = 'localhost',
				  user = 'XXXXXX',
				  passwd = 'XXXXXXXX',
				  db = dbschema,
				  )
		cursor = conn.cursor (MySQLdb.cursors.DictCursor)

		infile_tbl = "/home/andrea/Dropbox/HSR/Martinelli/MS/Linkage/output/20120503/data/out.ms5.full.d02.chr8.nplexp.r15-nonparametric.tbl"
		outfile_tbl = "/home/andrea/Dropbox/HSR/Martinelli/MS/Linkage/output/20120503/data/out.ms5.full.d02.chr8.nplexp.r15-nonparametric.annotated.tbl"
		annotateTbl(cursor, infile_tbl, outfile_tbl)
		cursor.close()
		conn.close()
		#########################################################################################
	"""
	def __init__(self, ):
		"""
		Init class
		"""
		self.pedfile = None ### -> output of Illumina Genome Studio
		self.mapfile = None ### -> output of Illumina Genome Studio
		self.datfile = None ### -> output of Illumina Genome Studio

	def annotateTbl(cursor, tblfile, tblfile_annotated, ):
		"""
		From file tbl to tblfile_annotated file (output) using DB, UCSC HG19
		"""
		writer = csv.writer(open(tblfile_annotated, 'wb'), delimiter = "\t")
		with open(tblfile, 'rb') as inf:
			reader = csv.reader(inf, delimiter = "\t")
			for row in reader:
				if row[0] not in ['na', 'CHR'] :
					marker = row[2].strip()
					print "Annotating chr %s marker %s" %(row[0], marker)
					marker_data = annotateSnp(cursor, marker)
					annotation = []
					for n in ['chrom', 'chromStart', 'chromEnd', 'strand', 'refNCBI', 'refUCSC', 'observed', 'genesymbol']:
						annotation.append(marker_data[n])
					writer.writerow( row + annotation )
					#if not marker_result.has_key(maker):
					#	marker_result[marker] = {'results': (row) }
					#else:
					#	print "\nATTENTION: marker %s already exist!!\n" %(marker)
				else:
					annotation = ['chrom', 'chromStart', 'chromEnd', 'strand', 'refNCBI', 'refUCSC', 'observed', 'genesymbol']
					writer.writerow( row + annotation )
	
	def annotateSnp(cursor, rs, ):
		"""
		from ucsc hg19 returns dictionary of snp annotation
		"""
		query = "select chrom, chromStart, chromEnd, strand, refNCBI, refUCSC, observed from ucsc_hg19.snp132 where name = '%s' limit 1" %(rs)
		cursor.execute(query)
		results = cursor.fetchall()
		snp_data = None
		for r in results:
			snp_data = {'chrom': r['chrom'], 'chromStart': r['chromStart'], 'chromEnd': r['chromEnd'], 'strand': r['strand'], 'refNCBI': r['refNCBI'], 'refUCSC': r['refUCSC'], 'observed': r['observed'], 'genesymbol': ''}
			geneq = "select name2 from ucsc_hg19.refGene where txStart < '%d' and txEnd > '%d' and chrom = '%s' limit 1" %(snp_data['chromStart'], snp_data['chromEnd'], snp_data['chrom'])
			# print geneq
			cursor.execute(geneq)
			generes = cursor.fetchall()
			for g in generes:
				snp_data['genesymbol'] = g['name2']
		return snp_data
	
	def sortByChr(map, ):
		"""
		Important: static structure of tuple array defined as:
			chr, markerId, cM, markerindex, genoindex
		"""
		map_data_tmp = [] # tuple of values
		with open(map, 'rb') as inf:
			reader = csv.reader(inf, delimiter = "\t")
			position_index = 0
			marker_index = 0
			for row in reader:
				if len(row) == 3:
					if row[2] != 'NaN': ### this IF is important to avoid invalid sorting due to NAN existance
						elem = row[0], row[1], float(row[2]), marker_index, position_index
					position_index += 1
					marker_index += 1
					map_data_tmp.append( elem )
		map_data_sortedbychr = sorted(map_data_tmp, key=itemgetter(0, 2), reverse = False)
		return map_data_sortedbychr
	
	def extractChr(maplist, subj_geno, chr, filter = False, filter_types = {'distance': 0.2, }, ):
		"""
		Main points:
			- given the sorted list of markers, filter target chr
			- extract only markers that follow these rules:
				* no NaN -> this rule is written in the method sortByChr, since affects sorting (python errors on float(nan))
				* no overlapping markers on the same position -> this rule is coded here, since is a direct filter to marker, not only a sorting
		"""
		# chr = id of target chr -> will be included as separated key into subj_geno
		## OBSOLETE: maplist_sorted = sorted(maplist, key=itemgetter(0, 2), reverse = False) # just to be safe, re-sort markers by chr and CM position
		for sbjid, datadict in subj_geno.iteritems():
			if datadict.has_key('geno'):
				print sbjid, "-> Acquiring genotypes by filtering only chromosome", chr
				# pair data! -> get genotype positions first, following rules
				markercm = 0.0 # this marker centimorgan: this is ther reference starting value of cM
				genopos = [] # list of positions to filter from ped -> this index refers to ped columns of only genotypes (no first subject columns, like fam id, sbj id, father, mother, sex, status)
				genopos_exact = [] # list of positions relative to ped columns: is the same of genopos but takes into account the case of doubling values if paired markers are not separated by space but by \t
				markerpos = [] # positions of markers relative to map file -> tuple array (NB: do not sort it anymore!)
				markerlist = [] # tuple array of marker ids, subset of overall marker tuple array
				# preceding_makercm = 0.0 # cm used for filters
				##### markerindex = 0 # index of marker position
				for snp in maplist:
					if str(chr) == str(snp[0]) and float(snp[2]) > markercm and snp[2] != 'NaN':
						if filter: ## if filter option is active, default is not!
							## parse filters and identify action: here only one active per time
							action = None
							for k, v in filter_types.iteritems():
								if v is not None:
									action = k
							step = filter_types[action] ## so far only distance available, thus step is distance between consecutive markers
							## filter markers by this filter
							if float(snp[2]) >= markercm + step:
								genopos.append(int(snp[-1])) # append this index position (relative to file), to geno list
								markercm = float(snp[2]) # to avoid overlapping cm, update reference cm to this marker position
								markerpos.append(int(snp[-2])) ## this position is static!
								markerlist.append(snp)
						else:
							genopos.append(int(snp[-1])) # append this index position (relative to file), to geno list
							markercm = float(snp[2]) # to avoid overlapping cm, update reference cm to this marker position
							markerpos.append(int(snp[-2])) ## this position is static!
							markerlist.append(snp)
						#### markerindex += 1
				## case 1: len genotypes = len marker list
				if len(datadict['geno']) == len(maplist):
					# now parse gnotypes and cut ped -> create key of chr with snp array
					for snp in genopos:
						if subj_geno[sbjid].has_key(chr):
							subj_geno[sbjid][chr].append(datadict['geno'][snp])
							genopos_exact.append(snp) # progressive update of geno excat positions 
						else:
							subj_geno[sbjid][chr] = [ datadict['geno'][snp] ]
							genopos_exact.append(snp) # progressive update of geno excat positions 
				## case 2: len geno is twice len markers -> reason: split CSV by \t on the paired alleles
				elif len(datadict['geno']) == len(maplist)*2:
					# now parse gnotypes and cut ped -> create key of chr with snp array
					for snp in genopos:
						if subj_geno[sbjid].has_key(chr):
							subj_geno[sbjid][chr].append( datadict['geno'][snp*2] )
							genopos_exact.append(snp*2) # progressive update of geno excat positions 
							subj_geno[sbjid][chr].append( datadict['geno'][snp*2+1] )
							genopos_exact.append(snp*2+1) # progressive update of geno excat positions 
						else:
							subj_geno[sbjid][chr] = [ datadict['geno'][snp*2], datadict['geno'][snp*2+1] ]
							genopos_exact.append(snp*2) # progressive update of geno excat positions 
							genopos_exact.append(snp*2+1) # progressive update of geno excat positions 
				else:
					print "Error, genotypes have different len: len(datadict['geno']) != len(maplist)::", len(datadict['geno']), len(maplist)
					#sys.exit()
					
				# update genopositions of this chr in the relative structure
				if subj_geno[sbjid].has_key('genopos'):
					if subj_geno[sbjid]['genopos'].has_key(chr):
						subj_geno[sbjid]['genopos'][chr] = genopos_exact
					else:
						subj_geno[sbjid]['genopos'] = {chr: genopos_exact}
				else:
					subj_geno[sbjid]['genopos'] = {chr: genopos_exact}
				
				# update markerpositions of this chr in the relative structure
				if subj_geno[sbjid].has_key('markerpos'):
					if subj_geno[sbjid]['markerpos'].has_key(chr):
						subj_geno[sbjid]['markerpos'][chr] = markerpos
					else:
						subj_geno[sbjid]['markerpos'] = {chr: markerpos}
				else:
					subj_geno[sbjid]['markerpos'] = {chr: markerpos}
				
				# update markerlist of this chr in the relative structure
				if subj_geno[sbjid].has_key('markerlist'):
					if subj_geno[sbjid]['markerlist'].has_key(chr):
						subj_geno[sbjid]['markerlist'][chr] = markerlist
					else:
						subj_geno[sbjid]['markerlist'] = {chr: markerlist}
				else:
					subj_geno[sbjid]['markerlist'] = {chr: markerlist}
					
			else:
				print sbjid, "-> No genotypes for this ID"
	
	def acquirePed(ped, subj_geno):
		"""
		subj_geno = {} -> this is the reference structure! composition:
			k = subject id
			v = dictionary:
				k = geno: genotype data, all markers, all chr -> v = list of markers
					fam: familiar data -> v = list of data
					chr 1..N: slice of genotypes for the chr key -> v = list of markers
					genopos: relative positions of slices for each chr -> v = dictionary:
						k = chr -> v = list of positions relative to geno markers
					markerpos: list of positions of only filtered markers for this chr
						k = chr -> v = list of positions of markers refered to the overall map file previously acquired as tuple array
					markerlist: list of ids
						k = chr -> v = tuple array
		"""
		with open(ped, 'rb') as inf:
			reader = csv.reader(inf, delimiter = "\t")
			for row in reader:
				if len(row) > 6: # -> add genotypes
					print "this file has genotypes, acquire them all for subject", row[1]
					## old not valid since presence of NULL data -> you must parse columns...
					#if subj_geno.has_key(row[1]):
						# subj_geno[row[1]]['geno'] = row[6:] # -> the right solution 
					#else:
						# subj_geno[row[1]] = {'geno': row[6:], }
					
					subj_geno[row[1]] = {'geno': []}
					for m in row[6:]:
						if m is not None and m != '':
							subj_geno[row[1]]['geno'].append(m.strip())
				
				elif len(row) > 1: # -> add familiar data
					print "this file has familiar data, acquire them all for subject", row[1]
					if subj_geno.has_key(row[1]):
						subj_geno[row[1]]['fam'] = row
					else:
						subj_geno[row[1]] = {'fam': row, }
				else:
					print "Skipping subject info: no match len."
					pass
	
	def writePed(ped, subj_geno, what = 'geno'):
		# subj_geno = {}
		# what is or geno or the chr number as int
		with open(ped, 'wb') as inf:
			writer = csv.writer(inf, delimiter = "\t")
			genolen = None
			# understand what to write and write it
			print "Writing data of::", what
			# get geno len to un-genotypes (all 0s)
			for sbjid, datadict in subj_geno.iteritems():
				if datadict.has_key(what):
					genolen = len(datadict[what])
				#print "\tgenolen:", genolen, ". keys:", datadict.keys()
			# write data (also for missing genotypes)
			if genolen is not None:
				print "\tnumber of markers:", genolen
				for sbjid, datadict in subj_geno.iteritems():
					print "Writing data for subject", sbjid
					if datadict.has_key(what):
						writer.writerow(datadict['fam'] + datadict[what])
					else:
						writer.writerow(datadict['fam'] + [0]*(genolen) )
			else:
				print "No markers!"
	
	def writeMapDat(map, dat, subj_geno, target_chr, diseasecol = 0 ):
		# input subj_geno
		## OBSOLETE: maplist_sorted = sorted(maplist, key=itemgetter(0, 2), reverse = False) # just to be safe
		## first: discovery what subject has got genotypes and acquire position for this chr (if available!)
		chrmarkerlist = None # tuple array
		for sbjid, datadict in subj_geno.iteritems():
			if datadict.has_key('markerlist'):
				if datadict['markerlist'].has_key(target_chr):
					chrmarkerlist = datadict['markerlist'][target_chr] ## NB: IDENTICAL for all subjects in this implementation!!!!!!!
				else:
					print "Missing CHR in marker positions -> Run first extractChr to acquire them all and re-run this method"
				pass
		## second: now write files
		## write MAP
		print "Wrinting MAP"
		with open(map, 'wb') as inf:
			writer = csv.writer(inf, delimiter = "\t")
			## OBSOLETE: maplist_sorted = sorted(maplist, key=itemgetter(0, 2), reverse = False) # just to be safe
			for marker in chrmarkerlist:
				chr, snp, cm, markerindex, genoindex = marker
				writer.writerow( [chr, snp, cm] )
		## write DAT
		print "Wrinting DAT"
		with open(dat, 'wb') as inf:
			writer = csv.writer(inf, delimiter = "\t")
			## OBSOLETE: maplist_sorted = sorted(maplist, key=itemgetter(0, 2), reverse = False) # just to be safe
			if diseasecol == 0:
				writer.writerow( ["A", "Condition"] )
			for marker in chrmarkerlist:
				chr, snp, cm, markerindex, genoindex = marker
				writer.writerow( ["M", snp] )
			if diseasecol > 0:
				writer.writerow( ["A", "Condition"] )		


class MemoryMonitor(object):
	def __init__(self, username):
		"""Create new MemoryMonitor instance."""
		self.username = username

	def usage(self):
		"""Return int containing memory used by user's processes."""
		self.process = subprocess.Popen("ps -u %s -o rss | awk '{sum+=$1} END {print sum}'" % self.username,
					    shell=True,
					    stdout=subprocess.PIPE,
					    )
		self.stdout_list = self.process.communicate()[0].split('\n')
		return int(self.stdout_list[0])
        
class MyDB:
	"""
	Andrea - DB Class
	"""
	def __init__(self, host, user, passwd, schema):
		"""
		Init DB class
		"""
		import MySQLdb
		self.conn = MySQLdb.connect( 
			host = host,
			user = user,
			passwd = passwd,
			db = schema )
		#self.cursor = self.conn.cursor ()
		self.cursor = self.conn.cursor (MySQLdb.cursors.DictCursor)
		
	#cursor = conn.cursor ()
	
		
	def close(self):
		"""
		Destroy all
		"""
		self.cursor.close()
		self.conn.close()
		
	def execute(self, query):
		try:
			self.cursor.execute("""
				%s
				""" %(query)
				)
		except MySQLdb.Error, e:
			print "DATABASE Error %d: %s" % (e.args[0], e.args[1])
			sys.exit (1)

class FastASequence:
	"""
	Basic class for sequence fasta management
	"""
	def __init__(self,):
		"""
		Initialize object FastA
		"""
		self.cursor = None
		self.schema = None
		self.table = None
		if self.cursor is not None and self.schema is not None and self.table is not None:
			self.initTable(self.cursor, self.schema, self.table)
		
		self.sequence_dictionary = {}
		
		pass
	
	def importAllSequencesFromFile(self, infile):
		"""
		input: fasta file, unix format
		"""
		handle = open(infile, "rU")
		self.sequence_dictionary = SeqIO.to_dict(SeqIO.parse(handle, "fasta")) # seq dictionary
		handle.close()
	
	def initTable(self, ):
		"""
		Create table of FastA sequences
		"""
		query = """
			CREATE TABLE IF NOT EXISTS `%(schema)s`.`%(table)s` (
			  `experiment_id` varchar(255) NOT NULL,
			  `result_id` int(255) NOT NULL AUTO_INCREMENT,
			  `fasta_id_query` varchar(255) NOT NULL COMMENT 'fasta_id',
			  `db_file` varchar(255) NOT NULL COMMENT 'file name',
			  `fasta_id_match` varchar(255) NOT NULL COMMENT 'fasta_id matched by query sequence',
			  `blast_details` text NOT NULL,
			  `alm_title` varchar(255) NOT NULL,
			  `alm_score` double NOT NULL,
			  `alm_bits` float NOT NULL,
			  `alm_expect` double NOT NULL,
			  `alm_identities` int(11) NOT NULL,
			  `alm_positives` varchar(255) NOT NULL,
			  `alm_gaps` varchar(255) NOT NULL,
			  `alm_align_length` int(25) NOT NULL,
			  `alm_strand` varchar(50) NOT NULL,
			  `alm_frame` varchar(150) NOT NULL,
			  `alm_query` varchar(5000) NOT NULL,
			  `alm_query_start` int(50) NOT NULL COMMENT '1-based',
			  `alm_query_end` int(50) NOT NULL COMMENT '1-based',
			  `alm_match` text NOT NULL,
			  `alm_sbjct` varchar(5000) NOT NULL,
			  `alm_sbjct_start` int(50) NOT NULL COMMENT '1-based',
			  `alm_sbjct_end` int(50) NOT NULL COMMENT '1-based',
			  PRIMARY KEY (`result_id`),
			  KEY `fasta_id_query` (`fasta_id_query`),
			  KEY `fasta_id_match` (`fasta_id_match`)
			) ENGINE=MyISAM  DEFAULT CHARSET=latin1 COMMENT='Results of BLASTN runs.' AUTO_INCREMENT=0 ;
		"""  %{
			'schema': self.schema,
			'table': self.table,
		      }
		self.cursor.execute(query)
 
class SequenceAnalysis:
	"""
	Utils for fasta handling
	"""
	def __init__(self, ):
		print "[AP]\tInitializing FastaSequence object"
	
	
	def getSetOfSimilarityElems(self, f_fasta_id, cursor, nucleothreshold, f_inputset_of_fasta_ids, blastnexpectthreshold):
		"""
		Query DB matrix of Blastn for similarities and returns a set of fasta_id.
		"""
		
		###### search for similarities in the matrix blastn
		query = Template("""
			SELECT  `fasta_id_query` ,  `fasta_id_match` ,  `alm_score` ,  `alm_bits` ,  `alm_expect` ,  `alm_identities` ,  `alm_positives` ,  `alm_align_length` ,   `alm_query_start` ,  `alm_query_end` ,  `alm_match` ,  `alm_sbjct` ,  `alm_sbjct_start` ,  `alm_sbjct_end` 
			FROM  `blastn` 
			WHERE 
			  (`fasta_id_query` LIKE  '$SEQFASTAID' OR  `fasta_id_match` LIKE  '$SEQFASTAID' )
			  AND 
			  ( `alm_query_start` <= $NUCLEOTHRESHOLD AND  `alm_sbjct_start` <= $NUCLEOTHRESHOLD )
			  AND
			  `alm_expect` <= $EXPECTBLASTNFLOAT
			ORDER BY  `blastn`.`alm_expect` ASC 
			""")
		outset = set() # final set
		outset.add( f_fasta_id ) # to be safe in case of empty set from query (even with the same sequence!!! because len(seq) = 20), add it
		q = query.substitute( SEQFASTAID = f_fasta_id, NUCLEOTHRESHOLD = nucleothreshold, EXPECTBLASTNFLOAT = float(blastnexpectthreshold))
		cursor.execute( q )
		similarity_elems = cursor.fetchall()
		for elem in similarity_elems:
			###### add fasta_id value in the cluster set iff present in the input set
			if elem[0] in f_inputset_of_fasta_ids:
				outset.add(elem[0]) # add elem query
			#else:
				#print "\tNote: found a similarity value in the matrix with a sequence that is not present in your input dataset: ", elem[0]
			if elem[1] in f_inputset_of_fasta_ids:
				outset.add(elem[1]) # add elem db
			#else:
				#print "\tNote: found a similarity value in the matrix with a sequence that is not present in your input dataset: ", elem[1]
		return outset
		
	def sortArrayOfFastaIdBySequenceLen(self, array_fasta_id, dic_seqs, reverse_option):
		"""
		From an input array of fasta_id (belonging to input sequences), sort fasta_id by sequence len using dictionary.
		Reverse parameter is boolean, if True, then sort DESCENDING. Else, ASCENDING.
		"""
		ret_array = []
		tmp_array = []
		tuple_array = []
		index = 0
		for k in array_fasta_id:
			k_seq_len = dic_seqs[k][1].__len__()
			tmp_array.append( "%s_%s" %(k_seq_len, k) )
			elem = k_seq_len, k
			tuple_array.append( elem )
		
		#print tuple_array, reverse_option
		tmp_array_sorted = sorted(tuple_array, key=itemgetter(0), reverse = int(reverse_option))
		for k in tmp_array_sorted:
			ret_array.append( k[1] )
		
		return ret_array
		
	def getArrayOfElemFromFastaIdArray(self, array_fasta_id, dic_seqs, query_position_index_0based):
		"""
		From an input array of fasta_id (belonging to input sequences), get array of field values in the same order. (=> like a table)
		"""
		ret_array = []
		for k in array_fasta_id:
			k_field = dic_seqs[k][query_position_index_0based]
			ret_array.append( k_field )
		return ret_array
		

class BlastN:
	def __init__(self, ):
		print "[AP]\tInitializing BlastN object"
	
	def indexDb(self, blastndb):
		"""
		Index blastn reference file -> the input one or the specified db!
		Returns bool
		"""
		if os.path.exists(blastndb):
			os.system("formatdb -i %s -p F" %(blastndb) )
			return True
		else:
			return False

	def runBlastn(self, blastnprogram, query, blastndb, blastouttmpxml, blastnidentity, blastncustomstring):
		"""
		Run a single (!) BLASTN thread.
		Uses NON global call to /usr/local/bin/blastn
		I.E.:
		blastn -query POOL9_CCAAGG.01seq.fa -db /home/andrea/Dropbox/TIGET/Input/POOL9/POOL9_CCAAGG.fa -out MY_Test_Out_01vsAll.xml -outfmt 5 -perc_identity 95
		"""
		print "[AP]\tBlasting data."
		#i.e. blastn: /opt/shared/apps/ncbi/ncbi-blast-2.2.25+/bin/./blastn
		call = "%(blastnprogram)s -query %(query)s -db %(blastndb)s -out %(blastouttmpxml)s -outfmt 5 -perc_identity %(blastnidentity)s -word_size 20 %(blastncustomstring)s" 	%{	
				'query': query,
				'blastnprogram': blastnprogram,
				'blastndb': blastndb,
				'blastouttmpxml': blastouttmpxml,
				'blastnidentity': blastnidentity,
				'blastncustomstring': blastncustomstring,
			}
		#print call
		os.system( call )
	
	def populateDatabase(self, infile):
		"""
		input: fasta file, unix format
		"""
		fastaSeqs = FastaSequence()
		fastaSeqs.importAllSequencesFromFile(infile)
	
	def initTable(self, ):
		"""
		Create table of FastA sequences
		"""
		query = """
			CREATE TABLE IF NOT EXISTS `%(schema)s`.`%(table)s` (
			  `experiment_id` varchar(255) NOT NULL,
			  `result_id` int(255) NOT NULL AUTO_INCREMENT,
			  `fasta_id_query` varchar(255) NOT NULL COMMENT 'fasta_id',
			  `db_file` varchar(255) NOT NULL COMMENT 'file name',
			  `fasta_id_match` varchar(255) NOT NULL COMMENT 'fasta_id matched by query sequence',
			  `blast_details` text NOT NULL,
			  `alm_title` varchar(255) NOT NULL,
			  `alm_score` double NOT NULL,
			  `alm_bits` float NOT NULL,
			  `alm_expect` double NOT NULL,
			  `alm_identities` int(11) NOT NULL,
			  `alm_positives` varchar(255) NOT NULL,
			  `alm_gaps` varchar(255) NOT NULL,
			  `alm_align_length` int(25) NOT NULL,
			  `alm_strand` varchar(50) NOT NULL,
			  `alm_frame` varchar(150) NOT NULL,
			  `alm_query` varchar(5000) NOT NULL,
			  `alm_query_start` int(50) NOT NULL COMMENT '1-based',
			  `alm_query_end` int(50) NOT NULL COMMENT '1-based',
			  `alm_match` text NOT NULL,
			  `alm_sbjct` varchar(5000) NOT NULL,
			  `alm_sbjct_start` int(50) NOT NULL COMMENT '1-based',
			  `alm_sbjct_end` int(50) NOT NULL COMMENT '1-based',
			  PRIMARY KEY (`result_id`),
			  KEY `fasta_id_query` (`fasta_id_query`),
			  KEY `fasta_id_match` (`fasta_id_match`)
			) ENGINE=MyISAM  DEFAULT CHARSET=latin1 COMMENT='Results of BLASTN runs.' AUTO_INCREMENT=0 ;
		"""  %{
			'schema': self.schema,
			'table': self.table,
		      }
		self.cursor.execute(query)

class StatAnalysis:
	"""
	A class for math transformations
	"""
	def __init__(self, ):
		#import math, numpy, scipy.stats
		print "[IL]\tNew object for analysis."
		self.default_val_for_negative_infinite = -99999999999999999
		self.default_val_for_positive_infinite = 99999999999999999
	
	def log2(self, val):
		"""
		Given an input score, compure LOGbase2(x/1-x))
		since logb(a) = ln(a)/ln(b), valid domain interval is (0, 1)
		"""
		if val == 0:
			### val = default_val_for_0
			return self.default_val_for_negative_infinite
		elif val == 1:
			### val = 0.9999999999999999 # no more!! it will not work! because computer round it at 1 adding only a single other 9
			return self.default_val_for_positive_infinite
		else:
			score = float(val)/ (1-val) # score is float and represents inner function
			return math.log( score ) / math.log( 2 )
	
	def log2xon1x(self, values_list):
		"""
		Given an input score, compure LOGbase2(x/1-x))
		since logb(a) = ln(a)/ln(b), valid domain interval is (0, 1)
		"""
		log2 = []
		for val in values_list:
			if val == 0:
				### val = default_val_for_0
				log2.append( self.default_val_for_negative_infinite )
			elif val == 1:
				### val = 0.9999999999999999 # no more!! it will not work! because computer round it at 1 adding only a single other 9
				log2.append( self.default_val_for_positive_infinite )
			else:
				score = float(val)/ (1-val) # score is float and represents inner function
				log2.append( math.log( score ) / math.log( 2 ) )
		return log2
	
	def zScore(self, values_list):
		"""
		Given a list of values in input, compute the Z-Score among data.
		Returns a list of values (same len of input) in which 
		"""
		zscore_list = []
		avg = float( sum(values_list) )/len(values_list)
		stddev = np.std( values_list )
		for v in values_list:
			zscore_list.append( (float(v) - avg ) / float(stddev) )
		return zscore_list
	
	######################################################
	### Z SCORE (from http://telliott99.blogspot.it/2010/01/z-scores.html)
	######################################################
	def z_score(self, row):
		## for z score computation
		from numpy import nan, isnan
		from numpy import array, mean, std, random
		L = [n for n in row if not isnan(n)]
		m = mean(L)
		s = std(L)
		zL = [1.0 * (n - m) / s for n in L]
		if len(L) == len(row):  return zL
		# deal with nan
		retL = list()
		for n in row:
			if isnan(n):
				retL.append(nan)
			else:
				retL.append(zL.pop(0))
		assert len(zL) == 0
		return retL
	
	def z_scores_by_row(self, A):
		## for z score computation
		from numpy import nan, isnan
		from numpy import array, mean, std, random
		retL = [z_score(r) for r in A]
		Z = array(retL)
		Z.shape = A.shape
		return Z
	######################################################
	
	def percentage_score(self, row):
		from numpy import nan, isnan
		from numpy import array, mean, std, random
		L = [n for n in row if not isnan(n)]
		s = sum(L)
		pL = [100.0 * (float(n)/s) for n in L]
		if len(L) == len(row):  return pL
		# deal with nan
		retL = list()
		for n in row:
			if isnan(n):
				retL.append(nan)
			else:
				retL.append(zL.pop(0))
		assert len(pL) == 0
		return retL

class BioAnnotation:
	"""
	General class for annotation. Needs Bio.Entrez and SeqIO
	"""
	def __init__(self, ):
		print "[IL]\tNew object for bio-annotation."
		self.chrom = None
		self.bpStart = None
		self.bpEnd = None
		self.strand = 1 # 1 or -1
		self.myemail = "semperme@example.com"
		self.sleeptime = 1
		self.outFormat = "fasta" # default
		self.seqRecord = None
		self.id = 9606
	
	def sliceGenomicRegion(self, id):
		"""
		returns object from Entrez
		"""
		if self.bpStart is not None and self.bpEnd is not None:
			Entrez.email = self.myemail    # Always tell NCBI who you are
			handle = Entrez.efetch(	db = "nucleotide", 
						id = 9606, ###"307603377", 
						rettype = self.outFormat, 
						strand = self.strand, 
						seq_start = self.bpStart, 
						seq_stop = self.bpEnd)
			record = SeqIO.read(handle, self.outFormat)
			handle.close()
			self.seqRecord = record
		else:
			print "[AP]\tBefore slicing genome, update object boundaries (start - end) and ID"
	
	def fetchGenomicSequence(self, ):
		"""
		returns genomics bases of a region from Entrez
		"""
		if self.seqRecord is not None:
			return self.seqRecord.seq
			
			

class IntegrationLocus:
	"""
	Class of integration loci
	"""
	def __init__(self, ):
		print "[IL]\tCreating Integration Locus object"
		self.dict_loci = {}
		self.chromosome = None
		self.locus_start = None
		self.locus_end = None
		
	
	def checkCollision(self, cursor, schema, table, chromosome, locus):
		"""
		Look for all sample data integratino for contaminations
		example of sql:
			SELECT * 
			FROM  MLD03.`no_redundant_MLD03_FREEZE_SIX_MONTHS` 
			WHERE  `chr` =  '1'
			      AND  `integration_locus` 
			      BETWEEN 1044619 
			      AND 1044625 
		"""
		#print "[DB]\tChecking locus position for collisions"
		outdata = 0 
		this_table = ""
		if table.startswith('no_'):
			this_table = '_'.join( str(x) for x in table.split('_')[1:] )
		else:
			this_table = table
		# per prendere i nomi UNIVOCI delle LAM devi passare per la REDUNDANT, invece che dalle NO:REDUNDANT (come era versione 1)
		query = """
			SELECT count( * ) as number_collisions
			FROM  `%s`.`%s`
			WHERE  `chr` =  '%s'
			      AND  `integration_locus` BETWEEN '%d' AND '%d' 
			;
			""" %(	schema, 
				this_table,
				chromosome,
				int(locus) - 3,
				int(locus) + 3,
				)
		#print query
		cursor.execute(query)
		results = cursor.fetchall()
		for r in results:
			outdata += int(r["number_collisions"])
		return outdata
	
	def getExon(self, cursor, schema, table):
		"""
		Annotation function: retrieve exon for this IL. If not available, it returns intron data
		Database data deal with ANNOTATION DB, not integration DB. i.e.: schema = ucsc_hg19; table = refGene
		"""
		query_template = Template( """
			SELECT `name` ,  `chrom` ,  `strand` ,  `txStart` ,  `txEnd` ,  `cdsStart` ,  `cdsEnd` ,  `exonCount` ,  `exonStarts` ,  `exonEnds` ,  `score` ,  `name2` 
			FROM  `%(schema)s`.`%(table)s`
			WHERE  `chrom` LIKE  '$chrom'
			  AND  `txStart` <= '$txStart' AND txEnd >= '$txEnd'
			  ORDER BY  `refGene`.`exonCount` DESC 
			  LIMIT 1
			  ;
			""" %{
			      'schema': schema, 
			      'table': table,
			      }
			)
		self.dict_loci = {}
		exon_annotation = []
		if not self.chromosome.startswith('chr'):
			self.chromosome = 'chr' + chrom
		
		# query for this gene
		query = query_template.substitute(chrom = self.chromosome, txStart = int(self.locus_start), txEnd = int(self.locus_end))
		cursor.execute( query )
		results = cursor.fetchall()
		if len(results) > 0:
			self.dict_loci[(chrom, locus_start, locus_end, locus_id)] = results[0]
		for res in results:
			reskeys = res.keys()
			#print "\t-->", res['name2'], row, res
			for exon in range(0, int(res['exonCount']) ):
				if self.locus_start >= int(res['exonStarts'].split(',')[exon]) -3 and self.locus_end <= int(res['exonEnds'].split(',')[exon]) +3:
					#you have a valid exon
					#now check for strand and reverse count
					if res['strand'] == '+':
						return row.strip().split()[0:4] + ['Exon', exon + 1, res['exonStarts'].split(',')[exon], res['exonEnds'].split(',')[exon], res['name2'], res['strand'], res['txStart'], res['txEnd'], res['exonCount']]
					elif res['strand'] == '-':
						reverse_strand_exon_index = int(res['exonCount']) - (exon + 1)
						return row.strip().split()[0:4] + ['Exon', reverse_strand_exon_index, res['exonStarts'].split(',')[reverse_strand_exon_index], res['exonEnds'].split(',')[reverse_strand_exon_index], res['name2'], res['strand'], res['txStart'], res['txEnd'], res['exonCount']]

	
	
	def close(self):
		"""
		Destroy all
		"""
		pass
		
	def query(self, cursor, schema, table, query):
		try:
			#print query
			cursor.execute("""
				%s
				""" %(query)
				)
		except MySQLdb.Error, e:
			print "DATABASE Error %d: %s" % (e.args[0], e.args[1])
			sys.exit (1)

class Experiment:
	"""
	Questa e' la classe principale per gli esperimenti.
	Ci sono le funzioni per le no redundant e altre manipolazioni dati.
	
	Aggiornamento a giugno 2012:
		ho creato un nuovo approccio di uso dati: nella funzione getUniqueIntegrationLociImproved genero le strutture dati primario integrationLoci_by_TST_uniqueDetails che contengono come chiave l'IS univoco riportato come no redundant e come valori i dettagli che caratterizzano quel is con anche i loci, distinti per ogni TST (tissue, sample, timepoint).
		Da qui e' possibile fare ogni genere di operazione:
			- filtraggio IS da rimuovere [filterUniqueIntegrationSites]
			- accorpamento di elementi TST [groupElementsOfTST]
			- creazione dei valori di summary con aggiunta di zscore e percentuali
		La cosa importante e' che con questa struttura e' sempre possibile risalire al singolo IS redundant e al suo riferimento no redundant; anche la struttura integrationLoci_by_TST_wrtRed serve a questo scopo: traduce l'IS no redundant negli altri redundant (memorizzato per ogni TST).
	"""
	def __init__(self, ):
		self.name = ""
		self.integrationLoci = {}
		self.integrationLoci_data = ()
		self.lam_data = {} # k = complete_name, v = dictionary where k1 = column name, v1 = column value
		self.sequencecount_by_genes = {} # k = patient, v = dict:: k = sample, v = dict:: k = gene name refSeqSymbol, 'allgenes'; v = dict:: k = {seqcount, treatment, tissue, chr, orientation, lenght, agv_position }
		self.sequencecount_by_iss_genes = {} # k = patient, v = dict:: k = sample, v = dict:: k = (gene name refSeqSymbol, chr, position); v = dict:: k = {seqcount, treatment, tissue, chr, orientation, lenght, agv_position }
		self.integrationLoci_by_TS = {} # k = patient, v = dict:: k = (tissue, sample), v = dict:: k = (chr, position); v = dictionary of locus data
		self.integrationLoci_by_TS_only = {} # k = (tissue, sample), v = dict:: k = (chr, position); v = dictionary of locus data
		self.integrationLoci_by_TST_wrtRed = {} # k = (tissue, sample, timepoint), v = dict:: k = (chr, position) from cell line; v = (chr, position) of reference IS from no-redundant
		self.integrationLoci_by_TST_uniqueDetails = {} # k = (tissue, sample, timepoint), v = dict:: k = (chr, position) of reference IS from no-redundant; v = dict:: details with specific keys: 'loci' as tuple of (chr, locus), 'refSeqSymbol' as gene name, 'sequence_count' = sequence_count
		
	def getDistinctNamesFromDB(self, cursor, schema, table, field):
		"""
		Look for all sample data in the field SAMPLE of the no_redundant table.
		example of sql:
			SELECT DISTINCT  `sample` 
			FROM  `MLD01`.`no_redundant_MLD01_ALL_NEW_NAMES` 
		"""
		print "[DB]\tRetrieving all distinct %s names from DB" %(field)
		outdata = set()
		query = """
			SELECT DISTINCT  `%s` 
			FROM  `%s`.`%s`;
			""" %(	field, schema, table )
		
		cursor.execute(query)
		results = cursor.fetchall()
		for r in results:
			outdata.add(r[field])
		return outdata


	def getUniqueIntegrationLoci(self, cursor, schema, table, select, keysBedFormat = False):
		"""
		Look for all real IS data:
		- get sorted data -> sorted tuple
		- for each IS, count +- 3 rule and create a dictionary: given N integrations for a single IS, k = chr, locus; v = tuple of N integrations
		
		Alternative query:
			SELECT DISTINCT `chr` , `integration_locus` , `sample` , `tissue` , `treatment` , sum( sequence_count ) AS sequence_count
			FROM sequence_was1001.`redundant_WAS1001_FREEZE_2012_newTags`
			WHERE chr != '0'
			GROUP BY `chr` , `integration_locus`
			ORDER BY `chr` , `integration_locus`
		"""
		print "[DB]\tComputing unique integration loci extraction (+/- 3 bp) from DB"
		dict_is = {}
		query = """
			SELECT %(select)s 
			FROM  `%(schema)s`.`%(table)s`
			WHERE `chr` not like '0'
			ORDER BY  `%(table)s`.`chr` ASC,  `%(table)s`.`integration_locus` ASC 
			-- LIMIT 1000
			;
			""" %{
			      'select': select,
			      'schema': schema, 
			      'table': table
			      }
#		NEW!!	
# 		query = """
# 			SELECT DISTINCT `chr` , `integration_locus` , `sample` , `tissue` , `treatment` , sum( sequence_count ) AS sequence_count, refSeqSymbol
# 			FROM  `%(schema)s`.`%(table)s`
# 			WHERE `chr` not like '0'
# 			GROUP BY `%(table)s`.`chr` , `%(table)s`.`integration_locus`
# 			ORDER BY  `%(table)s`.`chr` ASC,  `%(table)s`.`integration_locus` ASC 
# 			-- LIMIT 1000
# 			;
# 			""" %{
# 			      'select': select,
# 			      'schema': schema, 
# 			      'table': table
# 			      }
			      
		cursor.execute(query)
		results = cursor.fetchall()
		self.integrationLoci_data = results
		
		integration_locus = None
		for i in range(0,len(results)):
			# se la diff tra integration loci e' <= 3 E l'attuale locus in analisi ha una differenza <= 3 con il primo allora, inseriscilo nel dizionario
			if i > 0 and integration_locus is not None and ( abs(results[i]["integration_locus"] - results[i-1]["integration_locus"]) <= 3 and abs(integration_locus - results[i]["integration_locus"]) <= 3 ):
				# check for keys layout (bed format or not) -> default is not related to strand!!!!!!!
				if keysBedFormat:
					dict_is[( str(results[i]["chr"]), integration_locus, integration_locus+1 )].append( results[i] )
				else:
					dict_is[( str(results[i]["chr"]), integration_locus )].append( results[i] )
				
			else:
				if keysBedFormat:
					integration_locus = int(results[i]["integration_locus"])
					dict_is[( str(results[i]["chr"]), integration_locus, integration_locus+1 )] = [results[i]]
				else:
					integration_locus = int(results[i]["integration_locus"])
					dict_is[( str(results[i]["chr"]), integration_locus )] = [results[i]]
			
		return dict_is, results
	
	def getUniqueIntegrationLociImproved(self, cursor, schema, table, select, keysBedFormat = False, ):
		"""
		Look for all real IS data:
		- get sorted data -> sorted tuple
		- for each IS, count +- 3 rule and create a dictionary: given N integrations for a single IS, k = chr, locus; v = tuple of N integrations
		 
		Output:
		- dict_is: k = IS (chr, locus); v = dict:: k = 'loci' as tuple of (chr, locus), 'refSeqSymbol' as gene name, 'sequence_count' = sequence_count
		
		minimum fields: `integration_locus` , `sample` , `tissue` , `treatment` , strand, sum( sequence_count ) AS sequence_count, refSeqSymbol, orientation_ucsc 
		
		Alternative query:
			SELECT DISTINCT `chr` , `integration_locus` , `sample` , `tissue` , `treatment` , sum( sequence_count ) AS sequence_count
			FROM sequence_was1001.`redundant_WAS1001_FREEZE_2012_newTags`
			WHERE chr != '0'
			GROUP BY `chr` , `integration_locus`
			ORDER BY `chr` , `integration_locus`
		"""
		print "[DB]\tComputing unique integration loci extraction (+/- 3 bp) from DB"
		dict_is = {}
		query = """
			-- SELECT DISTINCT `chr` , `integration_locus` , `sample` , `tissue` , `treatment` , sum( sequence_count ) AS sequence_count, refSeqSymbol
			-- FROM  `%(schema)s`.`%(table)s`
			-- WHERE `chr` not like '0'
			-- GROUP BY `%(table)s`.`chr` , `%(table)s`.`integration_locus`, `treatment`
			-- ORDER BY  `%(table)s`.`chr` ASC,  `%(table)s`.`integration_locus` ASC 
			
			SELECT %(select)s 
			FROM  `%(schema)s`.`%(table)s`
			WHERE `chr` not like '0'
			ORDER BY  `%(table)s`.`chr` ASC,  `%(table)s`.`integration_locus` ASC 
			-- LIMIT 10000
			;
			""" %{
			      'select': select,
			      'schema': schema, 
			      'table': table
			      }
			      
		cursor.execute(query)
		results = cursor.fetchall()
		
		# visto che self.integrationLoci_by_TST_uniqueDetails e self.integrationLoci_by_TST_wrtRed HANNO IDENTICHE STRUTTURE ma differenti valori finali, posso compilarli in parallelo
		
		integration_locus = None
		for i in range(0,len(results)):
			# se la diff tra integration loci e' <= 3 E l'attuale locus in analisi ha una differenza <= 3 con il primo allora, inseriscilo nel dizionario
			if i > 0 and integration_locus is not None and ( abs(results[i]["integration_locus"] - results[i-1]["integration_locus"]) <= 3 and abs(integration_locus - results[i]["integration_locus"]) <= 3 ):
				# check for keys layout (bed format or not) -> default is not related to strand!!!!!!!
				if keysBedFormat:
					#dict_is[( str(results[i]["chr"]), integration_locus, integration_locus+1 )].append( results[i] )
					
					# update a global dictionary of data that pairs single experiment IS with reference no-redundant IS
					if self.integrationLoci_by_TST_wrtRed.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						# update self.integrationLoci_by_TST_* structures
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1)] = (str(results[i]["chr"]), integration_locus, integration_locus+1)
					else:
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
								(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1): (str(results[i]["chr"]), integration_locus, integration_locus+1)
							}
				else:
					#dict_is[( str(resuts[i]["chr"]), integration_locus )].append( results[i] )
					
					# update a global dictionary of data that pairs single experiment IS with reference no-redundand IS
					if self.integrationLoci_by_TST_wrtRed.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), int(results[i]["integration_locus"]))] = (str(results[i]["chr"]), integration_locus)
					else:
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
								(str(results[i]["chr"]), int(results[i]["integration_locus"])): (str(results[i]["chr"]), integration_locus)
							}
				
			else: # questo else avvia il conto di OGNI NUOVO IS, quindi e' il primo ad essere passato!!!!
				if keysBedFormat:
					integration_locus = int(results[i]["integration_locus"])
					#dict_is[( str(results[i]["chr"]), integration_locus, integration_locus+1 )] = [results[i]]
					
					# update a global dictionary of data that pairs single experiment IS with reference no-redundand IS
					if self.integrationLoci_by_TST_wrtRed.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1)] = (str(results[i]["chr"]), integration_locus, integration_locus+1)
					else:
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
							(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1): (str(results[i]["chr"]), integration_locus, integration_locus+1)
							}
				else:
					integration_locus = int(results[i]["integration_locus"])
					#dict_is[( str(results[i]["chr"]), integration_locus )] = [results[i]]
					
					# update a global dictionary of data that pairs single experiment IS with reference no-redundand IS
					if self.integrationLoci_by_TST_wrtRed.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), int(results[i]["integration_locus"]))] = (str(results[i]["chr"]), integration_locus)
					else:
						self.integrationLoci_by_TST_wrtRed[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
							(str(results[i]["chr"]), int(results[i]["integration_locus"])): (str(results[i]["chr"]), integration_locus)
							}
							
		integration_locus = None
		for i in range(0,len(results)):
			#print results[i]["tissue"], results[i]["sample"], results[i]["treatment"], results[i]["chr"], results[i]["integration_locus"]
			# se la diff tra integration loci e' <= 3 E l'attuale locus in analisi ha una differenza <= 3 con il primo allora, inseriscilo nel dizionario come reference
			if i > 0 and integration_locus is not None and ( abs(results[i]["integration_locus"] - results[i-1]["integration_locus"]) <= 3 and abs(integration_locus - results[i]["integration_locus"]) <= 3 ): # se non sei il primo elemento e la differenza tra te e il precedente <=3, ti assegno alla reference
				#print "\t...lo assegno al precendete IS (ref) al locus:", integration_locus
				# check for keys layout (bed format or not) -> default is not related to strand!!!!!!!
				
				# se esiste la chiave superiore (TST), verifica che il locus esista gia'
				# 	se questo stesso locus c'e' gia', incrementa sequence count
				# 	altrimenti se non c'è allora per questo TST, inseriscilo come nuovo con il suo sc relativo, ponendo come chiave il ref locus
				# altrimenti se non esiste la chiave TST, creala inserendo come prima chiave locus la ref 
				
				if keysBedFormat:
					# update a global dictionary of data that pairs single experiment IS with reference no-redundant IS
					if self.integrationLoci_by_TST_uniqueDetails.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						#print "\t\tquesto TST esiste gia', ora valuto il locus [if2, dopo if 1]"
						# aggiorna la struttura dati per ogni TST, dove raccogli le info di sequence_count e geni (details)
						if self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])].has_key(
							(str(results[i]["chr"]), integration_locus, integration_locus+1)
							):
							#print "\t\tquesto locus esiste gia', quindi incremento i valori della ref seq count e loci [if 3, fopo if 2, 1]"
							self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus, integration_locus+1)]['sequence_count'] += int(results[i]['sequence_count'])
							self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus, integration_locus+1)]['loci'].add( (str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1) )
						else:
							#print "\t\tquesto locus e' NOUVO per questo TST, quindi inserisco la ref in TST [else 3, Dopo if 2, 1]"
							self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus, integration_locus+1)] = {
								  'sequence_count': int(results[i]['sequence_count']),
								  'refSeqSymbol': results[i]['refSeqSymbol'],
								  'strand': results[i]['strand'],
								  'gene_orientation': results[i]['orientation_ucsc'],
								  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1)]),
								  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
								  }
					else:
						#print "\t\tquesto TST e' nuovo, inserisco il nuovo TST con la ref locus [if 3, fopo if 2, 1]"
						self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
							(str(results[i]["chr"]), integration_locus, integration_locus+1): 
							  {
							  'sequence_count': int(results[i]['sequence_count']),
							  'refSeqSymbol': results[i]['refSeqSymbol'],
							  'strand': results[i]['strand'],
							  'gene_orientation': results[i]['orientation_ucsc'],
							  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1)]),
							  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
							  }
							}
				else: # else del formato BED
					# update a global dictionary of data that pairs single experiment IS with reference no-redundant IS
					if self.integrationLoci_by_TST_uniqueDetails.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						#print "\t\tquesto TST esiste gia', ora valuto il locus [if2, dopo if 1]"
						# aggiorna la struttura dati per ogni TST, dove raccogli le info di sequence_count e geni (details)
						if self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])].has_key(
							(str(results[i]["chr"]), integration_locus)
							):
							#print "\t\tquesto locus esiste gia', quindi incremento i valori della ref seq count e loci [if 3, fopo if 2, 1]"
							self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus)]['sequence_count'] += int(results[i]['sequence_count'])
							self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus)]['loci'].add( (str(results[i]["chr"]), int(results[i]["integration_locus"]) ) )
						else:
							#print "\t\tquesto locus e' NOUVO per questo TST, quindi inserisco la ref in TST [else 3, Dopo if 2, 1]"
							self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus)] = {
								  'sequence_count': int(results[i]['sequence_count']),
								  'refSeqSymbol': results[i]['refSeqSymbol'],
								  'strand': results[i]['strand'],
								  'gene_orientation': results[i]['orientation_ucsc'],
								  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]))]),
								  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
								  }
					else:
						#print "\t\tquesto TST e' nuovo, inserisco il nuovo TST con la ref locus [if 3, fopo if 2, 1]"
						self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
							(str(results[i]["chr"]), integration_locus): 
							  {
							  'sequence_count': int(results[i]['sequence_count']),
							  'refSeqSymbol': results[i]['refSeqSymbol'],
							  'strand': results[i]['strand'],
							  'gene_orientation': results[i]['orientation_ucsc'],
							  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]))]),
							  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
							  }
							}
			else: # questo else avvia il conto di OGNI NUOVO IS, quindi e' il primo ad essere passato!!!! [else 1]
				# il mio IS di riferimento ora e' questo
				integration_locus = int(results[i]["integration_locus"])
				#print "\tNEW ref IS:", integration_locus
				
				if keysBedFormat:
					# update a global dictionary of data that pairs single experiment IS with reference no-redundand IS
					if self.integrationLoci_by_TST_uniqueDetails.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						#print "\t\tquesto TST esiste gia', quindi inserisco la ref come nuova [else 1]"
						self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus, integration_locus+1)] = {
							  'sequence_count': int(results[i]['sequence_count']),
							  'refSeqSymbol': results[i]['refSeqSymbol'],
							  'strand': results[i]['strand'],
							  'gene_orientation': results[i]['orientation_ucsc'],
							  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1)]),
							  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
							}
					else:
						#print "\t\tquesto TST e' NUOVO, quindi prima creo TST e poi inserisco la ref come nuova [else 2, dentro else 1]"
						self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
							(str(results[i]["chr"]), integration_locus, integration_locus+1): 
							  {
							  'sequence_count': int(results[i]['sequence_count']),
							  'refSeqSymbol': results[i]['refSeqSymbol'],
							  'strand': results[i]['strand'],
							  'gene_orientation': results[i]['orientation_ucsc'],
							  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]), int(results[i]["integration_locus"])+1)]),
							  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
							  }
							}
				else: # else del formato BED
					# update a global dictionary of data that pairs single experiment IS with reference no-redundand IS
					if self.integrationLoci_by_TST_uniqueDetails.has_key( (results[i]["tissue"], results[i]["sample"], results[i]["treatment"]) ):
						#print "\t\tquesto TST esiste gia', quindi inserisco la ref come nuova [else 1]"
						self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])][(str(results[i]["chr"]), integration_locus)] = {
							  'sequence_count': int(results[i]['sequence_count']),
							  'refSeqSymbol': results[i]['refSeqSymbol'],
							  'strand': results[i]['strand'],
							  'gene_orientation': results[i]['orientation_ucsc'],
							  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]))]),
							  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
							}
					else:
						#print "\t\tquesto TST e' NUOVO, quindi prima creo TST e poi inserisco la ref come nuova [else 2, dentro else 1]"
						self.integrationLoci_by_TST_uniqueDetails[(results[i]["tissue"], results[i]["sample"], results[i]["treatment"])] = {
							(str(results[i]["chr"]), integration_locus): 
							  {
							  'sequence_count': int(results[i]['sequence_count']),
							  'refSeqSymbol': results[i]['refSeqSymbol'],
							  'strand': results[i]['strand'],
							  'gene_orientation': results[i]['orientation_ucsc'],
							  'loci': set([(str(results[i]["chr"]), int(results[i]["integration_locus"]))]),
							  'gene_lenght_refseq': results[i]["gene_lenght_refseq"],
							  }
							}
					
							
		#return dict_is, results
	
	def addSummaryCountsToTSTISs(self, summary_field = 'alliss'):
		"""
		Given input structure [integrationLoci_by_TST_uniqueDetails] of ISS, parse them all and extract summaries that will be filled into the same structure. 
		Summary also adds:
			- zscore
			- percentage
		The fields order into this data tuple will be: chr, locus, seqreads, is strand, refSeqSymbol, gene orientation, gene lenght, zscore_seqcount, percentage_sequencecount
		
		NB: use this function only at the end of each filtering procedure that you want to use!!!!!! Only in this way all counts will be consistent! For example, if you have contamination data, remove them all before running this method.
		"""
		tmp_summary_dict = {}
		### step 1: collect data
		for tst, loci in self.integrationLoci_by_TST_uniqueDetails.iteritems():
			# start updating tmp dictionary
			tmp_summary_dict[tst] = []
			
			############## z score #################
			# init vars to compute z score
			zscore_list = [[],[]] # simple array of size = len(dict_is.keys()), made by array of array: [[genes], [sequence count], [zscore - only after processing seqeunce count it is appended]]
			stata = StatAnalysis()
			# process zscore and append values
			for i, det in loci.iteritems():
				zscore_list[0].append(i)
				zscore_list[1].append(det['sequence_count'])
			zscore_outlist = stata.z_score(zscore_list[1])
			zscore_list.append(zscore_outlist)
			# create a new dictionary to store comparable zscores
			dict_is_zscores = {}
			for i in range(0, len(zscore_list[0])):
				dict_is_zscores[zscore_list[0][i]] = zscore_list[2][i]
			if len(set(dict_is_zscores.keys()) - set(loci.keys())) is not 0: # check at least len consistency in keys
				print "[AP]\tError on dictionary keys for Zscore: data not reliable!!"
				sys.exit()
			
			############## percentage #################
			# init vars to compute percentage score
			percscore_list = [[],[]] # simple array of size = len(dict_is.keys()), made by array of array: [[genes], [sequence count], [perc - only after processing seqeunce count it is appended]]
			# process percentage and append values
			for i, det in loci.iteritems():
				percscore_list[0].append(i)
				percscore_list[1].append(det['sequence_count'])
			percscore_outlist = stata.percentage_score(percscore_list[1])
			percscore_list.append(percscore_outlist)
			# create a new dictionary to store comparable zscores
			dict_is_percscores = {}
			for i in range(0, len(percscore_list[0])):
				dict_is_percscores[percscore_list[0][i]] = percscore_list[2][i]
			if len(set(dict_is_percscores.keys()) - set(loci.keys())) is not 0: # check at least len consistency in keys
				print "[AP]\tError on dictionary keys for percentage score: data not reliable!!"
				sys.exit()
				
			############# now insert summary ###########
			# locus, seqreads, is strand, refSeqSymbol, gene orientation, gene lenght, zscore_seqcount, percentage_sequencecount
			for locus, isdata in loci.iteritems():
				tmp_summary_dict[tst].append(
					locus + tuple( [ isdata['sequence_count'], isdata['strand'], isdata['refSeqSymbol'], isdata['gene_orientation'], isdata['gene_lenght_refseq'], dict_is_zscores[locus], dict_is_percscores[locus] ] )
					)
		### step 2: write data into global structure
		for tst, tupledata in tmp_summary_dict.iteritems():
			self.integrationLoci_by_TST_uniqueDetails[tst][summary_field] = tupledata
		
	
	def filterUniqueIntegrationSites(self, filter_iss, tuple_keysBedFormat = False, stored_keysBedFormat = False):
		"""
		Given a tuple of ISS as input (bed or not), filter the global data [integrationLoci_by_TST_wrtRed, integrationLoci_by_TST_uniqueDetails] of ISS (unique) by removing ISS reported in the file. Pass Through filter.
		Same filter will be done on self.integrationLoci_by_TST_wrtRed dictionary (of course using values! here)
		This function is useful when you have to remove collisions or contaminations or specific ISs by input data.
		"""
		if stored_keysBedFormat: # if INPUT data are in bed format (and of course also existing data)
			# cicla si tutti i loci da eliminare
			if tuple_keysBedFormat: # if file data are also in BED format (rare):
				for (chr, start, end ) in filter_iss:
					## da ogni voce del main dictionary
					#if (chr, int(start), int(end) ) in dictionary_iss.keys(): # se K e' nel dizionario, quindi esiste come unico IS, eliminalo
					#	del dictionary_iss[(chr, int(start), int(end) )]
					# da ogni singolo dizionario di TST
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						if (chr, int(start), int(end) ) in loci_dict.values(): # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
							for locus, reflocus in loci_dict.iteritems():
								if reflocus == (chr, int(start), int(end)):
									tst_k_toprune[tst].append(locus) # append values (keys) to prune here
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]
							
					# dal dizionario dei dettagli
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						if (chr, int(start), int(end)) in loci_dict.keys(): # se tra i valori del dizionario, ovvero proprio gli IS unici elimina l'elemento chiave da loci_dict
							tst_k_toprune[tst].append((chr, int(start), int(end))) # append values (keys) to prune here
							print "[AP]\tFound ISs to prune:", chr, start, " into", tst
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]

			else: # =>file data of filters have only chr, locus -> cambia riferimenti di filter_iss 
				for (chr, start ) in filter_iss:
					## da ogni voce del main dictionary
					#if (chr, int(start), int(end) ) in dictionary_iss.keys(): # se K e' nel dizionario, quindi esiste come unico IS, eliminalo
					#	del dictionary_iss[(chr, int(start), int(end) )]
					# da ogni singolo dizionario di TST
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						for (tst_chr, tst_start, tst_end), (ref_chr, ref_start, ref_end) in loci_dict.iteritems():
							if chr == ref_chr and int(start) == ref_start: # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
								tst_k_toprune[tst].append((tst_chr, tst_start, tst_end)) # append values (keys) to prune here
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]
							
					# dal dizionario dei dettagli
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						for (ref_chr, ref_start, ref_end), details in loci_dict.iteritems():
							if chr == ref_chr and int(start) == ref_start: # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
								tst_k_toprune[tst].append((ref_chr, ref_start, ref_end)) # append values (keys) to prune here
								print "[AP]\tFound ISs to prune:", chr, start, " into", tst
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]
							
		else:
			# cicla si tutti i loci da eliminare
			if tuple_keysBedFormat: # se solo i dati da prunare sono in BED 
				for (chr, start, end ) in filter_iss:
					## da ogni voce del main dictionary
					#if (chr, int(start), int(end) ) in dictionary_iss.keys(): # se K e' nel dizionario, quindi esiste come unico IS, eliminalo
					#	del dictionary_iss[(chr, int(start), int(end) )]
					# da ogni singolo dizionario di TST
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						if (chr, int(start) ) in loci_dict.values(): # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
							for locus, reflocus in loci_dict.iteritems():
								if reflocus == (chr, int(start)):
									tst_k_toprune[tst].append(locus) # append values (keys) to prune here
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]
							
					# dal dizionario dei dettagli
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						if (chr, int(start)) in loci_dict.keys(): # se tra i valori del dizionario, ovvero proprio gli IS unici elimina l'elemento chiave da loci_dict
							tst_k_toprune[tst].append((chr, int(start))) # append values (keys) to prune here
							print "[AP]\tFound ISs to prune:", chr, start, " into", tst
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]
			else:
				for (chr, start ) in filter_iss:
					## da ogni voce del main dictionary
					#if (chr, int(locus) ) in dictionary_iss.keys(): # se K e' nel dizionario, quindi esiste come unico IS, eliminalo
					#	del dictionary_iss[(chr, int(locus) )]
					# da ogni singolo dizionario di TST
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						if (chr, int(start) ) in loci_dict.values(): # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
							for locus, reflocus in loci_dict.iteritems():
								if reflocus == (chr, int(start)):
									tst_k_toprune[tst].append(locus) # append values (keys) to prune here
					for tst, loci_dict in self.integrationLoci_by_TST_wrtRed.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]
							
					# dal dizionario dei dettagli
					tst_k_toprune = {} # chiavi da prunare specifiche per questo dizionario
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems(): # cicla TUTTI i diversi dizionari di ogni TST, dove k = tst  e v = loci_dict
						tst_k_toprune[tst] = []
						if (chr, int(start) ) in loci_dict.keys(): # se tra i valori del dizionario, ovvero proprio gli IS unici elimina l'elemento chiave da loci_dict
							tst_k_toprune[tst].append((chr, int(start))) # append values (keys) to prune here
							print "[AP]\tFound ISs to prune:", chr, start, " into", tst
					for tst, loci_dict in self.integrationLoci_by_TST_uniqueDetails.iteritems():
						for k in tst_k_toprune[tst]:
							del loci_dict[k]
				
	def removeIsBySequenceCount(self, sc_threshold, group_ids = None):
		"""
		Input: the threshold, the group ids (as list of keys relative to integrationLoci_by_TST_uniqueDetails); default is None -> it means to run filter to all groups of integrationLoci_by_TST_uniqueDetails
		Output: updated integrationLoci_by_TST_wrtRed and integrationLoci_by_TST_uniqueDetails
		Logics: for each locus, evaluate sequence count: if this number is below threshold (strict: seq_count < threshold) then discard it.
		"""
		tmp_prune_dict = {}
		if group_ids is not None: # controlla se vuoi lanciare su tutti i gruppi o solo su alcuni (qui ramo degli eletti)
			for tst in group_ids:
				if tst in self.integrationLoci_by_TST_uniqueDetails.keys():
					# start updating tmp dictionary
					tmp_prune_dict[tst] = [] # array degli elementi da eliminare (tuple array!)
					for i, det in self.integrationLoci_by_TST_uniqueDetails[tst].iteritems():
						if det['sequence_count'] < sc_threshold:
							tmp_prune_dict[tst].append(i)
							#print "trovato questo locus da tagliare:", i, "con sc =", det['sequence_count']
				############# now delete items ###########
					# from UNIQUE
					for locus in tmp_prune_dict[tst]:
						del self.integrationLoci_by_TST_uniqueDetails[tst][locus]
						print "[AP]\tObject integrationLoci_by_TST_uniqueDetails updated (removed elements)"
					# and from UNIQUE wrt redundant 
						tst_k_toprune = [] # array of local ISs associated to redundant loci for each tst
						if tst not in self.integrationLoci_by_TST_wrtRed.keys(): # nel caso non trovi la chiave (che puo' essere se fai prima l'espansione), segnalalo
							print "[AP]\tWarning: from self.integrationLoci_by_TST_wrtRed[tst] I cannot update", tst, "key not found; seems right?"
						else:
							if locus in self.integrationLoci_by_TST_wrtRed[tst].values(): # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
								for locallocus, reflocus in self.integrationLoci_by_TST_wrtRed[tst].iteritems():
									if reflocus == locus:
										tst_k_toprune.append(locus) # append values (keys) to prune here
							for k in tst_k_toprune:
								del self.integrationLoci_by_TST_wrtRed[tst][k]
								print "[AP]\tObject integrationLoci_by_TST_wrtRed updated (removed elements)"
		else:
			for tst, loci in self.integrationLoci_by_TST_uniqueDetails.iteritems():
				# start updating tmp dictionary
				tmp_prune_dict[tst] = []
				for i, det in loci.iteritems():
					if det['sequence_count'] < sc_threshold:
						tmp_prune_dict[tst].append(i)
			############# now delete items ###########
			for tst, tupledata in tmp_prune_dict.iteritems():
				# from UNIQUE
				for locus in tupledata:
					del self.integrationLoci_by_TST_uniqueDetails[tst][locus]
					print "[AP]\tObject integrationLoci_by_TST_uniqueDetails updated (removed elements)"
				# and from UNIQUE wrt redundant 
					tst_k_toprune = [] # array of local ISs associated to redundant loci for each tst
					if tst not in self.integrationLoci_by_TST_wrtRed.keys(): # nel caso non trovi la chiave (che puo' essere se fai prima l'espansione), segnalalo 
						print "[AP]\tWarning: from self.integrationLoci_by_TST_wrtRed[tst] I cannot update", tst, "key not found; seems right?" 
					else:
						if locus in self.integrationLoci_by_TST_wrtRed[tst].values(): # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
							for locallocus, reflocus in self.integrationLoci_by_TST_wrtRed[tst].iteritems():
								if reflocus == locus:
									tst_k_toprune.append(locus) # append values (keys) to prune here
						for k in tst_k_toprune:
							del self.integrationLoci_by_TST_wrtRed[tst][k]
							print "[AP]\tObject integrationLoci_by_TST_wrtRed updated (removed elements)"
	
	
	def removeIsByRatioOnSequenceCount(self, group_ids, ratio_threshold = 0.1):
		"""
		Input: the threshold, the group ids (as list of keys relative to integrationLoci_by_TST_uniqueDetails) to be used TOGETHER to fill the formula; default is None -> it will prompt to input data and select groups
		Output: updated integrationLoci_by_TST_wrtRed and integrationLoci_by_TST_uniqueDetails
		Logics: for each locus, evaluate sequence count: if this number is below ratio then discard it.
		Formula/Ratio: 
			[Excel based] =IF(AND(MAX($AT2:$AV2)>0,AT2>0),(IF(MAX(AU2,AV2)>0,(IF(AT2/MAX(AU2,AV2)<0.1,"",AT2)),AT2)),"")
			[Traslated]
				- ISs must have at least 1 value for all input groups
				- for each group, evaluate ratio: IS_sc_thisG / IS_max_sc_allOtherG < 0.1 [ratio_threshold]
					if True -> IS is accepted
					else discarded from THIS group [because for example is a contamination]
		"""
		loci_toremove = {} # per ogni gruppo di input, lista dei loci da eliminare da SOLO questo gruppo (perche' ratio sotto la soglia)
		if group_ids is not None and len(group_ids) > 1 : # controlla se vuoi lanciare su tutti i gruppi o solo su alcuni (qui ramo degli eletti)
			for group in group_ids:
				if group in self.integrationLoci_by_TST_uniqueDetails.keys(): # looppa tra gruppi
					loci_toremove[group] = []
					### primo controllo: tutti gli altri seq count > 0 (degli altri gruppi)
					other_group_ids = list(set(group_ids) - set([group])) # all other groups
					print "[AP]\tEvaluating", group, "versus", other_group_ids
					for locus, details in self.integrationLoci_by_TST_uniqueDetails[group].iteritems(): # per ogni locus, verifica che esista almeno negli altri
						locus_singleton_group = False
						other_loci_set = set() # set di tutti gli altri loci che ora riempio
						for otherg in other_group_ids: # cerca questo locus negli altri gruppi
							other_loci_set |= set(self.integrationLoci_by_TST_uniqueDetails[otherg].keys())
						if locus not in other_loci_set: # se non c'e' nel set dei loci globale (1 sola volta perche' set)
							locus_singleton_group = True
						if not locus_singleton_group: # solo se NON stai trattando un locus che e' unico di questo gruppo, allora calcola il ratio, else vai pure avanti e tienilo per questo lineage
							other_group_sc = [] # contiene i valori degli altri gruppi come sequence count
							for otherg in other_group_ids: # cerca questo locus negli altri gruppi
								if self.integrationLoci_by_TST_uniqueDetails[otherg].has_key(locus): # non essendo detto che c'e' per tutti i gruppi, controlla in quali e poi appendilo all'array
									other_group_sc.append( self.integrationLoci_by_TST_uniqueDetails[otherg][locus]['sequence_count'] ) # appendi valori sc nella var array
							ratio = details['sequence_count']/float(max(other_group_sc)) # compute ratio
							#print locus, ratio
					### secondo controllo: ratio > threshold
							if ratio < ratio_threshold: # se il rapporto e' inferiore alla soglia, eliminalo da QUESTO!! per farlo, prima memorizza valore e poi lo elimini (else runtime error su dizionario)
								loci_toremove[group].append(locus)
# 								if locus == ('1', 23446003, 23446004):
# 									print "ratio below threshold for is", locus, "adding locus to prune list"
			# ora esegui il tagli degli elementi
			# from UNIQUE
			for group, loci2prune in loci_toremove.iteritems():
				print "[AP]\tElements to remove for group %s due to ratio < threshold: %d" %(group, len(loci2prune))
				for locus in loci2prune:
					del self.integrationLoci_by_TST_uniqueDetails[group][locus]
# 					if locus == ('1', 23446003, 23446004) and locus not in self.integrationLoci_by_TST_uniqueDetails[group].keys():
# 						print "locus", locus, "removed from UNIQUE"
			# from UNIQUE wrt RED
					tst_k_toprune = [] # array of local ISs associated to redundant loci for each tst
					if group not in self.integrationLoci_by_TST_wrtRed.keys(): # nel caso non trovi la chiave (che puo' essere se fai prima l'espansione), segnalalo
						print "[AP]\t\tWarning: from self.integrationLoci_by_TST_wrtRed[group] I cannot update", group, "key not found; seems right?"
					else:
						if locus in self.integrationLoci_by_TST_wrtRed[group].values(): # se tra i valori del dizionario, ovvero proprio gli IS unici (le chiavi di dictionary_iss), trovi questo locus, allora accedi alla chiave di questo dizionario ed elimina l'elemento chiave da loci_dict
							for locallocus, reflocus in self.integrationLoci_by_TST_wrtRed[group].iteritems():
								if reflocus == locus:
									tst_k_toprune.append(locus) # append values (keys) to prune here
						for k in tst_k_toprune:
							del self.integrationLoci_by_TST_wrtRed[group][k]
							#print "[AP]\tObject integrationLoci_by_TST_wrtRed updated (removed elements)"
		else:
			print "[AP]\tYou MUST specify at least 2 elements to compare (no default values). Try again."
							
	
	def getSharedISs(self, group_ids):
		"""
		Input: groups ids to parse into integrationLoci_by_TST_uniqueDetails
		Output: dictionary of ISs with sequence count: k = ISs, v = dictionary:: k = groups, v = sequence count from this group
		Logics: for each group, extract only ISs shared, then collate data into tuple array
		"""
		shared_loci = {}
		shared_loci_set = set() # set of shared loci
		if group_ids is not None and len(group_ids) > 1 : # controlla se vuoi lanciare su tutti i gruppi o solo su alcuni (qui ramo degli eletti)
			loopindex = 0
			for group in group_ids:
				if loopindex == 0: # se primo ciclo, riempi set, else inizi a fare intersezione
					if group in self.integrationLoci_by_TST_uniqueDetails.keys(): # looppa tra gruppi
						shared_loci_set |= set(self.integrationLoci_by_TST_uniqueDetails[group].keys())
					else:
						print "[AP]\tWarning: this group does not exist into integrationLoci_by_TST_uniqueDetails", group
					loopindex += 1
				else: # loopindex > 0
					if group in self.integrationLoci_by_TST_uniqueDetails.keys(): # looppa tra gruppi
						shared_loci_set &= set(self.integrationLoci_by_TST_uniqueDetails[group].keys())
					else:
						print "[AP]\tWarning: this group does not exist into integrationLoci_by_TST_uniqueDetails", group
			# ora che hai tutti i loci intersezione stretta nel set, mettili in output
			for locus in shared_loci_set:
				shared_loci[locus] = {}
				for group in group_ids:
					if group in self.integrationLoci_by_TST_uniqueDetails.keys(): # looppa tra gruppi
						shared_loci[locus][group] = self.integrationLoci_by_TST_uniqueDetails[group][locus]['sequence_count']
			return shared_loci
		else:
			print "[AP]\tYou MUST specify at least 2 elements to compare (no default values). Try again."
			
	
	def groupElementsOfTST(self, grouping_dictionary):
		"""
		Input: dictionary of groups: k = new object to create; v = list of keys to group belonging to integrationLoci_by_TST_uniqueDetails
		Output: updates of integrationLoci_by_TST_uniqueDetails
		Logics: for each k of input dictionary, look at the values into integrationLoci_by_TST_uniqueDetails as keys and merge datasets
		"""
		test = ('1', 26767110, 26767111)
		tmp_outdict = {} # struttura dict temp che poi verra' aggiunta a UNIQUE
		for group_id, list_tst in grouping_dictionary.iteritems(): # loop sui singoli gruppi 
			print "[AP]\tGrouping elements of", group_id
			tmp_outdict[group_id] = {} # inizializza dizionario di gruppo
			iss_set = set() # set degli iss del gruppo (perche' posso avere iss comuni nell'unione e devo cambiare i valori dei sequence count!)
			for tst in list_tst: # per ogni tst della lista
				#print "[AP]\tAnayzing tst", tst, "for group", group_id
				if self.integrationLoci_by_TST_uniqueDetails.has_key(tst): # mero check su consistenza, else ERROR
					if len(iss_set) > 0:
						for locus, details in self.integrationLoci_by_TST_uniqueDetails[tst].iteritems():
# 							if locus == test:
# 									print locus, tst, "ha sc:", details['sequence_count']
							if locus in iss_set: # se questo iss e' gia' nel set corrente (quindi ho una sovrapposizione di iss tra insiemi tst), allor adevo aumentare il sequence count e il valore loci (merging data)							
								tmp_outdict[group_id][locus]['sequence_count'] += details['sequence_count'] # aggiungi sc
								tmp_outdict[group_id][locus]['loci'] |= details['loci'] # unisci loci
# 								if locus == test:
# 									print locus, "nel gruppo", group_id, "esiste gia' -> aggiorno sc tmp_outdict a:", tmp_outdict[group_id][locus]['sequence_count']
							else:
								tmp_outdict[group_id][locus] = details.copy() # copia direttamente questo locus
								iss_set.add(locus) # aggiorna set degli iss visitati!
# 								if locus == test:
# 									print locus, group_id, "nel gruppo e' NUOVO -> inserico sc attuale:", tmp_outdict[group_id][locus]['sequence_count']
					else:
						iss_set |= set(self.integrationLoci_by_TST_uniqueDetails[tst].keys()) # aggiungi questi iss al set
						tmp_outdict[group_id] = {}
						for locus, details in self.integrationLoci_by_TST_uniqueDetails[tst].iteritems():
							tmp_outdict[group_id][locus] = details.copy() # aggiungi TUTTO il dataset corrente al gruppo, as is
# 						if test in iss_set:
# 							print test, group_id, "\tnel gruppo ho appena inserito test, con sc:", tmp_outdict[group_id][test]['sequence_count']
				else:
					print "[AP]\tERROR: no such element found in integrationLoci_by_TST_uniqueDetails to be grouped! you asked to group it into the group id", group_id
		
		# ora aggiorna self.integrationLoci_by_TST_uniqueDetails con questi valori di gruppo
		for k, v in tmp_outdict.iteritems():
			if k in self.integrationLoci_by_TST_uniqueDetails.keys():
				print "[AP]\tATTENTION: this key already exist in the structure:", k, ". It will be updated!"
			print "[AP]\tUpdating integrationLoci_by_TST_uniqueDetails with new group:", k
			self.integrationLoci_by_TST_uniqueDetails[k] = v.copy()
		tmp_outdict = None # delete vars
		
	
	def getUniqueIntegrationLociByTissueSample(self, cursor, schema, table, experiment_loci_dict, patient, tissue, sample, treatment_list, keysBedFormat = False):
		"""
		IN: 
			experiment_loci_dict = dictionary obtained running getUniqueIntegrationLoci
		Logics: given experiment_loci_dict, extract for each tissue/sample the specific extracted loci by comparing query data with experiment_loci_dict keys.
		
		Look for all real IS data filtering by tissue and sample:
		- get sorted data -> sorted tuple
		- for each IS, count +- 3 rule and create a dictionary: given N integrations for a single IS, k = chr, locus; v = tuple of N integrations
		- self.integrationLoci_by_TS will be updated with this data
		"""
		# check input experiment_loci_dict
		if len(experiment_loci_dict.keys()) == 0:
			print "[AP]\tAcquire UNIQUE GLOBAL loci first, using getUniqueIntegrationLoci"
		else:
			treatment_in_clause = "', '".join( str(x) for x in treatment_list)
			print "[DB]\tComputing unique integration loci extraction (+/- 3 bp) from DB for", tissue, sample
			
			query = """
				SELECT DISTINCT `chr` , `integration_locus` , `sample` , `tissue` , `treatment` , strand, sum( sequence_count ) AS sequence_count, refSeqSymbol, orientation_ucsc 
				FROM  `%(schema)s`.`%(table)s`
				WHERE `chr` not like '0'
					AND tissue = '%(tissue)s'
					AND sample = '%(sample)s'
					AND treatment in ('%(treatment_in_clause)s')
				GROUP BY `chr` , `integration_locus` , `sample` , `tissue` , `treatment`
				ORDER BY  `%(table)s`.`chr` ASC,  `%(table)s`.`integration_locus` ASC 
				-- LIMIT 1000
				;
				""" %{
					  'select': '',
					  'schema': schema, 
					  'table': table,
					  'tissue': tissue,
					  'sample': sample,
					  'treatment_in_clause': treatment_in_clause,
					  }
			#print query
			cursor.execute(query)
			results = cursor.fetchall()
			
			# insert patient if not exists
			if not self.integrationLoci_by_TS.has_key(patient):
				self.integrationLoci_by_TS[patient] = {}
			
			# problema da risolvere: se uso solo i dati non ridondanti, mi perdo il possibile clustering dei siti +/-3, quindi mi serve espanderli qui e raccogliere le info ma mettendo come riferimento if IS reference (quello del dizionario, quello riportato nel file!!!!)
			if len(self.integrationLoci_by_TST_wrtRed.keys()) > 0: # nel caso in cui tu abbia lanciato la forma improved: getUniqueIntegrationLociImproved, tu hai gia' tutto!!!!! qui aggiungi proprio i dettagli, come seq count e geni, ecc
				print "[AP]\n[AP]\tUse the variable that you have!!!  self.integrationLoci_by_TST_wrtRed\n[AP]\n"
			else: # devi fare una pessima ricerca completa su dizionario.... molto inefficiente
				integration_locus = None
				for res in results:
					# insert data into global structure, based on tissue and sample
					integration_locus = int(res['integration_locus'])
					if keysBedFormat: # if format is BED like of keys
						for (chr, start, end), seq_array in experiment_loci_dict.iteritems():
							for iss_datails in seq_array:
								if iss_datails['chr'] == str(res["chr"]) and int(iss_datails['integration_locus']) == integration_locus: 
									# add values into output data structure
									if self.integrationLoci_by_TS[patient].has_key((tissue, sample)):
										# anche se non e' proprio il primo questo IS, mettilo tra i dati!!!! solo cosi' funza
										self.integrationLoci_by_TS[patient][(tissue, sample)][(chr, start, end)] = res['sequence_count']
									else:
										self.integrationLoci_by_TS[patient][(tissue, sample)] = {
												(chr, start, end): res['sequence_count']
											}
					else:
						for (chr, start, end), seq_array in experiment_loci_dict.iteritems():
							for iss_datails in seq_array:
								if iss_datails['chr'] == str(res["chr"]) and int(iss_datails['integration_locus']) == integration_locus: 
									if self.integrationLoci_by_TS[patient].has_key((tissue, sample)):
										# anche se non e' proprio il primo questo IS, mettilo tra i dati!!!! solo cosi' funza
										self.integrationLoci_by_TS[patient][(tissue, sample)][(chr, start)] = res['sequence_count']
									else:
										self.integrationLoci_by_TS[patient][(tissue, sample)] = {
												(chr, start): res['sequence_count']
											}		
								
#				## OLD, obsolete			
# 				if keysBedFormat: # if format is BED like of keys				
# 					if experiment_loci_dict.has_key( (str(res["chr"]), integration_locus, integration_locus + 1) ): # only if locus is in the global dictionary, take it!!
# 						if self.integrationLoci_by_TS[patient].has_key((tissue, sample)):
# 							# check if locus format must be BED like or not
# 							self.integrationLoci_by_TS[patient][(tissue, sample)][( str(res["chr"]), integration_locus, integration_locus + 1 )] = res
# 						else:
# 							self.integrationLoci_by_TS[patient][(tissue, sample)] = {
# 									( str(res["chr"]), integration_locus, integration_locus +1 ): res
# 								}		
# 				else:
# 					if experiment_loci_dict.has_key( (str(results[i]["chr"]), integration_locus) ):
# 						if self.integrationLoci_by_TS[patient].has_key((tissue, sample)):
# 							self.integrationLoci_by_TS[patient][(tissue, sample)][( str(res["chr"]),integration_locus )] = res
# 						else:
# 							self.integrationLoci_by_TS[patient][(tissue, sample)] = {
# 									( str(res["chr"]),integration_locus ): res
# 								}			

	def getLAMNamesFromDB(self, cursor, schema, table,):
		"""
		Look for all sample data in the joined fields N_LAM AND POOL of the no_redundant table.
		example of sql:
			SELECT DISTINCT `n_LAM`, `pool`, `tag`
			FROM  `MLD01`.`no_redundant_MLD01_ALL_NEW_NAMES` 
		"""
		#pass
		print "[DB]\tRetrieving details of LAMs (n_LAM, POOL and TAG) from DB"
		outdata = {} # k = tuple of values, v = string derived from tuple
		this_table = ""
		if table.startswith('no_'):
			this_table = '_'.join( str(x) for x in table.split('_')[1:] )
		else:
			this_table = table
		# per prendere i nomi UNIVOCI delle LAM devi passare per la REDUNDANT, invece che dalle NO:REDUNDANT (come era versione 1)
		query = """
			-- SELECT DISTINCT  `n_LAM`, `pool`, `tag`, `complete_name` as lam_header
			SELECT DISTINCT  `n_LAM` ,  `pool` ,  `tag` ,  `tissue` ,  `treatment` ,  `sample` , `enzyme`, `complete_name` as lam_header
			FROM  `%s`.`%s`;
			""" %(	schema, this_table )
		cursor.execute(query)
		results = cursor.fetchall()
		for r in results:
			outdata[(r["n_LAM"], r["pool"], r["tag"])] = r["lam_header"].strip()
		return outdata
	
	def getLAMDataFromDB(self, cursor, schema, table,):
		"""
		Look for all LAM data. NB: it is safe if used with REDUNDANT!!!! not working with no-redundant
		example of sql:
			SELECT DISTINCT  `n_LAM` ,  `pool` ,  `tag` ,  `tissue` ,  `treatment` ,  `sample` , `enzyme`, `complete_name` as lam_header
			FROM  `MLD01`.`no_redundant_MLD01_ALL_NEW_NAMES` 
		"""
		#pass
		print "[DB]\tRetrieving all important details of LAMs from DB (`n_LAM` ,  `pool` ,  `tag` ,  `tissue` ,  `treatment` ,  `sample` , `enzyme`, `complete_name`)"
		outdata = {} # k = tuple of values, v = dict of values from db
		query = """
			SELECT DISTINCT  `n_LAM` ,  `pool` ,  `tag` ,  `tissue` ,  `treatment` ,  `sample` , `enzyme`, `complete_name` as lam_header
			FROM  `%s`.`%s`;
			""" %(	schema, table )
		cursor.execute(query)
		results = cursor.fetchall()
		for r in results:
			outdata[(r["n_LAM"], r["pool"], r["tag"])] = {		'lam_header': r["lam_header"].strip(), 
										'tissue': r["tissue"].strip(), 
										'treatment': r["treatment"].strip(), 
										'sample': r["sample"].strip(), 
										'enzyme': r["enzyme"].strip(), 
										'n_LAM': r["n_LAM"].strip(), 
										'pool': r["pool"].strip(), 
										'tag': r["tag"].strip(), 
									}
		return outdata


	def mapLamToTCT(self, cursor, schema, table,):
		"""
		Map LAM to tissue, sample (cell line-type), timepoint.
		This function output 2 dictionaries (map): k = (tissue, cell type, tp), v = complete_name of LAM
		Exec query and map all names: each 'complete_name' is already mapped to TCT:
			SELECT DISTINCT  `n_LAM` ,  `pool` ,  `tag` ,  `tissue` ,  `treatment` ,  `sample` ,  `complete_name` 
			FROM  `MLD01`.`redundant_MLD01_ALL_NEW_NAMES_bis` 
			ORDER BY  `tissue` ,  `treatment` ,  `sample`
		
		"""
		print "[DB]\tMapping LAMs to Tissue/CellType/TimePoint -> creating groups from DB"
		outdata_tct_lam = {} # k = tuple of values, v = set of lams; string derived from tuple
		outdata_lam_tct = {} # k = lam id, v = tuple (t, c, t)
		this_table = ""
		if table.startswith('no_'):
			this_table = '_'.join( str(x) for x in table.split('_')[1:] )
		else:
			this_table = table
		# per prendere i nomi UNIVOCI delle LAM devi passare per la REDUNDANT, invece che dalle NO:REDUNDANT (come era versione 1)
		query = """
			SELECT DISTINCT  `n_LAM` ,  `pool` ,  `tag` ,  `tissue` ,  `treatment` ,  `sample` , `enzyme`, `complete_name` as lam_header
			FROM  `%s`.`%s`
			ORDER BY  `tissue` ,  `treatment` ,  `sample` ;
			""" %(	schema, this_table )
		cursor.execute(query)
		results = cursor.fetchall()
		for r in results:
			if outdata_tct_lam.has_key((r["tissue"], r["sample"], r["treatment"])):
				outdata_tct_lam[(r["tissue"], r["sample"], r["treatment"])].add( r["lam_header"].strip() )
			else:
				outdata_tct_lam[(r["tissue"], r["sample"], r["treatment"])] = set( [ r["lam_header"].strip() ] )
			outdata_lam_tct[ r["lam_header"].strip() ] = (r["tissue"], r["sample"], r["treatment"])
		return outdata_tct_lam, outdata_lam_tct
		
	def mapTCTpPoolTag(self, cursor, schema, table, filterBadLam = False, sequence_schema = None, fileOfCounts = None):
		"""
		Map LAM to tissue, sample (cell line-type), timepoint.
		This function output 2 dictionaries (map): k = (tissue, cell type, tp), v = dict:: 'pool' list(), 'tag' list()
		Exec query and map all names: each 'complete_name' is already mapped to TCT:
			SELECT DISTINCT  `n_LAM` ,  `pool` ,  `tag` ,  `tissue` ,  `treatment` ,  `sample` ,  `complete_name` 
			FROM  `MLD01`.`redundant_MLD01_ALL_NEW_NAMES_bis` 
			ORDER BY  `tissue` ,  `treatment` ,  `sample`
		
		"""
		print "[DB]\tMapping LAMs (pool, tag) to Tissue/CellType/TimePoint -> creating groups from DB"
		outdata_tct_pooltag = {} # k = tuple of values, v = set of lams; string derived from tuple
		this_table = ""
		if table.startswith('no_'):
			this_table = '_'.join( str(x) for x in table.split('_')[1:] )
		else:
			this_table = table
		# per prendere i nomi UNIVOCI delle LAM devi passare per la REDUNDANT, invece che dalle NO:REDUNDANT (come era versione 1)
		query = """
			SELECT DISTINCT `n_LAM` , `pool` , `tag` , `sample` , `vector` , `tissue` , `treatment` , `enzyme`
			FROM  `%s`.`%s`
			ORDER BY `sample` , `tissue` , `treatment` ;
			""" %(	schema, this_table )
		cursor.execute(query)
		results = cursor.fetchall()
		for r in results:
			# cicla sui dati per estrarre le LAM; se il filtro bedLam e' attivo allora controlla il sequence count della LAM superi il filtro e la relativa soglia
			usethislam = True # check variable; default TRUE!!!
			this_lam_seq_counts = {}
			if filterBadLam and fileOfCounts is not None:
				#getTableCountsFromFile(self, count_file, column_name, pool, tag)
				mapping_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'mapping_seq_counter', r['pool'], r['tag'])
				mapping_unique_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'mapping_unique_seq_counter', r['pool'], r['tag'])
				trimmed_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'trimmed_seq_counter', r['pool'], r['tag'])
				repeats_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'repeats_seq_counter', r['pool'], r['tag'])
				repeats_unique_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'repeats_unique_seq_counter', r['pool'], r['tag'])
				#excluded_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'excluded_seq_counter', r['pool'], r['tag'])
				raw_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'raw_seq_counter', r['pool'], r['tag'])
				# formula to understand if LAM is good or bad
				#print "\tConditions IF:", (mapping_seq_counter + repeats_seq_counter)/float(raw_seq_counter), trimmed_seq_counter/float(raw_seq_counter), mapping_seq_counter, repeats_seq_counter, trimmed_seq_counter, raw_seq_counter
				if (mapping_seq_counter)/float(trimmed_seq_counter) > 0.1:
					usethislam = True
				else:
					usethislam = False
			elif filterBadLam and sequence_schema is not None:
				mapping_seq_counter, mapping_unique_seq_counter = self.getTableCountsFromDB(cursor, sequence_schema, 'mapping', r['pool'], r['tag'])
				trimmed_seq_counter, trimmed_unique_seq_counter = self.getTableCountsFromDB(cursor, sequence_schema, 'trimmed', r['pool'], r['tag'])
				repeats_seq_counter, repeats_unique_seq_counter = self.getTableCountsFromDB(cursor, sequence_schema, 'repeats', r['pool'], r['tag'])
				#excluded_seq_counter, excluded_unique_seq_counter = self.getTableCountsFromDB(cursor, sequence_schema, 'excluded', r['pool'], r['tag'])
				raw_seq_counter, raw_unique_seq_counter = self.getTableCountsFromDB(cursor, sequence_schema, 'raw', r['pool'], r['tag'])
				# formula to understand if LAM is good or bad
				if (mapping_seq_counter)/float(trimmed_seq_counter) > 0.1:
					usethislam = True
				else:
					usethislam = False
			else:
				print "[AP]\tNo quality filters set for LAMs: all LAMs are going to be used."
			# if this LAM passed quality filter
			if usethislam:
				if outdata_tct_pooltag.has_key((r["tissue"], r["sample"], r["treatment"])):
					outdata_tct_pooltag[(r["tissue"], r["sample"], r["treatment"])]['pool'].append( r["pool"].strip() )
					outdata_tct_pooltag[(r["tissue"], r["sample"], r["treatment"])]['tag'].append( r["tag"].strip() )
				else:
					outdata_tct_pooltag[(r["tissue"], r["sample"], r["treatment"])] = { 'pool': [ r["pool"].strip() ], 'tag': [ r["tag"].strip() ], }
			else:
				print "[AP]\tDiscarded LAM for quality reasons: pool %s, tag %s" %( r['pool'], r['tag'] )
			
		return outdata_tct_pooltag
	
	def addCountsToMapTSTDictionaryFromFile(self, tst_dictionary, fileOfCounts):
		"""
		Given the dictionary of counts, add the field 'couts' as dictionary of counts that you extract from LAM details file (here as input to method).
		"""
		if fileOfCounts is not None and os.path.isfile(fileOfCounts):
			for (tissue, sample, timepoint), details in tst_dictionary.iteritems():
				this_lam_seq_counts = {}
			
				#getTableCountsFromFile(self, count_file, column_name, pool, tag)
				mapping_seq_counter = self.getTSTTableCountsFromFile(fileOfCounts, 'mapping_seq_counter', tissue, sample, timepoint)
				mapping_unique_seq_counter = self.getTSTTableCountsFromFile(fileOfCounts, 'mapping_unique_seq_counter', tissue, sample, timepoint)
				trimmed_seq_counter = self.getTSTTableCountsFromFile(fileOfCounts, 'trimmed_seq_counter', tissue, sample, timepoint)
				repeats_seq_counter = self.getTSTTableCountsFromFile(fileOfCounts, 'repeats_seq_counter', tissue, sample, timepoint)
				repeats_unique_seq_counter = self.getTSTTableCountsFromFile(fileOfCounts, 'repeats_unique_seq_counter', tissue, sample, timepoint)
				#excluded_seq_counter = self.getTableCountsFromFile(fileOfCounts, 'excluded_seq_counter', tissue, sample, timepoint)
				raw_seq_counter = self.getTSTTableCountsFromFile(fileOfCounts, 'raw_seq_counter', tissue, sample, timepoint)
				# add counts to tmp dict that will fill out data
				if this_lam_seq_counts.has_key( (tissue, sample, timepoint) ):
					this_lam_seq_counts[(tissue, sample, timepoint)]['mapping_seq_counter'] += mapping_seq_counter
					this_lam_seq_counts[(tissue, sample, timepoint)]['mapping_unique_seq_counter'] += mapping_unique_seq_counter
					this_lam_seq_counts[(tissue, sample, timepoint)]['repeats_seq_counter'] += repeats_seq_counter
					this_lam_seq_counts[(tissue, sample, timepoint)]['repeats_unique_seq_counter'] += repeats_unique_seq_counter
				else:
					this_lam_seq_counts[(tissue, sample, timepoint)] = {
							'mapping_seq_counter': mapping_seq_counter,
							'mapping_unique_seq_counter': mapping_unique_seq_counter,
							'repeats_seq_counter': repeats_seq_counter,
							'repeats_unique_seq_counter': repeats_unique_seq_counter
						}
				# update input data structure
				details['counts'] = this_lam_seq_counts[(tissue, sample, timepoint)]
		else:
			print "[AP]\tNo valid file in input for counts in the phase of editing data."
	
	def getTableCountsFromDB(self, cursor, schema, tablename, pool, tag):
		"""
		With this method you get the number of sequences (redundant and unique) for each of the tables in sequejce_* set.
		This function is used in LAM details.
		"""
		summary_query_template = Template( """
			SELECT count(*) as seq_counter
			FROM  `%(schema)s`.`$table`
			WHERE `pool` = '$pool'
			  AND `tag` = '$tag' ;
			  """ %{
				'schema': schema, 
				'table': tablename,
				}
			)
		uniquereads_summary_query_template = Template( """
			SELECT COUNT( DISTINCT  `read` ) as unique_reads 
			FROM  `%(schema)s`.`$table`
			WHERE `pool` = '$pool'
			  AND `tag` = '$tag' ;
			  """ %{
				'schema': schema, 
				'table': tablename,
				}
			)
		table_summary_query = summary_query_template.substitute(
			table = tablename,
			pool = pool, 
			tag = tag,
			)
		#print table_summary_query
		cursor.execute(table_summary_query)
		table_query_results = cursor.fetchall()
		table_seq_counter = 0
		for s in table_query_results:
			table_seq_counter += int(s['seq_counter'])
		
		table_uniquereads_summary_query = uniquereads_summary_query_template.substitute(
			table = tablename,
			pool = pool, 
			tag = tag,
			)
		#print table_uniquereads_summary_query
		cursor.execute(table_uniquereads_summary_query)
		table_uniquereads_summary_results = cursor.fetchall()
		table_unique_seq_counter = 0
		for s in table_uniquereads_summary_results:
			table_unique_seq_counter += int(s['unique_reads'])
		print "[AP]\t\tPool %s, Tag %s\t%s numbers: #reads = %s, #unique_reads = %s" %(pool, tag, tablename, table_seq_counter, table_unique_seq_counter)
		
		# check consistency!!
		if table_seq_counter < table_unique_seq_counter:
			print "\n\n[AP]\t\tWARNING!! #unique_reads is greater than #seqs... THIS MUST BE IMPOSSIBLE!\n\n"
			sys.exit()
		return table_seq_counter, table_unique_seq_counter
		
	def getTableCountsFromFile(self, count_file, column_name, pool, tag):
		"""
		With this method you get the number of sequences (redundant and unique) for each of the tables that you RUN IN ADVANCE; so we are speeding up access!
		This function is used in LAM details.
		count_file must be CSV with \t separators
		column_name = column from which to extract reads number
		Output: only 1 element!!!! the 
		"""
		def getColumnPositionInHeader(row_list, column_name):
			#print row_list
			for k in range(0, len(row_list)):
				if row_list[k] == column_name:
					return k
		
		row_index = 0
		with open(count_file, 'rb') as inf:
			reader = csv.reader(inf, delimiter = "\t")
			try:
				for row in reader:
					if row_index == 0:
						counton_column_position = getColumnPositionInHeader(row, column_name)
						tag_colindex = getColumnPositionInHeader(row, 'tag')
						pool_colindex = getColumnPositionInHeader(row, 'pool')
					# find sequence counts for this column_name
					if pool == row[pool_colindex] and tag == row[tag_colindex]:
						return float(row[counton_column_position])
					row_index += 1
			except IOError:
				print '[AP]\tERROR: Cannot open', count_file
	
	def getTSTTableCountsFromFile(self, count_file, column_name, tissue, sample, treatment):
		"""
		With this method you get the number of sequences (redundant and unique) for each of the tables that you RUN IN ADVANCE; so we are speeding up access!
		This function is used in LAM details.
		count_file must be CSV with \t separators
		column_name = column from which to extract reads number
		Output: only 1 element!!!! the 
		"""
		def getColumnPositionInHeader(row_list, column_name):
			#print row_list
			for k in range(0, len(row_list)):
				if row_list[k] == column_name:
					return k
		
		row_index = 0
		count = 0
		with open(count_file, 'rb') as inf:
			reader = csv.reader(inf, delimiter = "\t")
			try:
				for row in reader:
					if row_index == 0:
						counton_column_position = getColumnPositionInHeader(row, column_name)
						tissue_colindex = getColumnPositionInHeader(row, 'tissue')
						sample_colindex = getColumnPositionInHeader(row, 'sample')
						treatment_colindex = getColumnPositionInHeader(row, 'treatment')
					# find sequence counts for this column_name
					if tissue == row[tissue_colindex] and sample == row[sample_colindex] and str(treatment) == row[treatment_colindex]:
						count += float(row[counton_column_position])
					row_index += 1
			except IOError:
				print '[AP]\tERROR: Cannot open', count_file
		return count
	
	
	def getGeneSequenceCount(self, cursor, schema, table, patient, tissue, sample, treatment):
		"""
		Get genes with sequence count for a SPECIFIC tissue, sample and treatment.
		Output: update of the structure self.sequencecount_by_genes -> # main data structure: dictionary
			k = patient, v = dict:: k = (tissue, sample, treatment), v = dict:: k = gene name refSeqSymbol, 'allgenes'; v = dict:: k = {seqcount, treatment, tissue, chr, orientation, lenght, agv_position }
			
			NOTE: allgenes has this tuple order:: refSeqSymbol, seqreads, treatment, chr, orientation, lenght, avg_position
			
		Query example to extract genes:
		
			SELECT refSeqSymbol AS gene, sum( sequence_count ) AS seqreads
			FROM `redundant_WAS1001_FREEZE_2012_newTags`
			WHERE sample = 'CD34'
			    AND treatment = '0'
			    AND tissue = 'BM'
			GROUP BY refSeqSymbol
			ORDER BY seqreads DESC
			
		   ** alternative considering integration sites:
		    
			SELECT refSeqSymbol AS gene, sum( sequence_count ) AS seqreads, `integration_locus`
			FROM `redundant_WAS1001_FREEZE_2012_newTags`
			WHERE sample = 'CD34'
			    AND treatment = '0'
			    AND tissue = 'BM'
			GROUP BY integration_locus
			ORDER BY `redundant_WAS1001_FREEZE_2012_newTags`.`refSeqSymbol` ASC

		"""
		query_genes_template = Template("""
			SELECT refSeqSymbol, sum( sequence_count ) AS seqreads, treatment, chr, strand as orientation, gene_lenght_refseq as lenght, AVG(integration_locus) as avg_position
			FROM `$schema`.`$table`
			WHERE `sample` = '$sample' AND `tissue` = '$tissue' AND `treatment` = '$treatment'
			GROUP BY refSeqSymbol
			ORDER BY seqreads DESC
		      """)
		query_genes = query_genes_template.substitute(
		    schema = schema,
		    table = table,
		    sample = sample,
		    tissue = tissue,
		    treatment = treatment,
		  )
		print "[DB]\tLooking for genes and sequence counts of", patient, sample, tissue, treatment
		cursor.execute(query_genes)
		results = cursor.fetchall()
		tuple_array = [] # bulk data structure:: refSeqSymbol, seqreads, treatment, chr, orientation, lenght, avg_position
		for r in results:
			if r['refSeqSymbol'] != "No gene":
				if self.sequencecount_by_genes.has_key(patient):
					if self.sequencecount_by_genes[patient].has_key( (tissue, sample, treatment) ):
						self.sequencecount_by_genes[patient][ (tissue, sample, treatment) ][r['refSeqSymbol']] = {
						    'seqcount': int(r['seqreads']), 
						    'treatment': r['treatment'], 
						    'tissue': tissue,
						    'chr': r['chr'], 
						    'orientation': r['orientation'], 
						    'lenght': r['lenght'],
						    'avg_position': int(r['avg_position']),
						  }
					else:
						print "\t\tAdd (tissue, sample, treatment) to DB genes"
						self.sequencecount_by_genes[patient][ (tissue, sample, treatment) ] = {
						    r['refSeqSymbol']: {
							'seqcount': int(r['seqreads']), 
							'treatment': r['treatment'], 
							'tissue': tissue,
							'chr': r['chr'], 
							'orientation': r['orientation'], 
							'lenght': r['lenght'],
							'avg_position': int(r['avg_position']),
							}
						    }
				else:
					print "\t\tAdd patient to DB genes"
					self.sequencecount_by_genes[patient] =  {
					      (tissue, sample, treatment): {
						  r['refSeqSymbol']: {
						      'seqcount': int(r['seqreads']), 
						      'treatment': r['treatment'], 
						      'tissue': tissue,
						      'chr': r['chr'], 
						      'orientation': r['orientation'], 
						      'lenght': r['lenght'],
						      'avg_position': int(r['avg_position']),
						    }
						 }
					      }
					      
					      
				tuple_row = r['refSeqSymbol'], int(r['seqreads']), r['treatment'], r['chr'], r['orientation'], r['lenght'], int(r['avg_position'])
				tuple_array.append( tuple_row )
		# store all data into a matrix (tuple) sortable
		self.sequencecount_by_genes[patient][ (tissue, sample, treatment) ]['allgenes'] = tuple_array 
	
	#def convertToBedFormat(self, tupleIss, extend_endposition):
		#"""
		#Given an input dictionary of ISs, returns a bed file of ISs. If extend_endposition is activated, it returns a bed with end position = start position + 3.
		#Per ogni elemento della lista, convertilo in formato BED. La lista deve avere già la disposizione BED: chr, start position, notes
		#"""
		#for i in 
		#pass
	
	def getGeneSequenceCountForSingleIS(self, cursor, schema, table, patient, tissue, sample, treatment):
		"""
		Get genes with sequence count BUT distinct for specific INTEGRATION SITES for a SPECIFIC tissue, sample and treatment. This is the difference between getGeneSequenceCount and getGeneSequenceCountForSingleIS. This meas that here I also count unique integration sites, in the SAME logics of the existing function getUniqueIntegrationLoci. The zscore is also computed.
		
		Output: update of the structure self.sequencecount_by_iss_genes -> # main data structure: dictionary
			k = patient, v = dict:: k = (tissue, sample, treatment), v = dict:: k = (gene name refSeqSymbol, chr, position); v = dict:: k = {seqcount, treatment, tissue, chr, orientation, lenght, agv_position, zscore_seqcount }
			
			NOTE: allgenes has this tuple order:: refSeqSymbol, seqreads, treatment, chr, orientation, lenght, avg_position, zscore_seqcount, percentage_sequencecount
			
		Query example to extract genes:
		
			SELECT refSeqSymbol AS gene, sum( sequence_count ) AS seqreads, `integration_locus`
			FROM `redundant_WAS1001_FREEZE_2012_newTags`
			WHERE sample = 'CD34'
			    AND treatment = '0'
			    AND tissue = 'BM'
			GROUP BY integration_locus
			ORDER BY `redundant_WAS1001_FREEZE_2012_newTags`.`refSeqSymbol` ASC

		"""
		query_genes_template = Template("""
			SELECT refSeqSymbol, sum( sequence_count ) AS seqreads, treatment, chr, strand as orientation, gene_lenght_refseq as lenght, integration_locus, integration_locus as avg_position
			FROM `$schema`.`$table`
			WHERE `sample` = '$sample' 
			  AND `tissue` = '$tissue' 
			  AND `treatment` = '$treatment' 
			  AND `chr` not like '0' 
			  AND refSeqSymbol not like 'No gene'
			GROUP BY integration_locus
			ORDER BY `chr` ASC,  `integration_locus` ASC 
		      """)
		query_genes = query_genes_template.substitute(
		    schema = schema,
		    table = table,
		    sample = sample,
		    tissue = tissue,
		    treatment = treatment,
		  )
		print "[DB]\tLooking for genes and sequence counts of", patient, sample, tissue, treatment
		cursor.execute(query_genes)
		results = cursor.fetchall()
		tuple_array = [] # bulk data structure:: refSeqSymbol, seqreads, treatment, chr, orientation, lenght, position, zscore_seqcount, percentage_sequencecount
		tuple_array_bed = [] # bulk data structure for bed file handlers:: chr, start, end, name, seq count, orientation 
		# first of all, compute UNIQUE integration loci identification
		dict_is = {} # k = (chr, position, gene); v = sequence count (global, thus merging all not unique ISs (+- 3 rule)
		unique_loci = []
		integration_locus = None
		for i in range(0,len(results)):
			# se la diff tra integration loci e' <= 3 E l'attuale locus in analisi ha una differenza <= 3 con il primo allora, inseriscilo nel dizionario
			if i > 0 and integration_locus is not None and ( abs(results[i]["integration_locus"] - results[i-1]["integration_locus"]) <= 3 and abs(integration_locus - results[i]["integration_locus"]) <= 3 ):
				dict_is[( str(results[i]["chr"]), int(integration_locus), results[i]["refSeqSymbol"] )] += int(results[i]['seqreads'])
			else:
				integration_locus = int(results[i]["integration_locus"])
				dict_is[( str(results[i]["chr"]), int(integration_locus), results[i]["refSeqSymbol"] )] = int(results[i]['seqreads'])
		
		# init vars to compute z score
		zscore_list = [[],[]] # simple array of size = len(dict_is.keys()), made by array of array: [[genes], [sequence count], [zscore - only after processing seqeunce count it is appended]]
		stata = StatAnalysis()
		# process zscore and append values
		for i, sc in dict_is.iteritems():
			zscore_list[0].append(i)
			zscore_list[1].append(sc)
		zscore_outlist = stata.z_score(zscore_list[1])
		zscore_list.append(zscore_outlist)
		# create a new dictionary to store comparable zscores
		dict_is_zscores = {}
		for i in range(0, len(zscore_list[0])):
			dict_is_zscores[zscore_list[0][i]] = zscore_list[2][i]
		if len(set(dict_is_zscores.keys()) - set(dict_is.keys())) is not 0: # check at least len consistency in keys
			print "[AP]\tError on dictionary keys for Zscore: data not reliable!!"
			sys.exit()
		
		# init vars to compute percentage score
		percscore_list = [[],[]] # simple array of size = len(dict_is.keys()), made by array of array: [[genes], [sequence count], [perc - only after processing seqeunce count it is appended]]
		# process percentage and append values
		for i, sc in dict_is.iteritems():
			percscore_list[0].append(i)
			percscore_list[1].append(sc)
		percscore_outlist = stata.percentage_score(percscore_list[1])
		percscore_list.append(percscore_outlist)
		# create a new dictionary to store comparable zscores
		dict_is_percscores = {}
		for i in range(0, len(percscore_list[0])):
			dict_is_percscores[percscore_list[0][i]] = percscore_list[2][i]
		if len(set(dict_is_percscores.keys()) - set(dict_is.keys())) is not 0: # check at least len consistency in keys
			print "[AP]\tError on dictionary keys for percentage score: data not reliable!!"
			sys.exit()
			
		for r in results:
			if ( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] ) in dict_is.keys(): # se ho questo sito di integrazione nella mia lista di siti unici, ovvero nel dizionario appena creato (da cui prendo SOLO sequence count, il resto lo posso prendere ancora da results)
				if self.sequencecount_by_iss_genes.has_key(patient):
					if self.sequencecount_by_iss_genes[patient].has_key( (tissue, sample, treatment) ):
						self.sequencecount_by_iss_genes[patient][ (tissue, sample, treatment) ][( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )] = {
						    'seqcount': dict_is[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )], # sequence count COMULATO!!!! per siti UNICI
						    'treatment': r['treatment'], 
						    'tissue': tissue,
						    'chr': r['chr'], 
						    'orientation': r['orientation'], 
						    'lenght': r['lenght'],
						    'avg_position': int(r['avg_position']),
						    'integration_locus': r['integration_locus'],
						    'zscore_seqcount': dict_is_zscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )],
						    'perc_seqcount': dict_is_percscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )],
						  }
					else:
						print "\t\tAdd (tissue, sample, treatment) to DB genes"
						self.sequencecount_by_iss_genes[patient][ (tissue, sample, treatment) ] = {
						    ( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] ): {
							'seqcount': dict_is[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )], # sequence count COMULATO!!!! per siti UNICI
							'treatment': r['treatment'], 
							'tissue': tissue,
							'chr': r['chr'], 
							'orientation': r['orientation'], 
							'lenght': r['lenght'],
							'avg_position': int(r['avg_position']),
							'integration_locus': r['integration_locus'],
							'zscore_seqcount': dict_is_zscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )],
							'perc_seqcount': dict_is_percscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )],
							}
						    }
				else:
					print "\t\tAdd patient to DB genes"
					self.sequencecount_by_iss_genes[patient] =  {
					      (tissue, sample, treatment): {
						  ( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] ): {
						      'seqcount': dict_is[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )], # sequence count COMULATO!!!! per siti UNICI
						      'treatment': r['treatment'], 
						      'tissue': tissue,
						      'chr': r['chr'], 
						      'orientation': r['orientation'], 
						      'lenght': r['lenght'],
						      'avg_position': int(r['avg_position']),
						      'integration_locus': r['integration_locus'],
						      'zscore_seqcount': dict_is_zscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )],
						      'perc_seqcount': dict_is_percscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )],
						    }
						  }
					      }
				
				# now create sortable global structures
				# 'allgenes' field -> for counts
				tuple_row = r['refSeqSymbol'], dict_is[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )], r['treatment'], r['chr'], r['orientation'], r['lenght'], int(r['integration_locus']), dict_is_zscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )], dict_is_percscores[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )] # refSeqSymbol, seqreads, treatment, chr, orientation, lenght, avg_position, zscore, percentage
				tuple_array.append( tuple_row )
				# 'bedtuple' field -> for bed data management
				# if orientatino is - -> end will be in -, else (also if it does not exist) is +
				if str(r['orientation']) == '2':
					tuple_row_bed = r['chr'], int(r['integration_locus']), int(r['integration_locus']) - 1, tissue + " " + sample + " " + str(treatment) , dict_is[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )], r['orientation'] # chr, start, end, name, seq count, orientation 
				else:
					tuple_row_bed = r['chr'], int(r['integration_locus']), int(r['integration_locus']) + 1, tissue + " " + sample + " " + str(treatment) , dict_is[( str(r["chr"]), int(r['integration_locus']), r["refSeqSymbol"] )], r['orientation'] # chr, start, end, name, seq count, orientation 
				tuple_array_bed.append( tuple_row_bed )
		
		# store all data into a matrix (tuple) sortable
		self.sequencecount_by_iss_genes[patient][ (tissue, sample, treatment) ]['allgenes'] = tuple_array 
		self.sequencecount_by_iss_genes[patient][ (tissue, sample, treatment) ]['bedtuple'] = tuple_array_bed
	

class CSVFileManager:
	"""
	Acquire and handle files in CSV format.
	"""
	def __init__(self, infile):
		self.infile = infile
		self.csvdataset = set() # out set of data from CSV, 1 for each row
		self.csvdatatuple = [] # out tuple array of data from CSV, 1 for each row
		
	def acquireDataAsSet(self, delimiter = "\t"):
		"""
		Given the input file, acquire data as set.
		"""
		if self.infile is not None and os.path.isfile(self.infile):
			row_index = 0
			with open(self.infile, 'rb') as inf:
					reader = csv.reader(inf, delimiter = delimiter)
					try:
						for row in reader:
							self.csvdataset.add( tuple([x for x in row]) )
							row_index += 1
					except IOError:
						print '[AP]\tERROR: Cannot open', infile
		else:
			print "[AP]\tError: problems reading your CSV file. Is the path correct?"
	
	def acquireDataAsTuple(self, delimiter = "\t"):
		"""
		Given the input file, acquire data as tuple.
		"""
		if self.infile is not None and os.path.isfile(self.infile):
			row_index = 0
			with open(self.infile, 'rb') as inf:
					reader = csv.reader(inf, delimiter = delimiter)
					try:
						for row in reader:
							self.csvdatatuple.append( tuple([x for x in row]) )
							row_index += 1
					except IOError:
						print '[AP]\tERROR: Cannot open', infile
		else:
			print "[AP]\tError: problems reading your CSV file. Is the path correct?"
			
	def updateLocationsAsNumbers(self, ):
		"""
		Avoiding possible errors, in this class I also added a method to convert CSV data in form of (chr, locus) as string into location number.
		In case of BED files/tuple, USE BEDFILEMANAGER class!!!! Not this one!
		"""
		if self.csvdatatuple is not None and self.csvdatatuple != []:
			tmp_csvdatatuple = []
			for k in self.csvdatatuple:
				if len(k) == 2: # location is the last value; if it would be BED file, you MUST use BEDFileManager instead!!!! This is only for CSV data
					tmp_csvdatatuple.append( (k[0], int(k[1])) )
			print "[AP]\tUpdating self.csvdatatuple with numbers as locations"
			self.csvdatatuple = None
			self.csvdatatuple = tmp_csvdatatuple
		else:
			print "[AP]\tFirst acquire data!! Run acquireDataAsXYZ"

class JobMonitor:
	"""
	Andrea - JobMonitor Class for controlling performances of jobs
	"""
	def __init__(self, jobId, cursor, schema, table):
		"""
		Init class
		"""
		self.cursor = cursor
		self.table = table
		self.schema = schema
		self.jobId = jobId
		print "[AP]\tJob Monitoring System is now active. Tracking job", self.jobId
		
	def insert(self, jobStatus, jobOutfile, programName, commandline, patientId, notes):
		"""
		Using connection.cursor, insert data into JMS table
		"""
		query = """
			INSERT INTO `%s`.`%s` (
			  jobId ,
			  startedAt ,
			  jobStatus ,
			  jobOutfile_basename ,
			  jobOutfile_path,
			  programName,
			  commandline ,
			  server ,
			  patientId ,
			  notes,
			  systemUsername,
			  systemEnv
			  )
			VALUES (
			  '%s',  now(),  '%s', '%s', '%s', '%s', '%s',  '%s',  '%s',  '%s', '%s',  '%s'
			  );
		      """ %(self.schema,
			    self.table, 
			    self.jobId, 
			    #strftime("%Y-%m-%d %H:%M:%S", gmtime()),
			    jobStatus, 
			    os.path.basename(jobOutfile),
			    os.path.dirname(jobOutfile),
			    programName + " [%s]" %(os.environ['_']),
			    commandline, 
			    '\n'.join( str(s) for s in os.uname() ),
			    patientId, 
			    notes,
			    os.environ['USER'],
			    '\n'.join( str(k)+" -> "+str(v) for k, v in os.environ.iteritems() ),
			    )
		#print query
		self.run(query)
	
	def update(self, jobStatus):
		"""
		Update JMS
		"""
		query = """
			UPDATE  `%s`.`%s`
			SET  finishedAt = now(),
			  jobStatus =  '%s'
			WHERE jobId = '%s';
		  """ %(    self.schema,
			    self.table, 
			    # strftime("%Y-%m-%d %H:%M:%S", gmtime()),
			    jobStatus,
			    self.jobId, 
			)
		self.run(query)
	
	def close(self):
		"""
		Destroy all
		"""
		pass
	
	def createTable(self, ):
		"""
		Create table if not exists
		"""
		query = """
		      CREATE TABLE IF NOT EXISTS `%s`.`%s` (
			`jobId` bigint(15) NOT NULL COMMENT 'date +"\%Y\%m\%d\%H\%M\%S"',
			`startedAt` datetime NOT NULL,
			`finishedAt` datetime DEFAULT NULL,
			`jobStatus` enum('running','finished','aborted','queued') NOT NULL,
			`jobOutfile_basename` varchar(255) DEFAULT NULL,
			`jobOutfile_path` text,
			`programName` varchar(255) DEFAULT NULL,
			`commandline` text NOT NULL,
			`server` varchar(255) NOT NULL,
			`patientId` varchar(100) NOT NULL,
			`notes` text NOT NULL,
			`systemUsername` varchar(50) DEFAULT NULL,
			`systemEnv` text,
			PRIMARY KEY (`jobId`)
		      ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
		""" %(	self.schema,
			self.table,
			    )
		self.run(query)
	
	def truncateTable(self, ):
		query = """
		      TRUNCATE TABLE `%s`.`%s`;
		""" %(	self.schema,
			self.table,
			    )
		self.run(query)
		
	def run(self, query):
		try:
			#print query
			self.cursor.execute(""" %s """ %(query) )
		except MySQLdb.Error, e:
			print "DATABASE Error", e.args, e
			sys.exit (1)


class CIS:
	"""
	Class for Common Insertion Sites Analysis
	  - parses Manfred R out files
	  
	Flow example:
		c=CIS()
		data="/home/andrea/Dropbox/TIGET/Workbench/R/CIS/Manfred/MLD03_cis.nocollisions.R.out"
		c.parseManfredROut(data)
		c.adjustChrXYM()
		c.writeCISfileFromManfredResults("/home/andrea/Dropbox/TIGET/Workbench/R/CIS/Manfred/manfred1206.csv")
	
	Steps:
		1. acquire data from R file (out file) -> remember that X and Y MUST BE replaced by 23 and 24 BEFORE running R
		2. adjust data chr X and Y
		3. write data into file tsv -> then you will import file into DB

	"""
	def __init__(self, ):
		"""
		Init class CIS
		"""
		self.cis_set = set()
		self.cis_app_dictionary = {}
		self.tmp_chr_sections = {}
		self.cis_data = {} # k = chrom, v = dict with k= order, v = ISs (loci) as set [NB: ALL VALUES ARE STRING!!!]
	
	def parseManfredROut(self, ROutfile):
		"""
		Input:
		  - R output file of CIS analysis
		Logics:
		  - look for cisOrder string "Ordnung" "ORDER" 
		  - extract all points and chromosomes
		"""
		
		def getChrs(rfileout):
			chr_list = re.findall(r"\"Chromosom\" \"(\d+)\"", rfileout)
			return chr_list
		
		if os.path.isfile(ROutfile):
			rfile = open(ROutfile, 'r').read()
			chromosomes = getChrs(rfile) # array
			chromosomes.sort()
				
			
			
			chr_sections = {}
			chromosome = '0'
			string_chrinit = "[1] \"Chromosom\" \""
			string_before = "[1] \"Auflistung der IS in CIS nach jeweils h\\xf6chster Ordnung<=30\""
			string_after = "[1] \"Simulationserg. f\\xfcr Anz. IS in CIS: MW, SD, p-Werte (eins., zweis,.)\""
			string_order = "[1] \"Ordnung\" \""
			
			cis_start_row = 0
			cis_end_row = 0
			chr_number = 'NULL'
			index = 0
			#start_listing_iss = False
			for row in open(ROutfile, 'r'):
				if row.startswith(string_chrinit):
					chr_numbers = re.findall(r"\"Chromosom\" \"(\d+)\"", row)
					if chr_numbers[0] != '':
						chr_number = chr_numbers[0]
					chr_sections[chr_number] = {'cis_start_row': '', 'cis_end_row': '', 'cis_content': '', 'order_list': [], 'iss_list': []}
				if row.startswith(string_before):
					#print "string before"
					cis_start_row = index + 1
				if row.startswith(string_after):
					#print "string after"
					cis_end_row = index
					if cis_end_row - cis_start_row >= 0:
						chr_sections[chr_number]['cis_start_row'] = cis_start_row
						chr_sections[chr_number]['cis_end_row'] = cis_end_row
				index += 1
			
			rows_ROutfile = open(ROutfile, 'r').read().split('\n')
			for k in chr_sections.keys():
				chr_sections[k]['cis_content'] = '\n'.join( str(x) for x in rows_ROutfile[chr_sections[k]['cis_start_row'] : chr_sections[k]['cis_end_row']] ) # slice of the file
			
			
			### look for data within ICS data section
			string_order = "[1] \"Ordnung\" \""
			for k in chr_sections.keys():
				self.cis_data[k] = {}
				order = None
				for row in chr_sections[k]['cis_content'].split('\n'):
					if row.startswith(string_order):
						order_found = re.findall(r"\"Ordnung\" \"(\d+)\"", row)
						if order_found[0] != '':
							order = order_found[0]
							self.cis_data[k][order] = set()
					else:
						if order is not None:
							self.cis_data[k][order] |= set( re.findall(r"\s+(\d{3,30})", row) )
		else:
			print "[ERROR]\tR file not existing; check it please."
		
		      
		
	
	def adjustChrXYM(self, ):
		"""
		If needed, R program accepts only chr with numbers thus X == 23, Y == 24, M == 25.
		Now change it back.
		"""
		if '23' in self.cis_data.keys():
			self.cis_data["X"] = self.cis_data['23']
			del self.cis_data['23']
		if '24' in self.cis_data.keys():
			self.cis_data["Y"] = self.cis_data['24']
			del self.cis_data['24']
		if '25' in self.cis_data.keys():
			self.cis_data["M"] = self.cis_data['25']
			del self.cis_data['25']
	
	def writeCISfileFromManfredResults(self, outfile, fileformat = None):
		"""
		Goal:
		  - write a file of different file format 
		    1. csv (default): chr, order, locus
		    2. bed: BEd file format with details of order: chr ("chr{1..Y}", locus, locus + 1, Order, strand (default +, this is NOT true strand)
		Preconditions:
		  - self.cis_data not null -> run parseManfredROut before!
		"""
		if self.cis_data is not None and len( self.cis_data.keys() ) > 0:
			#file format dependent section
			if fileformat == 'bed':
				# write a rich bed file
				cis_csvfile = csv.writer(open(outfile, 'wb'), delimiter='\t')
				for chr, v in self.cis_data.iteritems():
					for order, iss_set in v.iteritems():
						for il in iss_set:
							cis_csvfile.writerow( ["chr"+chr, il, int(il)+1, "Order_"+order, "999", "+"] )
			else: # default TSV (CSV with delimiter \t) -> suitable for db data inport!
				cis_csvfile = csv.writer(open(outfile, 'wb'), delimiter='\t')
				for chr, v in self.cis_data.iteritems():
					for order, iss_set in v.iteritems():
						for il in iss_set:
							cis_csvfile.writerow( [chr, order, il] )
		else:
			print "[ERROR]\tYou must have the variable cis_data of this object filled."
	
	def intersectCISwithAnnotationFile(self, infile, outfile, order):
		"""
		Given the input file of annotation (created by clinical_merge_annotation program), intersect data using the first 2 FIXED columns: chr, locus
		"""
		pass

class SaturationAnalysis:
	"""
	The class SaturationAnalysis contains a set of operation to perform bootstrap apprpoach on sequence dataset.
	"""
	def __init__(self, ):
		"""
		Init class
		"""
		print "[AP]\tSaturation Analysis: object instantiated."
		self.bootstrap_inputsets = {} # k = bin id; v = dictionary:: k2 = set number; v2 = set of elements (headers)
		self.bootstrap_loci = {} # k = bin id; v = dictionary:: k2 = set number; v2 = dictionary:: k3 = loci::Tuple(chr, start, end); v3 = number of retrieved sequences
		self.bootstrap_overlapping_loci_withbed = {} # ONLY overlpping reads!!! same structure of bootstra_loci:::: k = bin id; v = dictionary:: k2 = set number; v2 = dictionary:: k3 = loci::Tuple(chr, start, end); v3 = number of retrieved sequences
		self.bootstrap_stats = {} # k = bin id; v = dictionary:: k2 = set number; v2 = dictionary:: k3s = 'Avg_Reads_PerLocus', 'CoveredLoci_onBinLoci', Avg_Reads_PerLocus
		self.bootstrap_statistics = {} # k = bin id; v = dictionary:: k2 = set number; v2 = dictionary:: k3 = target table (e.g. : mapping, raw, trimmed, ...)
		self.experiment_loci = set() # set([ (chr, start)* ])
		self.extended_experiment_loci = set() # set([ (chr, start)*, (chr, start +1 +2 +3)* ])
		self.maxthreads = 2 # threading: max number of threads:: if this number (user can edit it) become > max number of cpu, then resize it
		if self.maxthreads > multiprocessing.cpu_count():
			self.maxthreads = multiprocessing.cpu_count() - 1
		self.maxwhereelements = 3000 # ATTENTION!!!! this number highly affetcs efficiency, scalability and stability. Highier is BETTER (i.e. 3000 or 5000), but not too much!
		self.current_table = None
		self.current_schema = None
		
	def fetchData(self, cursor, schema, table, select = "*", where = '1'):
		"""
		Method to retrieve and fetch sequences as tuple
		Output: tuple of data
		"""
		query = """
			SELECT %s
			FROM  `%s`.`%s`
			WHERE %s ;
			""" %(	select, schema, table, where )
		#print "[DB]\tFetching data form DB, ", schema, table
		#print query
		#pass
		cursor.execute(query)
		results = cursor.fetchall()
		return results
	
	def fetchGeneralData(self, cursor, query):
		"""
		Method to retrieve and fetch sequences as tuple
		Output: tuple of data
		"""
		#print query
		#pass
		cursor.execute(query)
		results = cursor.fetchall()
		return results
		
	def fetchSummaryDataAndUpdateBootstrap(self, cursor, schema, table, select, where, N_setId, B_binId):
		"""
		Method to retrieve and fetch count and update value into bootstrap method
		Output: -
		"""
		query = """
			SELECT %s
			FROM  `%s`.`%s`
			WHERE %s ;
			""" %(	select, schema, table, where )
		#print "[DB]\tFetching data form DB, ", schema, table
		#print query
		#pass
		cursor.execute(query)
		results = cursor.fetchall()
		count = 0
		for c in results:
			count += c['count_summary']
		if N_setId in self.bootstrap_statistics[B_binId].keys():
			self.bootstrap_statistics[B_binId][n][table] = count
		else:
			self.bootstrap_statistics[B_binId][n] = { table: count }

	def getListOfHeaders(self, cursor, schema, table, where = '1'):
		"""
		Method to retrieve and get header ids only
		Output: list of headers
		"""
		query = """
			SELECT `header`
			FROM  `%s`.`%s`
			WHERE %s ;
			""" %(	schema, table, where )
		print "[DB]\tFetching headers into list form DB, ", schema, table
		cursor.execute(query)
		results = cursor.fetchall()
		outdata = [j['header'] for j in results ] # list of header ids
		return outdata
		
	def bootstrap(self, cursor, schema, table, N_sets_number, K_elements_number, B_binId, wherefilter = '1'):
		"""
		INPUT:
		  - N_set = intero
		  - K_elements = intero
		  - B_binId = identifier of the Bin. Bin is for example a percentage of samplesize. Example: given samplesize I of 500 reads, I want divide I into growing percentages, from 10 to 99 at step 5, obtaining these bins: 10, 15, 20, 25, ..., 90, 95, 99.
		LOGICA:
		  dato il set iniziale di dati, suddividi il campione in N sotto gruppi di K elementi. Per ogni n in N, fai conteggi delle reads (come lam details)
		"""
		# 1. get initial sample set (all raw data) and header list
		#raw_data = self.fetchData( cursor, schema, table, "*", '1' )
		list_header = self.getListOfHeaders(cursor, schema, table, wherefilter)
		# 2. divide sample set into N sets of K elements, assigning to it a B id
		samplesize = len(list_header)
		print "[AP]\tDividing sample set into N subsets."
		self.bootstrap_inputsets[B_binId] = {}
		for n in range(0, N_sets_number):
			self.bootstrap_inputsets[B_binId][n] = random.sample(list_header, K_elements_number)
	
	def getLociCounts(self, cursor, schema, table, B_binId, count_field_selection = "*"):
		"""
		INPUT:
		  - self.bootstrap_inputsets
		LOGICS:
		  - for the specific B_binId, get count summary
		  
		QUERY to retrieve positions:
		    SELECT DISTINCT `chr` , `integration_locus`, COUNT(sequence_count) AS locus_seq_count
		    FROM `redundant_FREEZE_TWELVE_MONTHS`
		    WHERE `header` IN ('header1', 'header2', ...)
		    GROUP BY `chr` , `integration_locus`
		    ORDER BY `chr` ASC, `integration_locus` ASC
		    
		   alternative query for BED file export:
		    SELECT DISTINCT `chr` , `integration_locus` as start , `integration_locus` + 1 as end, COUNT( sequence_count ) AS locus_seq_count
		    FROM `redundant_FREEZE_SIX_MONTHS`
		    WHERE chr >0
		    GROUP BY `chr` , `integration_locus`
		    ORDER BY locus_seq_count ASC , `chr` ASC , `integration_locus` ASC
		    
		"""
		def sighandler(num, frame):
			"""
			Function to manage signal terminal/interrupt from user
			"""
			global sigint
			sigint = True
		
		def doMyWork(N_setId):
			"""
			Combination of steps for the LOOP:
			1. combine header for while/select statement
			2. run query
			3. manage result -> update self.bootstrap_loci variable
			"""
			header_list = self.bootstrap_inputsets[B_binId][N_setId]
			select = " DISTINCT `chr` , `integration_locus`, COUNT(sequence_count) AS locus_seq_count " 
			print "[DB]\tB_binId =", B_binId, "; SetId Number =", N_setId, "; schema and table =", self.current_schema, self.current_table
			
			##---------------------- working fine but very slow ------------------------------------
			where = " `header` IN ( '"
			where += "', '".join(header_list)
			where += "' ) GROUP BY `chr` , `integration_locus` ORDER BY `chr` ASC, `integration_locus` ASC "
			schema = self.current_schema
			table = self.current_table
			results = self.fetchData(cursor, schema, table, select, where) # returns a single entry
			# init N_setId on self.bootstra_loci var
			if N_setId not in self.bootstrap_loci[B_binId].keys():
				self.bootstrap_loci[B_binId][N_setId] = { } 
			# update self.bootstrap_loci var 
			for k in results:
				#print k
				self.bootstrap_loci[B_binId][N_setId][( k["chr"].strip().split('chr')[-1], int(k["integration_locus"]) )] = int(k["locus_seq_count"])

			
		######### init data ######
		self.current_table = table
		self.current_schema = schema
		
		if B_binId not in self.bootstrap_loci.keys():
			self.bootstrap_loci[B_binId] = {} # insert bin id into output set
		print "[AP]\tGet summary of reads from (schema table): ", self.current_schema, self.current_table
		
		sigint = False # keybord interrupt, now is false
		signal.signal(signal.SIGINT, sighandler)
		
		if len(self.bootstrap_inputsets) > 0 and self.bootstrap_inputsets is not None:
			#for each bin B, compute query counts and update variables
			map(doMyWork, self.bootstrap_inputsets[B_binId].keys()) # -> is like: for x in self.bootstrap_inputsets[B_binId].keys(): doMyWork(k)
				
		else:
			print "[WARNING]\tBefore running summary you must run bootstrap method."
			
	
	def acquireExperimentLoci(self, locations_from_bedfile, extend = False):
		"""
		Given BED file (BedFileHandler object .locations), convert it into a set using only chr and start position as values for each locus.
		"""
		# bed locations have 3 values as index (chr, start, end) while db locations only 2 (chr, start) -> convert bed locations into a set cutting unuseful info by creating a new var:: set([ (chr, start)* ])
		if locations_from_bedfile is not None:
			# create the set of experiment loci from bed file
			print "\n[AP]\tAcquiring Experiment Loci from BedFIleHandler object."
			for (chr, start, end), vals in locations_from_bedfile.iteritems():
				self.experiment_loci.add( (chr, start) )
				# add values to extended set: create a set of expanded experiment locations set -> this is because of the +/-3 issue for unique locus count
				self.extended_experiment_loci.add( (chr, start) )
				if extend: ## only if you set extension, consider neighboring (+ side!!) positions
					self.extended_experiment_loci.add( (chr, start + 1) )
					self.extended_experiment_loci.add( (chr, start + 2) )
					self.extended_experiment_loci.add( (chr, start + 3) )
			print "[AP]\t\t|\n[AP]\t\t+--> Available loci:", len(self.experiment_loci), "; extended to", len(self.extended_experiment_loci)
		else:
			print "[WARNING]\tBEDFileObject.locations must be not empty."
			
	
	def getSaturationDetails(self, ):
		"""
		Given the locations_from_bedfile (BedFileHandler object .locations), intersect data from BINS (bootstrap on DB) with data from BED (final integration loci from experiment)
		"""
		
		def mean(numberList):
			if len(numberList) == 0:
			    return float(0)
		    
			floatNums = [float(x) for x in numberList]
			return sum(floatNums) / len(numberList)
		
		def ratioLen(listN, listD):
			if len(listD) == 0:
			    return float(0) # TODO: SHOULD BE NAN, NOT A NUMBER AS NOW IS 
		    
			return len(listN) / float(len(listD))
			
		if len(self.experiment_loci) > 0 and len(self.extended_experiment_loci) > 0:
			if self.bootstrap_loci is not None:
				### remember var type: 
				### self.bootstrap_loci = {} # k = bin id; v = dictionary:: k2 = set number; v2 = dictionary:: k3 = loci::Tuple(chr, start, end); v3 = number of retrieved sequences
				### self.bootstrap_stats = {} # k = bin id; v = dictionary:: k2 = set number; v2 = dictionary:: k3s = 'Avg_Reads_PerLocus', 'CoveredLoci_onBinLoci', Avg_Reads_PerCoveredLocus
				for B_binId, dic_nsets in self.bootstrap_loci.iteritems():
					if B_binId not in self.bootstrap_stats.keys():
						self.bootstrap_stats[B_binId] = {} # insert bin id into output set
						self.bootstrap_overlapping_loci_withbed[B_binId] = {}
					for N_setId, dict_loci in dic_nsets.iteritems():
						print "[AP]\t\tBinID", B_binId, "; SetID", N_setId
						#compute counts
						Avg_Reads_PerLocus = None
						Avg_Reads_PerCoveredLocus = None
						CoveredLoci_onBinLoci = 0
						present_lociset_seqreads = [] # array of only existing loci seq count in experiment_loci from bin
						self.bootstrap_overlapping_loci_withbed[B_binId][N_setId] = {} # update overlapping loci with bed file
						all_lociset_seqreads = [] # array of all seq counts from all loci of bin
						total_loci_persetId = len( dict_loci.keys() )
						lindex = 1
						for locus, seqcount in dict_loci.iteritems():
							all_lociset_seqreads.append( seqcount )
							#print "[AP]\t\t+--> locus (index %d) %s; seqcount %d" %(lindex, locus, seqcount)
							if locus in self.extended_experiment_loci:
								present_lociset_seqreads.append( seqcount )
								self.bootstrap_overlapping_loci_withbed[B_binId][N_setId][locus] = seqcount
								#print "[AP]\t\t|\t#!! locus in the experiment BED list"
							lindex += 1
						
						Avg_Reads_PerLocus = mean(all_lociset_seqreads)
						Avg_Reads_PerCoveredLocus = mean(present_lociset_seqreads)
						if len(all_lociset_seqreads)>0:
							CoveredLoci_onBinLoci = len(present_lociset_seqreads) / float(len(all_lociset_seqreads))
						
						# now update dictionary stats
						if N_setId not in self.bootstrap_stats[B_binId].keys():
							self.bootstrap_stats[B_binId][N_setId] = {
								'Avg_Reads_PerBinSetLocus': Avg_Reads_PerLocus, 
								'CoveredLoci_onBinSetLoci': CoveredLoci_onBinLoci, 
								'Avg_Reads_PerCoveredLocus': Avg_Reads_PerCoveredLocus,
								'BinSet_Retrieved_Loci_onRequestedReads': total_loci_persetId,
								'BinSet_Retrieved_BEDLoci': len(present_lociset_seqreads),
								'BinSet_len_all_lociset_seqreads': len(all_lociset_seqreads),
								'BinSet_RequestedReads_Kelems': len(self.bootstrap_inputsets[B_binId][N_setId]),
								'SumOfSeqCounts_onBinSetLoci': sum(all_lociset_seqreads),
								'Experiment_AllLociFromBED': len(self.extended_experiment_loci),
								'Ratio_IS_presentbed_on_availabe': ratioLen(present_lociset_seqreads, all_lociset_seqreads),
								'Ratio_IS_presentbed_on_allexperimentIss': ratioLen(present_lociset_seqreads, self.extended_experiment_loci),
								}
						print "[AP]\t\t=> Summary: ", '; '.join( str(x + ": " + str(y)) for x, y in self.bootstrap_stats[B_binId][N_setId].iteritems() )
						#print "[AP]\t\t=> all_lociset_seqreads =", all_lociset_seqreads
						
			else:
				print "[WARNING]\tGet Bootstrap loci first. Loci var is: ", self.bootstrap_loci, "."
		else:
			print "[WARNING]\tGet Experiment loci first."
			#self.acquireExperimentLoci(locations_from_bedfile)
	
	
	def getCountSummary(self, cursor, schema, table, B_binId, count_field_selection = "*"):
		"""
		INPUT:
		  - self.bootstrap_inputsets
		LOGICS:
		  - for the specific B_binId, get count summary
		NOTE:
		  - here I use THREADS and SEMAPHORE (bounded for pools) [http://stackoverflow.com/questions/936933/boundedsemaphore-hangs-in-threads-on-keyboardinterrupt]
		"""
		
		def sighandler(num, frame):
			"""
			Function to manage signal terminal/interrupt from user
			"""
			global sigint
			sigint = True
		
		def optimizeQueryByUnion(select, header_list, max_where_ids = 500):
			"""
			Given a long query with many ids in where condition in OR, split it and merge them as Union
			"""
			query_template = Template( """
				SELECT $select
				FROM  `%(schema)s`.`%(table)s`
				WHERE  $where
				  
				""" %{
				      'schema': self.current_schema, 
				      'table': self.current_table,
				      }
				)
			# divide where clauses into W sub statements
			query_list = []
			where_counter = 0
			startpos = 0
			endpos = None
			# the for loop will divide headers into W groups based on steps given by max_where_ids; remaining elements go in the last group with the next IF
			for group in range(0, len(header_list), max_where_ids):
				startpos = group
				endpos = startpos + max_where_ids
				#print "\tsets from %d to %d" %(startpos, endpos)
				query_list.append( query_template.substitute(select = select, where = " `header` IN ( '" + "', '".join(header_list[startpos:endpos]) + "' ) " ) ) 
			
			query = " UNION ALL ".join( str(uq) for uq in query_list)
			return query
		
		def doMyWork(N_setId):
			"""
			Combination of steps for the LOOP:
			1. combine header for while/select statement
			2. run query
			3. manage result
			"""
			header_list = self.bootstrap_inputsets[B_binId][N_setId]
			select = " count(%s) as count_summary " %(count_field_selection)
			print "[DB]\tB_binId =", B_binId, "; SetId Number =", N_setId, "; schema and table =", self.current_schema, self.current_table
			
			##---------------------- working best, alternative solution -----------------------------------
			query_forcount = optimizeQueryByUnion(select, header_list, self.maxwhereelements) # this number is the maximum number of where ids per select
			#print query_forcount
			results = self.fetchGeneralData(cursor, query_forcount)
			##---------------------- working fine but very slow ------------------------------------
			#where = " `header` IN ( '"
			#where += "', '".join(header_list)
			#where += "' ) "
			#schema = self.current_schema
			#table = self.current_table
			#results = self.fetchData(cursor, schema, table, select, where) # returns a single entry
			
			
			count = None
			count = 0
			# print results
			for c in results:
				count += c['count_summary']
			print "[DB]\t\t..retrieved reads =>", count
			if N_setId in self.bootstrap_statistics[B_binId].keys():
				self.bootstrap_statistics[B_binId][N_setId][table] = count
			else:
				self.bootstrap_statistics[B_binId][N_setId] = { table: count }
			
			##----------------------------- not working parallel with threads ---------------------
			#pool_semaphore.release() # release a thread space
			
		######### init data ######
		self.current_table = table
		self.current_schema = schema
		
		if B_binId not in self.bootstrap_statistics.keys():
			self.bootstrap_statistics[B_binId] = {} # insert bin id into output set
		print "[AP]\tGet summary of reads from (schema table): ", self.current_schema, self.current_table
		
		sigint = False # keybord interrupt, now is false
		signal.signal(signal.SIGINT, sighandler)
		
		if len(self.bootstrap_inputsets) > 0 and self.bootstrap_inputsets is not None:
			##----------------------------- working sequentially ---------------------
			##0. FOR LOOP
			#for n, hl in self.bootstrap_inputsets[B_binId].iteritems(): # n = set number id, hl = header list
				##map(doMyWork, n)
				#doMyWork( n ) ## THIS IS WORKING FINE
			
			##1. MAP -> working fine ;)
			map(doMyWork, self.bootstrap_inputsets[B_binId].keys())
			##----------------------------- not working parallel ---------------------
			
			##1. THREADS and SEMAPHORE (as queue)
			#pool_semaphore = threading.BoundedSemaphore(value = self.maxthreads)
			#for n, hl in self.bootstrap_inputsets[B_binId].iteritems(): # n = set number id, hl = header list
				#pool_semaphore.acquire() # controlling threading, acquire a new thread
				#thread = threading.Thread(target = doMyWork, args = ([n])) # start thread on specific function (NOTE: called function HAS NOT CALLBACK!!! -> I do not know how to use callback into threads [into multiprocessing is easy, see http://stackoverflow.com/questions/884650/python-spawn-parallel-child-processes-on-a-multi-processor-system-use-multipr])
				#thread.start()
				## remenber to release semaphore of pool in the function!
				
				##p = Process(target = doMyWork, args = (cursor, schema, table, B_binId, count_field_selection, n) )
				##p.start()
				##p.join()
			##2. PROCESS POOL
			#processes_pool = multiprocessing.Pool(processes = self.maxthreads)
			#n = self.bootstrap_inputsets[B_binId].keys()
			#n.sort()
			#print "[AP]\tParallelizing jobs:: Processing these set ids:", n
			#presult = processes_pool.apply_async(doMyWork, n)
			#print(presult.get(timeout=1))
			##print(processes_pool.map(doMyWork, n))
			
			##thisres = processes_pool.map_async(doMyWork, n) 
			#print "[AP]\tProcesses are now running"
			##thisres.wait() # Wait on the results
			##print "[AP]\tWaiting response"
				
		else:
			print "[WARNING]\tBefore running summary you must run bootstrap method."

class AnalyzeBEDData:
	"""
	Compare BED data class
	"""
	def __init__(self, ):
		"""
		Init class
		"""
		self.beddata = {} # dictionary of BED data: k = 'BED data ID'; v = tuple of values 
		self.beddata_out_sequencecount = {} # dictionary of BED data: k = 'BED data ID'; v = tuple of values - this var is edited during object run, noot as input (like self.beddata)
		#print "[AP]\tAnalysis of BED data handler object initialized."
	
	def acquire(self, bedId, tupleBedData):
		"""
		Add to self.beddata an element, where bedId is the key and value is a set of tupleBedData (all points, i.e. the keys of the bed file handler).
		"""
		if self.beddata.has_key(bedId):
			print "[AP]\tError: this key already exists; change bedId."
		else:
			self.beddata[bedId] = tupleBedData
			#print "[AP]\tData for this key acquired: ", bedId
	
	def getCommonElements(self, key_list, ):
		"""
		Given an input array of keys, returns intersection of locations as set of tuple.
		NB: values of the set for each key MUST be of the same type and shape (<- i.e. chr, start, end only).
		"""
		all_valid_points = True
		# check key existance
		for k in key_list:
			if not self.beddata.has_key(k):
				all_valid_points = False
		# once everything is ok, continue
		if all_valid_points:
			retset = set()
			for k in key_list:
				if len(retset) == 0:
					retset = set(self.beddata[k])
				else:
					retset &= set(self.beddata[k])
			return retset
		else:
			print "[AP]\tErrors in key list (give: %s): some of then have not been previously acquired, please run first self.acquire." %(key_list)
	
	def getOverlappingElements(self, key_list, strictOverlappingAmongAll = False, keepBothRegions = False, orientationBased = False):
		"""
		Given an input array of keys, returns only the overlapping elements by position as set of tuple.
		NB: values of the set for each key MUST be of the same type and shape (<- i.e. chr, start, end only).
		WARNING: with input list > 2 elements, the overlapping output is now a NON strict overlapping -> to be considered as overlapping region, a point must be shared at least between 2 input sets, not necessary among all input sets.
		"""
		all_valid_points = True
		# check key existance
		for k in key_list:
			if not self.beddata.has_key(k):
				all_valid_points = False
		# once everything is ok, continue
		if all_valid_points:
			if len(key_list) == 1:
				print "[AP]\tNot so useful counting overlapping regions using only 1 input set, right? You will get all elements of the single input set!"

			#controlla che start_1 sia compreso tra start_2 e end_2, o end _2 sia tra start_2 ed end_2 per ogni set a confronto
			regionset = set()
			outregionset = set()
			for k in key_list: # per ogni lista BED da confrontare
				if len(regionset) == 0: # se e' il primo cliclo, inserisci tutto
					regionset = set(self.beddata[k])
				else: # altrimenti inizia a fare confronti sulla lunghezza
					for region_2 in self.beddata[k]:
						chr_2 = region_2[0]
						start_2 = region_2[1]
						end_2 = region_2[2]
						for region_1 in regionset:
							chr_1 = region_1[0]
							start_1 = region_1[1]
							end_1 = region_1[2]
							# compare and decide
							if chr_1 == chr_2: # chr first of all
								if orientationBased: # if you have to consider strands
									if (start_2 <= start_1 <= end_2 or end_2 <= start_1 <= start_2) or (start_2 <= end_1 <= end_2 or end_2 <= end_1 <= start_2):
										# report me regions; so far only 1, else both
										if keepBothRegions:
											outregionset.add(region_2)
											outregionset.add(region_1)
										else:
											outregionset.add(region_1)
								else: # in case of considering only single strand positions
									if (start_2 <= start_1 <= end_2) or (start_2 <= end_1 <= end_2):
										# report me regions; so far only 1, else both
										if keepBothRegions:
											outregionset.add(region_2)
											outregionset.add(region_1)
										else:
											outregionset.add(region_1)
					if not strictOverlappingAmongAll: # caso in cui NON richieda sovrapposizione stretta (ovvero: tutti gli input sets condividono/si sovrappongono quella regione), procedi qui unendo a regionset questi dati; se invece si richiede la sovrapposizione STRETTA, non aggiungere nulla a regionset
						# add all elements of this set to the comparison set
						regionset |= set(self.beddata[k]) # this is the point for STRICT analysis!!!
			return outregionset
		else:
			print "[AP]\tErrors in key list (give: %s): some of then have not been previously acquired, please run first self.acquire." %(key_list)
	
	def getDistributionOfElements(self, key_list, ):
		"""
		For each element of the bed file, count how many time it is shared among datasets and store values into a returning dictionary of points: k = point val from bed tuple; v = count (<- this count MUST be at most len(key_list)!!)
		Returns a dictionary with counts: max count MUST be equal to the number of IDs (len(key_list))
		"""
		all_valid_points = True
		# check key existance
		for k in key_list:
			if not self.beddata.has_key(k):
				all_valid_points = False
		# once everything is ok, continue
		if all_valid_points:
			locations_count = {}
			for k in key_list:
				for point in self.beddata[k]:
					if locations_count.has_key(point):
						locations_count[point] += 1
					else:
						locations_count[point] = 1
			return locations_count
		else:
			print "[AP]\tErrors in key list (give: %s): some of then have not been previously acquired, please run first self.acquire." %(key_list)
	
	def getElementsWithSequenceCountAbove(self, key_list, sequencecount_threshold = 1):
		"""
		For each element of the bed file, update out var self.beddata_out_sequencecount.
		"""
		all_valid_points = True
		# check key existance
		for k in key_list:
			if not self.beddata.has_key(k):
				all_valid_points = False
		# once everything is ok, continue
		nothing_changed = False # tell me if you found sequence counts or not
		if all_valid_points:
			for k in key_list: # for each input key
				# check output var, if already exists or not
				if not self.beddata_out_sequencecount.has_key(k): # put it in first
					self.beddata_out_sequencecount[k] = {}
				for point, details in self.beddata[k].iteritems(): # loop all points
					if details.has_key('options'): # find optional sequence count var
						if str(details['options'][1]).isdigit() and details['options'][1] > sequencecount_threshold: # if this number is digit and above threshold
							self.beddata_out_sequencecount[k][point] = details
					else:
						nothing_changed = True
		else:
			print "[AP]\tErrors in key list (give: %s): some of then have not been previously acquired, please run first self.acquire." %(key_list)
		# print warning
		if nothing_changed:
			print "[AP]\tWarning: in key list (give: %s) some of the points have not options available so sequence count NOT found!" %(key_list)
		

class BedFileHandler:
	"""
	BED File handler class
	"""
	def __init__(self, ):
		"""
		Init class
		"""
		self.bedfile = None # bed file [INPUT]
		self.beddata = None # TSV data, i.e. saved into DB or read from file [INPUT]
		self.bedtuple = None # tuple array of locations in bed format -> chr, start, end, optional fields
		self.locations = {} # k = (chr number, start, end); v = dict: k1 = 'field id'; v1 = this data (chr:: string, start:: Number, end:: Number, options[:: this is TUPLE]) [OUTPUT]
		self.data2write = None # iterable element (list or set) of arrays (values to write)
		#print "[AP]\tBED file handler object initialized."
	
	def getLocations(self, ):
		"""
		From file BED format (http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
		get locations and put them all in the dictionary self.locations
		"""
		if self.bedfile is not None:
			with open(self.bedfile, 'rb') as inf:
				reader = csv.reader(inf, delimiter = "\t")
				try:
					for row in reader:
						if len(row) > 3:
							self.locations[( row[0].strip().split('chr')[-1], int(row[1].strip()), int(row[2].strip()) )] = {
									'chr': row[0].strip().split('chr')[-1], 
									'start': int(row[1].strip()), 
									'end': int(row[2].strip()), 
									'options': row[3:],
									}
						else:
							self.locations[( row[0].strip().split('chr')[-1], int(row[1].strip()), int(row[2].strip()) )] = {
									'chr': row[0].strip().split('chr')[-1], 
									'start': int(row[1].strip()), 
									'end': int(row[2].strip()), 
									'options': (),
									}
				except csv.Error, e:
					sys.exit("[AP]\tError while reading file %s, line %d: %s" % (infile, reader.line_num, e))
		elif os.path.isfile(self.bedfile):
			print "[AP]\tWarning: you must have the object.bedfile with a valid path; try again.."
		else:
			print "[AP]\tWarning: you must have the object.bedfile not None; try again after data filled."
	
	
	def getLocationsFromData(self, ):
		"""
		From read file data (or DB data in TSV format) in BED format (http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
		get locations and put them all in the dictionary self.locations
		This function could be used also to UPDATE locations!
		"""
		if self.beddata is not None and len( self.beddata.strip().split('\n') ) > 0:
			for row in self.beddata.strip().split('\n'):
				fields = row.strip().split('\t')
				if len(fields) > 3:
					self.locations[( fields[0].strip().split('chr')[-1], int(fields[1].strip()), int(fields[2].strip()) )] = {
							'chr': fields[0].strip().split('chr')[-1], 
							'start': int(fields[1].strip()), 
							'end': int(fields[2].strip()), 
							'options': fields[3:],
							}
				elif len(fields) == 3:
					self.locations[( fields[0].strip().split('chr')[-1], int(fields[1].strip()), int(fields[2].strip()) )] = {
							'chr': fields[0].strip().split('chr')[-1], 
							'start': int(fields[1].strip()), 
							'end': int(fields[2].strip()), 
							'options': (),
							}
				else:
					print "[AP]\tWarning: This ID has no BED values (at least 3 columns)."
					#pass
		else:
			print "[AP]\tWarning: you must have the object.beddata not None; try again after data filled."
	
	def getLocationsFromTuple(self, ):
		"""
		From read data from tuple array of ISs where the tuple structure is bed-like: chr, start, end, opts;
		get locations and put them all in the dictionary self.locations
		This function could be used also to UPDATE locations!
		"""
		if self.bedtuple is not None and len( self.bedtuple ) > 0:
			for row in self.bedtuple:
				if len(row) > 3:
					chr = row[0]
					start = row[1]
					end = row[2]
					opts = tuple(row[3:])
					self.locations[( chr, start, end )] = {
							'chr': chr, 
							'start': start, 
							'end': end, 
							'options': opts,
							}
				elif len(row) == 3:
					chr = row[0]
					start = row[1]
					end = row[2]
					self.locations[( chr, start, end )] = {
							'chr': chr, 
							'start': start, 
							'end': end, 
							'options': (),
							}
				else:
					print "[AP]\tWarning: This location ID has no BED values."
					#pass
		else:
			print "[AP]\tWarning: you must have the object.beddata not None; try again after data filled."

	def editLocationsEndPosition(self, bpToEdit, orientationBased = False):
		"""
		Given the self.locations data, use this function to expand the end position of each ISs, ACCORDING TO ORIENTATION (if available)
		This function UPDATEs self.locations
		bpToEdit must be an integer (positive or negative!!!)
		"""
		if self.locations.keys() is not None:
			print "[AP]\tUpdating bed file locations: end points editing (according to orientation if available)."
			newlocations = {}
			for (chr, start, end), bedvals in self.locations.iteritems():
				# se esistono opzioni e tra queste c'e' orientamento, allora procedi concorde all'orientamento, else vale sempre +
				if len(bedvals['options']) is not 0:
					orientation = bedvals['options'][2] # orientation is always in 6th position
					if orientationBased: # if orientation affects extension
						if orientation == '2':
							newlocations[(chr, start, int(end) - bpToEdit)] = {
									'chr': chr, 
									'start': start, 
									'end': int(end) - bpToEdit, 
									'options': bedvals['options'],
									}
						else:
							newlocations[(chr, start, int(end) + bpToEdit)] = {
								'chr': chr, 
								'start': start, 
								'end': int(end) + bpToEdit, 
								'options': bedvals['options'],
								}
					else: # in case of avoiding orientation effect
						newlocations[(chr, start, int(end) + bpToEdit)] = {
							'chr': chr, 
							'start': start, 
							'end': int(end) + bpToEdit, 
							'options': bedvals['options'],
							}
				else: # default (with no orientation info) -> use +
					newlocations[(chr, start, int(end) + bpToEdit)] = {
							'chr': chr, 
							'start': start, 
							'end': int(end) + bpToEdit, 
							'options': bedvals['options'],
							}
			
			self.locations = dict(newlocations) # create a new instance of self.locations
		else:
			print "[AP]\tWarning: you must have the object.locations not None; try again after data filled."

	def write(self, outfile):
		"""
		Write BED file
		"""
		if self.data2write is not None:
			try:
				writer = csv.writer(open(outfile, 'wb'), delimiter = "\t")
				for d in self.data2write:
					writer.writerow(d)
					
			except csv.Error, e:
				sys.exit("[AP]\tError while writing CSV BED file %s. Check paths and names." % (outfile))
		else:
			print "[AP]\tWarning: you must have the object.data2write not None; try again after data filled."

	
	def bedToISs(self, ):
		"""
		From BED file to Integration site list of THIS bed ONLY!!
		"""
		print "[DB]\tComputing unique integration loci extraction (+/- 3 bp) from DB"
		dict_is = {}
		
		integration_locus = None
		i = 0
		# sort locations first
		issk = self.locations.keys()
		issk_sorted = sorted(issk, key=itemgetter(0,1), reverse = False)
					
		for i in range(0,len(issk_sorted)):
			# se la diff tra integration loci e' <= 3 E l'attuale locus in analisi ha una differenza <= 3 con il primo allora, inseriscilo nel dizionario
			if i > 0 and integration_locus is not None and ( abs(self.locations[issk_sorted[i]]["start"] - self.locations[issk_sorted[i-1]]["start"]) <= 3 and abs(integration_locus - self.locations[issk_sorted[i]]["start"]) <= 3 ):
				dict_is[issk_sorted[i]].append( self.locations[issk_sorted[i]] )
			else:
				integration_locus = int(self.locations[issk_sorted[i]]["start"])
				dict_is[ issk_sorted[i] ] = [self.locations[issk_sorted[i]]]
		
		return dict_is
				

class SaturationMonitor:
	"""
	Saturation monitor DB class: handles DB connections and queries
	"""
	def __init__(self, jobId, cursor, schema, table, comment = ""):
		"""
		Init class
		"""
		self.cursor = cursor
		self.table = table
		self.schema = schema
		self.jobId = jobId
		self.comment = comment
		print "[AP]\tSaturation DB is now active. Connected to Tracking job", self.jobId
		
	def insert(self, datafields = [], textdatafields = [], ):
		"""
		Using connection.cursor, insert data into table. [B_binId][N_setId]
		
		Old version:
			query =
				INSERT INTO  `%s`.`%s` (
					`jobId` ,
					`B_binId` ,
					`totalReads` ,
					`binReads` ,
					`N_setId` ,
					`targetTable` ,
					`readsCount`
					)
				VALUES (
					'%s', '%s', '%s', '%s', '%s',  '%s',  '%s'
					);

			      %(self.schema,
				    self.table, 
				    self.jobId, 
				    B_binId, 
				    totalReads, 
				    binReads, 
				    N_setId, 
				    targetTable, 
				    readsCount
				    )
		"""
		if len(textdatafields) > 0:
			query = """
				INSERT INTO  `%s`.`%s` 
				VALUES (
					'%s',  %s, %s
					);
			      """ %(self.schema,
				    self.table, 
				    self.jobId, 
				    "'" + "', '".join( str(f) for f in datafields ) + "'",
				    '"' + '", "'.join( str(f) for f in textdatafields ) + '"',
				    )
		else:
			query = """
				INSERT INTO  `%s`.`%s` 
				VALUES (
					'%s',  %s
					);
			      """ %(self.schema,
				    self.table, 
				    self.jobId, 
				    "'" + "', '".join( str(f) for f in datafields ) + "'",
				    )
		#print query
		self.run(query)
	
	def close(self):
		"""
		Destroy all
		"""
		pass
	
	def createTable(self, fields = [], textfields = []):
		"""
		Create table if not exists. [B_binId][N_setId]
		First version:
			query =
			      CREATE TABLE IF NOT EXISTS `%s`.`%s` (
				`jobId` BIGINT( 15 ) NOT NULL COMMENT  'date +"\%%Y\%%m\%%d\%%H\%%M\%%S"',
				`B_binId` INT( 4 ) NOT NULL COMMENT  'percentage on total ',
				`totalReads` BIGINT( 15 ) NOT NULL COMMENT  'number of total reads available from initial sample size (i.e. raw table)',
				`binReads` BIGINT( 15 ) NOT NULL COMMENT  'number of reads for this bin id, the fraction of total (=bin percentage)',
				`N_setId` INT( 10 ) NOT NULL COMMENT  'iD of subset n in N of this bootstrap set (bin id); default: 1..500',
				`targetTable` VARCHAR( 150 ) NOT NULL COMMENT  'target table name',
				`readsCount` BIGINT( 15 ) NOT NULL COMMENT  'count of reads for this n in N, bin ID and table.',
				PRIMARY KEY (  `B_binId` ,  `N_setId` ,  `targetTable` )
			      ) ENGINE = MYISAM CHARACTER SET utf8 COLLATE utf8_general_ci COMMENT =  'saturation analysis, las update %s';
			%(	self.schema,
				self.table,
				strftime("%Y-%m-%d %H:%M:%S", gmtime()),
		
		Second version (dynamic)
			query = 
			      CREATE TABLE IF NOT EXISTS `%s`.`%s` (
				`jobId` BIGINT( 15 ) NOT NULL COMMENT  'date +"\%%Y\%%m\%%d\%%H\%%M\%%S"',
				`B_binId` INT( 4 ) NOT NULL COMMENT  'percentage on total ',
				`totalReads` BIGINT( 15 ) NOT NULL COMMENT  'number of total reads available from initial sample size (i.e. raw table)',
				`binReads` BIGINT( 15 ) NOT NULL COMMENT  'number of reads for this bin id, the fraction of total (=bin percentage)',
				`N_setId` INT( 10 ) NOT NULL COMMENT  'iD of subset n in N of this bootstrap set (bin id); default: 1..500',
				[`NEWFIELD` BIGINT( 15 ) NOT NULL]*,
				PRIMARY KEY (  `B_binId` ,  `N_setId` ,  `targetTable` )
			      ) ENGINE = MYISAM CHARACTER SET utf8 COLLATE utf8_general_ci COMMENT =  'saturation analysis, las update %s';
			%(	self.schema,
				self.table,
				strftime("%Y-%m-%d %H:%M:%S", gmtime()),
		"""
		print "[AP]\tSaturation DB: table creation (schema, table):", self.schema, self.table
		# combine fields into a string
		# for real values
		field_template = Template("`$field` FLOAT NOT NULL,\n")
		newfield = ""
		for k in fields:
			newfield += field_template.substitute(field = k)
		# now for note/comments or other long data
		textfield_template = Template("`$field` LONGTEXT NULL,\n")
		newtextfield = ""
		for k in textfields:
			newtextfield += textfield_template.substitute(field = k)
		# create query
		query = """
		      CREATE TABLE IF NOT EXISTS `%s`.`%s` (
			`jobId` BIGINT( 15 ) NOT NULL COMMENT  'date +"\%%Y\%%m\%%d\%%H\%%M\%%S"',
			`B_binId` INT( 4 ) NOT NULL COMMENT  'percentage on total ',
			`N_setId` INT( 10 ) NOT NULL COMMENT  'iD of subset n in N of this bootstrap set (bin id); default: 1..500',
			%s %s
			INDEX ( `B_binId` ,  `N_setId` )
		      ) ENGINE = MYISAM CHARACTER SET utf8 COLLATE utf8_general_ci COMMENT =  '%s. %s';
		""" %(	self.schema,
			self.table,
			newfield,
			newtextfield,
			strftime("%Y-%m-%d %H:%M:%S", gmtime()),
			self.comment,
		      )
		self.run(query)
	
	def truncateTable(self, ):
		query = """
		      TRUNCATE TABLE `%s`.`%s`;
		""" %(	self.schema,
			self.table,
			    )
		self.run(query)
		
	def run(self, query):
		try:
			#print query
			self.cursor.execute(""" %s """ %(query) )
		except MySQLdb.Error, e:
			print "DATABASE Error", e.args, e
			sys.exit (1)
		

class eMail:
	"""
	Instance of email handler
	"""
	def __init__(self, jobId, cursor, schema, table):
		"""
		Init class
		"""
		from email.mime.image import MIMEImage
		from email.mime.multipart import MIMEMultipart
		from email.mime.text import MIMEText
		print "[AP]\tSending notice by email."
		
	def sendEmailMessage(self, app_name, sender_addr, destination_addr):
		# Create the container (outer) email message.
		msg = MIMEMultipart()
		msg['Subject'] = "App. %s Status Results Available (finished job)" %(app_name)
		# me == the sender's email address
		# family = the list of all recipients' email addresses
		COMMASPACE = ', '
		msg['From'] = sender_addr
		#msg['To'] = COMMASPACE.join(family)
		msg['To'] = destination_addr
		msg.preamble = """
		This is an automatic message from Alienware server.
		The application %(app_name)s finished. Now you can download results.
		For any support, ask to Andrea Calabria (reply to this email).
		Regards
		Andrea
		""" %{
		  'app_name': app_name,
		  'support': "Andrea Calabria",
		  }

		# Send the email via our own SMTP server.
		s = smtplib.SMTP('mail.hsr.it')
		s.sendmail(sender_addr, destination_addr, msg.as_string())
		s.quit()
