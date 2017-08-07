#!/usr/bin/python
# -*- coding: utf-8 -*-
import MySQLdb, sys, os, re, argparse, csv
import smtplib, numpy, scipy.stats
from operator import itemgetter, attrgetter
from string import Template
#from Bio import SeqIO
from time import gmtime, strftime

# my classes
#from Tiget import JobMonitor, Experiment, IntegrationLocus

header = """

+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		December 2012
 Contact:	andrea.calabria@hsr.it
 Revision:	0.2
+--------------------------------------------------+

 Description
  - This program is meant to be run for each TAG of LAM during the pipeline
    analysis. 
  - From BED file of each single LAM and from association file
    (from which I exctract info of tissue, sample and time point),
    it insert data into user defined redundant table (if it does 
    not exist, it creates the table and then inserts data, else 
    it updates table with new data).
  - Minimum required info to store into redundant DB (+ TST):
    'chr,integration_locus,group_name,pool,tag,strand,refSeqName,refSeqSymbol'
    where annotations also are not mandatory or could be updated later on.
  - Mandatory input data: patient, pool, barcode, target redundant table,
    and other DB options

 Note:
  - Association file: format must be TSV with these columns:
    <barcode_id> <barcode_sequence> <tissue> <sample> <time point> <lam_id> 
    <lam_fullname> <marker> <enzyme> <vector>
  
 TODO:
  - add data such as raw and trimmed sequence -> ad hoc script
  - add data of annotation (gene annotation from GREAT)

 Steps
 	1. parse association file to extract 'tissue, sample and time point data'
	2. parse BED to extract 'chr, inte site, strand' 
	   (this is the real starting point based on strand)
	3. create DB table (if required)
	4. import data into table
  
""" 

description = "This application will import BED ISs into redundant DB table."

usage_example = """
Examples of usage:
 (1) Extract data from MLD01 patient, pool12v2, tag ACCTAC: 
    APP -b /opt/NGS/results/MLD/MLD01/bed/pool12v2/MLD01_pool12v2_CD19PB_18_ACCTAC.noAB.noLTRLC.sorted.rel.pg.dupflag.tagfilter.iss.bed -a /home/andrea/ISAtools/NGS/association/asso.mld01.pool12.tsv --patient MLD01 --pool "pool12v2" --tag ACCTAC --dbschema sequence_mld01 --dbtable miseqtest

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ Import into DB redundant ISs from BED ] \n", description = description)
parser.add_argument('-b', '--bedfile', dest="bedfile", help="BED file to process. No default option.", action="store", required=True)
parser.add_argument('-a', '--associationfile', dest="associationfile", help="File of association: for each barcode this file (TSV format) specifies TAG sequence, tissue, sample and time point, lam id, lam name, marker, enzyme, vector. No default option.", action="store", required=True)
#parser.add_argument('--headerformat', dest="headerformat", help="Header format of the filter file. Cases: type 1:: '@M00174:25:000000000-A0T21:1:1:15041:1491 1:N:0:0'; type 2:: '@M00571:5:000000000-A267F:1:1101:14475:1610/1'. Select type number {1, .., N}. No default value assigned because you must explicitely be aware of what you are doing.", action="store", required=True)
parser.add_argument('-r', '--seqfileraw', dest="seqfileraw", help="File of raw sequences, extracted with a custom script of Andrea that formats for each read ID, as reported in the BED file, it stores raw sequences. Default is 'None'.", action="store", default=None)
parser.add_argument('-t', '--seqfiletrimmed', dest="seqfiletrimmed", help="File of trimmed sequences, extracted with a custom script of Andrea that formats for each read ID, as reported in the BED file, it stores trimmed sequences. Default is 'None'.", action="store", default=None)
parser.add_argument('--patient', dest="patient", help="Patient ID. No default option.", action="store", required=True)
parser.add_argument('--pool', dest="pool", help="Pool ID. No default option.", action="store", required=True)
parser.add_argument('--tag', dest="tag", help="TAG ID. No default option.", action="store", required=True)
parser.add_argument('-d', '--db', dest="userdb_choice", action="store", help="DB to be used; accepted values: {local, cluster, nas}. Default: local", default='local')
parser.add_argument('--dbschema', dest="dbschema", help="Target Database schema in which importing data. E.g.: MLD01. Mandatory.", action="store", required=True)
parser.add_argument('--dbtable', dest="dbtable", help="Target Database table in which importing data. E.g.: redundant_MLD01_ALL_NEW_NAMES. Mandatory.", action="store", required=True)
parser.add_argument('--writecsv', dest="writecsv", help="Write CSV file instead of importing data into MySQL. This will return a file with the same name of table. E.g.: redundant_MLD01_ALL_NEW_NAMES. Mandatory.", action="store_true", required=False, default=False)
parser.add_argument('--csvfilename', dest="csvfilename", help="If enabled the option writecsv, this file will be used as CSV. Default: datafile.csv or tablename", action="store", default="datafile.csv")

args = parser.parse_args()

########################################################################
####### GLOBAL VARS
########################################################################

# init db values
userdb = {}
if args.userdb_choice == 'local':
	userdb = {
		'host': "localhost",
		'user': "andrea",
		'passwd': "andrea",
		'db': args.dbschema,
		}
elif args.userdb_choice == 'cluster':
	userdb = {
		'host': "172.21.90.125",
		'user': "calabria",
		'passwd': "c4l4br14",
		'db': args.dbschema,
		}
elif args.userdb_choice == 'nas':
	userdb = {
		'host': "172.25.39.41",
		'user': "integration",
		'passwd': "5aSdXhS9Gpnndc24",
		'db': args.dbschema,
		}
else:
	print "ERROR: Wrong DB. You typed ", args.userdb_choice, " instead of valid options (see help for list)."
	sys.exit()

conn = MySQLdb.connect( 
	host = userdb["host"],
	user = userdb["user"],
	passwd = userdb["passwd"],
	db = userdb["db"]
	)
#cursor = conn.cursor ()
cursor = conn.cursor (MySQLdb.cursors.DictCursor)


########################################################################
####### FUNCTIONS
########################################################################

def parseAssociationfile(infile):
	"""
	IN: association file
	OUT: dictionary of TAGs (k=tag, v=dict: k={tissue, sample, timepoint}, v=values)
	FILE FORMAT: <barcode_id> <barcode_sequence> <tissue> <sample> <timepoint> <lam_id> <lam_fullname> <marker> <enzyme> <vector>
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
					}
			else:
				print "[AP] WARNING: this file contais a blank row!!!! Check it please."
	return assodict

def parseBEDfile(inbedfile):
	"""
	IN: bed file
	OUT: BED tuple array with chr, is, header, strand, score
	LOGICS: 	
		1. recognize start point from BED
		2. create tuple array of BED elements
	"""
	bedarray = [] # return var
	with open(inbedfile, 'rb') as inf:
		reader = csv.reader(inf, delimiter = "\t")
		for row in reader:
			chr = row[0].replace('chr', '')
			if row[5] == '+':
				bedarray.append( (chr, row[1], row[3], row[5], row[4]) )
			else:
				bedarray.append( (chr, row[2], row[3], row[5], row[4]) )
	results = sorted(bedarray, key=itemgetter(0, 1), reverse = False) # sorted by chr, start, end -> now do sort by orientation!!!
	return results
	
def queryIfTableExists(dbschema, dbtable):
	"""
	IN: db options
	OUT: boolean (true/false)
	LOGICS: 	
		1. query info schema to see if table exists
		2. return response
	"""
	exists = False
	query = """
		SELECT *
		FROM information_schema.tables
		WHERE table_schema = '%s'
			AND table_name = '%s';
		""" %(dbschema, dbtable)
	cursor.execute(query)
	results = cursor.fetchall()
	if len(results) == 1:
		exists = True
	return exists

def createTargetTable(dbschema, dbtable):
	"""
	IN: db options
	OUT: -
	LOGICS: if target table does NOT exist, create table
	"""
	query = """
		CREATE TABLE IF NOT EXISTS `%s`.`%s` (
		  `ID` int(11) NOT NULL AUTO_INCREMENT,
		  `group_name` varchar(200) DEFAULT NULL,
		  `n_LAM` varchar(200) DEFAULT NULL,
		  `pool` varchar(200) DEFAULT NULL,
		  `tag` varchar(20) DEFAULT NULL,
		  `sample` varchar(200) DEFAULT NULL,
		  `vector` varchar(50) DEFAULT NULL,
		  `tissue` varchar(50) DEFAULT NULL,
		  `treatment` varchar(20) DEFAULT NULL,
		  `enzyme` varchar(50) DEFAULT NULL,
		  `complete_name` varchar(200) DEFAULT NULL,
		  `header` varchar(200) DEFAULT NULL,
		  `sequence_raw` varchar(2000) DEFAULT NULL,
		  `sequence_trimmed` varchar(2000) DEFAULT NULL,
		  `chr` varchar(50) DEFAULT NULL,
		  `integration_locus` int(11) DEFAULT NULL,
		  `sequence_count` int(3) NOT NULL DEFAULT '1',
		  `start_seq` int(11) DEFAULT NULL,
		  `end_seq` int(11) DEFAULT NULL,
		  `score` double DEFAULT NULL,
		  `identity` double DEFAULT NULL,
		  `e_value` double DEFAULT NULL,
		  `strand` varchar(2) DEFAULT NULL,
		  `q_size` int(11) DEFAULT NULL,
		  `span` int(11) DEFAULT NULL,
		  `transcript_ucsc` varchar(20) DEFAULT NULL,
		  `orientation_ucsc` varchar(20) DEFAULT NULL,
		  `distance_to_TSS_ucsc` int(11) DEFAULT NULL,
		  `in_gene_ucsc` double DEFAULT NULL,
		  `upstream_end_ucsc` int(11) DEFAULT NULL,
		  `downstream_end_ucsc` int(11) DEFAULT NULL,
		  `gene_lenght_ucsc` int(11) DEFAULT NULL,
		  `ref_gene_ucsc` varchar(200) DEFAULT NULL,
		  `ref_seq_name_ucsc` varchar(200) DEFAULT NULL,
		  `entrez_gene` varchar(200) DEFAULT NULL,
		  `entrez_name` varchar(200) DEFAULT NULL,
		  `ensembl_gene` varchar(200) DEFAULT NULL,
		  `refSeqName` varchar(200) DEFAULT NULL,
		  `refSeqSymbol` varchar(200) DEFAULT NULL,
		  `gene_orientation_refseq` varchar(2) DEFAULT NULL,
		  `ingene_refseq` double DEFAULT NULL,
		  `distance_to_TSS_refseq` int(11) DEFAULT NULL,
		  `upstream_end_refseq` int(11) DEFAULT NULL,
		  `downstream_end_refseq` int(11) DEFAULT NULL,
		  `gene_lenght_refseq` int(11) DEFAULT NULL,
		  `miRNA` varchar(200) DEFAULT NULL,
		  `orientation_mirna` varchar(200) DEFAULT NULL,
		  `inmirna` double DEFAULT NULL,
		  `distance_to_TSS_mirna` int(11) DEFAULT NULL,
		  `upstream_end_mirna` int(11) DEFAULT NULL,
		  `downstream_end_mirna` int(11) DEFAULT NULL,
		  `lenght_mirna` int(11) DEFAULT NULL,
		  PRIMARY KEY (`ID`),
		  KEY `integration_locus` (`integration_locus`),
		  KEY `header` (`header`)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;
		""" %(dbschema, dbtable)
	cursor.execute(query)
	return True

def aquireSeqData(infile):
	"""
	IN: seqfile formatted as <id> <sequence>
	OUT: dictionary of data from CSV file
	LOGICS: 
	"""
	if not os.path.isfile(infile):
		print "[ERROR]\tNo valid input files. Check them please then run me again.\n"
		sys.exit()
	else:
		dicseq = {}
		with open(infile, 'rb') as inf:
			reader = csv.reader(inf, delimiter = "\t")
			for row in reader:
				if len(row)>1:
					dicseq[row[0]] = row[1].strip()
		return dicseq

######################################################					


print "[AP]\tChecking input files"
if not os.path.isfile(args.bedfile) or not os.path.isfile(args.associationfile):
	print "[ERROR]\tNo valid input files. Check them please then run me again.\n"
	sys.exit()

print "[AP]\tGet Association data"
dict_asso = parseAssociationfile(args.associationfile)

print "[AP]\tGet BED data"
beddata = parseBEDfile(args.bedfile)

# if you selected to write a CSV file instead of importing data into DB
if args.writecsv:
	print "[AP]\tCreate target FILE"
	destfile = None
	if not os.path.isfile(args.csvfilename):
		if args.csvfilename == "datafile.csv":
			destfile = args.dbschema + "." + args.dbtable + ".csv"
		else:
			destfile = args.csvfilename
	else:
		print "[AP]\tOverwriting existing data file!!!", args.csvfilename

	fo = open(destfile, 'w')

	## in case of having or not sequence file, change import method
	if args.seqfileraw is not None and args.seqfiletrimmed is not None:
		#print "[AP]\tGet SEQUENCE data (raw and trimmed)"
		#dict_seq_raw = aquireSeqData(args.seqfileraw)
		#dict_seq_trimmed = aquireSeqData(args.seqfiletrimmed)
		
		print "[AP]\tImport data into File"
		### now format data to import all and insert them all
		query_import = ""
		index = 1
		# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
		for (chr, locus, header, strand, score) in beddata:
			query_import += "\t".join( str(x) for x in  
				[index,
					args.patient, 
					dict_asso[args.tag]["lam_id"], 
					args.pool, 
					args.tag, 
					dict_asso[args.tag]["sample"], 
					dict_asso[args.tag]["vector"],
					dict_asso[args.tag]["tissue"], 
					dict_asso[args.tag]["timepoint"],
					dict_asso[args.tag]["enzyme"],
					dict_asso[args.tag]["lam_fullname"],
					header,
					"NULL","NULL",
					chr,
					locus,
					"NULL","NULL","NULL",
					score,
					"NULL","NULL",
					strand,
					"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL"
					]
			  ) + "\n"
			index += 1
		#print query_import.rstrip(', ')
		fo.write(query_import)

	else: # caso of not available sequence file
		print "[AP]\tImport data into target Table"
		### now format data to import all and insert them all
		step = 5000
		bedstart = 0 # starting point of bed data array for each step
		bedend = 0 # starting point of bed data array for each step
		index = 1
		# check max len of beddata to set bedend, no more than beddata len!
		if len(beddata) > step:
			bedend = step
		else:
			bedend = len(beddata)
		# now process data into steps of a maximum number of blocks (step)
		for elems in range(0, len(beddata), step):
			query_import = ""
			# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
			for (chr, locus, header, strand, score) in beddata[bedstart:bedend]:
				query_import += "\t".join( str(x) for x in  
					[index,
					args.patient, 
					dict_asso[args.tag]["lam_id"], 
					args.pool, 
					args.tag, 
					dict_asso[args.tag]["sample"], 
					dict_asso[args.tag]["vector"],
					dict_asso[args.tag]["tissue"], 
					dict_asso[args.tag]["timepoint"],
					dict_asso[args.tag]["enzyme"],
					dict_asso[args.tag]["lam_fullname"],
					header,
					"NULL","NULL",
					chr,
					locus,
					"1","NULL","NULL",
					score,
					"NULL","NULL",
					strand,
					"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL",	]
				  ) + "\n"
				index += 1
			#print query_import.rstrip(', ')
			fo.write(query_import)
			
			# increment bed indexes
			bedstart += step
			bedend += step

	# at the end, clos ethe file
	fo.close()


else: ## branch of args.writecsv, case write into DB
	print "[AP]\tCreate target Table if not exists"
	createTargetTable(args.dbschema, args.dbtable)

	## in case of having or not sequence file, change import method
	if args.seqfileraw is not None and args.seqfiletrimmed is not None:
		print "[AP]\tGet SEQUENCE data (raw and trimmed)"
		dict_seq_raw = aquireSeqData(args.seqfileraw)
		dict_seq_trimmed = aquireSeqData(args.seqfiletrimmed)
		
		print "[AP]\tImport data into target Table"
		### now format data to import all and insert them all
		query_import = """
			INSERT INTO `%s`.`%s` (`group_name`, `n_LAM`, `pool`, `tag`, `sample`, `tissue`, `treatment`, `enzyme`, `complete_name`, `header`, `chr`, `integration_locus`, `strand`, `score`, `vector`, `sequence_raw`, `sequence_trimmed`) 
			VALUES 
			""" %(args.dbschema, args.dbtable)
		# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
		for (chr, locus, header, strand, score) in beddata:
			query_import += " ('" + "', '".join( str(x) for x in  
				[args.patient, 
				dict_asso[args.tag]["lam_id"], 
				args.pool, 
				args.tag, 
				dict_asso[args.tag]["sample"], 
				dict_asso[args.tag]["tissue"], 
				dict_asso[args.tag]["timepoint"],
				dict_asso[args.tag]["enzyme"],
				dict_asso[args.tag]["lam_fullname"],
				header,
				chr,
				locus,
				strand,
				score,
				dict_asso[args.tag]["vector"],
				dict_seq_raw[header],
				dict_seq_trimmed[header]		 ]
			  ) + "'), "
		#print query_import.rstrip(', ')
		cursor.execute( query_import.rstrip(', ') )

	else: # caso of not available sequence file
		print "[AP]\tImport data into target Table"
		### now format data to import all and insert them all
		step = 5000
		bedstart = 0 # starting point of bed data array for each step
		bedend = 0 # starting point of bed data array for each step
		# check max len of beddata to set bedend, no more than beddata len!
		if len(beddata) > step:
			bedend = step
		else:
			bedend = len(beddata)
		# now process data into steps of a maximum number of blocks (step)
		for elems in range(0, len(beddata), step):
			query_import = """
				INSERT INTO `%s`.`%s` (`group_name`, `n_LAM`, `pool`, `tag`, `sample`, `tissue`, `treatment`, `enzyme`, `complete_name`, `header`, `chr`, `integration_locus`, `strand`, `score`, `vector`) 
				VALUES 
				""" %(args.dbschema, args.dbtable)
			# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
			for (chr, locus, header, strand, score) in beddata[bedstart:bedend]:
				query_import += " ('" + "', '".join( str(x) for x in  
					[args.patient, 
					dict_asso[args.tag]["lam_id"], 
					args.pool, 
					args.tag, 
					dict_asso[args.tag]["sample"], 
					dict_asso[args.tag]["tissue"], 
					dict_asso[args.tag]["timepoint"],
					dict_asso[args.tag]["enzyme"],
					dict_asso[args.tag]["lam_fullname"],
					header,
					chr,
					locus,
					strand,
					score,
					dict_asso[args.tag]["vector"] ]
				  ) + "'), "
			#print query_import.rstrip(', ')
			cursor.execute( query_import.rstrip(', ') )
			
			# increment bed indexes
			bedstart += step
			bedend += step

print "[AP]\tTask finished\n\tCheck target table\n\n"


######################## finalize objects ##################		
cursor.close()
conn.close()
