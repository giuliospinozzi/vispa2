#!/usr/bin/python
# -*- coding: utf-8 -*-
# to run system calls
from subprocess import call
import math
import sys, os, re, argparse, csv, time, random
from time import gmtime, strftime
# to process sequence data
from Bio import SeqIO
# to compute stats
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg') ## use background plot generation, avoid front end window opening
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pandas import DataFrame, read_csv

header = """

+--------------------------------------------------+
 Author:	Andrea Calabria
 Date:		September 2014
 Contact:	andrea.calabria@hsr.it
 Revision:	0.1
+--------------------------------------------------+

 Description
  - This program is aimed at optimizing pipeline steps, here the TRIMMING step.
  - Input: CSV file with columns: header, full sequence, ltr sequence, label [0,1]
  	 
 TODO:
  - Nothing yet

 STEPS:
  - acquire input dataset -> dictionary of elements with ground truth
  - compute trimming -> for each program: for each parameter configuration: run trimming application (output fasta/q), acquire results in memory, compute junction delta, append results to header dictionary as obtained results (delta will be considered in a custom interval as correctly retrieved or not, 0..K)
  - compute statistical measures -> parse dictionary of results vs ground truth to extract precision, recall, specificty, sensitivity, FDR
  - results -> return best results/solutions, plot stat metrics of best results (within or not all sulutions)

 NOTE [ITA]:
  - il path test e' splittato da quello stats e dipende dalla presenza del file --outMatrix (guarda checkArgs funtion)
""" 

description = "This application optimizes trimming procedure"

usage_example = """
Examples of usage:
 (1) #nohup /home/andrea/trimming_analysis/./do.py --fastaSeqs /home/andrea/trimming_analysis/32_LTR_trimming.fa --groundTruth /home/andrea/trimming_analysis/32_LTR_ground_truth.csv --genomicSeq AGAGATCAAGTCTCACTATGTTGCCCAGGTTGGTCTCGAACTCTTGGGCTCAAGCAATCCTCTCACCTCAGCCTCCCAAAGTGCTGGGATTACAGACATGAGCCACCCTGCTCGGCTAGAATT --adapterSeq GTGGAAAATCTCTAGCA --tmpfolder tmp_17bp --adapterFasta /home/andrea/trimming_analysis/validation/LTR.17bp.fa --fastqSeqs /home/andrea/trimming_analysis/32_LTR_trimming.fastq --doTests --doStats --outMatrix 17bp_Results.tsv 

"""

print header#, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ Optimize trimming ] \n", description = description)
parser.add_argument('--fastaSeqs', dest="fastaSeqs", help="Input sequences in FASTA format. No default option.", action="store", required=True)
parser.add_argument('--fastqSeqs', dest="fastqSeqs", help="Input sequences in FASTQ format. No default option.", action="store", required=True)
parser.add_argument('--groundTruth', dest="groundTruth", help="Input Ground Truth file with labeled sequences, format: header, full sequence, ltr sequence, label [0,1]. No default option.", action="store", required=True)
parser.add_argument('--genomicSeq', dest="genomicSeq", help="Genomic reference sequence used in the ground truth file. No default option. Write the sequence only, not FASTA format", action="store", required=True)
parser.add_argument('--adapterSeq', dest="adapterSeq", help="Adapter reference sequence used in the ground truth file, usually it is the LTR sequence. Write the sequence only. No default option.", action="store", required=True)
parser.add_argument('--adapterFasta', dest="adapterFasta", help="FASTA file of the adapter reference sequence, usually the LTR sequence. Write absolute paths. No default option.", action="store", required=True)
parser.add_argument('--deltaMax', dest="deltaMax", help="Max delta score interval (junction delta of post-trimming sequences from fasta) to compute, from 0 to this max value. Default: 6, that means simmetric +/-6.", action="store", default=6)
parser.add_argument('--outPrefix', dest="outPrefix", help="Output file prefix. Default: Results", action="store", default="Results")
parser.add_argument('--tmpfolder', dest="tmpfolder", help="Temporary folder where collect all results files from trimming programs. Default: /tmp", action="store", default="/tmp")
parser.add_argument('--outMatrix', dest="outMatrix", help="Output matrix file of all results as tab separated value. Default: Results.tsv", action="store", default="Results.tsv")
parser.add_argument('--doTests', dest="doTests", action="store_true", help="Do you want to perform tests? If yes, the program will run system calls to predefined programs (flexbar and eautils). Default: 'False'; alternative: 'True' or just activate the option.", default=False)
parser.add_argument('--doStats', dest="doStats", action="store_true", help="Do you want to perform statistical analyses of the performed tests? If yes, the program will run stat assessments. Default: 'False'; alternative: 'True' or just activate the option.", default=False)

args = parser.parse_args()

########################################################################
####### GLOBAL VARS
########################################################################
## plotting vars
cnames = {
	'aliceblue':            '#F0F8FF',
	'antiquewhite':         '#FAEBD7',
	'aqua':                 '#00FFFF',
	'aquamarine':           '#7FFFD4',
	'azure':                '#F0FFFF',
	'beige':                '#F5F5DC',
	'bisque':               '#FFE4C4',
	'black':                '#000000',
	'blanchedalmond':       '#FFEBCD',
	'blue':                 '#0000FF',
	'blueviolet':           '#8A2BE2',
	'brown':                '#A52A2A',
	'burlywood':            '#DEB887',
	'cadetblue':            '#5F9EA0',
	'chartreuse':           '#7FFF00',
	'chocolate':            '#D2691E',
	'coral':                '#FF7F50',
	'cornflowerblue':       '#6495ED',
	'cornsilk':             '#FFF8DC',
	'crimson':              '#DC143C',
	'cyan':                 '#00FFFF',
	'darkblue':             '#00008B',
	'darkcyan':             '#008B8B',
	'darkgoldenrod':        '#B8860B',
	'darkgray':             '#A9A9A9',
	'darkgreen':            '#006400',
	'darkkhaki':            '#BDB76B',
	'darkmagenta':          '#8B008B',
	'darkolivegreen':       '#556B2F',
	'darkorange':           '#FF8C00',
	'darkorchid':           '#9932CC',
	'darkred':              '#8B0000',
	'darksalmon':           '#E9967A',
	'darkseagreen':         '#8FBC8F',
	'darkslateblue':        '#483D8B',
	'darkslategray':        '#2F4F4F',
	'darkturquoise':        '#00CED1',
	'darkviolet':           '#9400D3',
	'deeppink':             '#FF1493',
	'deepskyblue':          '#00BFFF',
	'dimgray':              '#696969',
	'dodgerblue':           '#1E90FF',
	'firebrick':            '#B22222',
	'floralwhite':          '#FFFAF0',
	'forestgreen':          '#228B22',
	'fuchsia':              '#FF00FF',
	'gainsboro':            '#DCDCDC',
	'ghostwhite':           '#F8F8FF',
	'gold':                 '#FFD700',
	'goldenrod':            '#DAA520',
	'gray':                 '#808080',
	'green':                '#008000',
	'greenyellow':          '#ADFF2F',
	'honeydew':             '#F0FFF0',
	'hotpink':              '#FF69B4',
	'indianred':            '#CD5C5C',
	'indigo':               '#4B0082',
	'ivory':                '#FFFFF0',
	'khaki':                '#F0E68C',
	'lavender':             '#E6E6FA',
	'lavenderblush':        '#FFF0F5',
	'lawngreen':            '#7CFC00',
	'lemonchiffon':         '#FFFACD',
	'lightblue':            '#ADD8E6',
	'lightcoral':           '#F08080',
	'lightcyan':            '#E0FFFF',
	'lightgoldenrodyellow': '#FAFAD2',
	'lightgreen':           '#90EE90',
	'lightgray':            '#D3D3D3',
	'lightpink':            '#FFB6C1',
	'lightsalmon':          '#FFA07A',
	'lightseagreen':        '#20B2AA',
	'lightskyblue':         '#87CEFA',
	'lightslategray':       '#778899',
	'lightsteelblue':       '#B0C4DE',
	'lightyellow':          '#FFFFE0',
	'lime':                 '#00FF00',
	'limegreen':            '#32CD32',
	'linen':                '#FAF0E6',
	'magenta':              '#FF00FF',
	'maroon':               '#800000',
	'mediumaquamarine':     '#66CDAA',
	'mediumblue':           '#0000CD',
	'mediumorchid':         '#BA55D3',
	'mediumpurple':         '#9370DB',
	'mediumseagreen':       '#3CB371',
	'mediumslateblue':      '#7B68EE',
	'mediumspringgreen':    '#00FA9A',
	'mediumturquoise':      '#48D1CC',
	'mediumvioletred':      '#C71585',
	'midnightblue':         '#191970',
	'mintcream':            '#F5FFFA',
	'mistyrose':            '#FFE4E1',
	'moccasin':             '#FFE4B5',
	'navajowhite':          '#FFDEAD',
	'navy':                 '#000080',
	'oldlace':              '#FDF5E6',
	'olive':                '#808000',
	'olivedrab':            '#6B8E23',
	'orange':               '#FFA500',
	'orangered':            '#FF4500',
	'orchid':               '#DA70D6',
	'palegoldenrod':        '#EEE8AA',
	'palegreen':            '#98FB98',
	'paleturquoise':        '#AFEEEE',
	'palevioletred':        '#DB7093',
	'papayawhip':           '#FFEFD5',
	'peachpuff':            '#FFDAB9',
	'peru':                 '#CD853F',
	'pink':                 '#FFC0CB',
	'plum':                 '#DDA0DD',
	'powderblue':           '#B0E0E6',
	'purple':               '#800080',
	'red':                  '#FF0000',
	'rosybrown':            '#BC8F8F',
	'royalblue':            '#4169E1',
	'saddlebrown':          '#8B4513',
	'salmon':               '#FA8072',
	'sandybrown':           '#FAA460',
	'seagreen':             '#2E8B57',
	'seashell':             '#FFF5EE',
	'sienna':               '#A0522D',
	'silver':               '#C0C0C0',
	'skyblue':              '#87CEEB',
	'slateblue':            '#6A5ACD',
	'slategray':            '#708090',
	'snow':                 '#FFFAFA',
	'springgreen':          '#00FF7F',
	'steelblue':            '#4682B4',
	'tan':                  '#D2B48C',
	'teal':                 '#008080',
	'thistle':              '#D8BFD8',
	'tomato':               '#FF6347',
	'turquoise':            '#40E0D0',
	'violet':               '#EE82EE',
	'wheat':                '#F5DEB3',
	'white':                '#FFFFFF',
	'whitesmoke':           '#F5F5F5',
	'yellow':               '#FFFF00',
	'yellowgreen':          '#9ACD32', 
	}
my_cnames = ['red', 'blue', 'green', 'orange', 'pink', 'gold', 'violet', 'cyan', 'forestgreen', 'silver', 'navy']
markers = {u's': u'square', u'^': u'triangle_up', u'o': u'circle', u'v': u'triangle_down', u'<': u'triangle_left', u'>': u'triangle_right'}
my_markers = markers.keys()
##

genoRefLen = len(args.genomicSeq)
adapterSeq = len(args.adapterSeq)
application_dictionary = {
	'flexbar': {
		'exec_name': 'flexbar2.5', # how the software is called in the local server
		'prototype': ['ai', 'ao', 'at'],
		'ai': [-1, -2, -4, -5, -8, -10], # steps
		'ao': [adapterSeq, adapterSeq-2, adapterSeq-4, adapterSeq-6, adapterSeq-8, adapterSeq-10], 
		'at': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.5, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4],
		},
	'eautils': {
		'exec_name': 'fastq-mcf',
		'prototype': ['m', 'p', 'l'],
		'm': [adapterSeq, adapterSeq-2, adapterSeq-4, adapterSeq-6, adapterSeq-8, adapterSeq-10], # Minimum clip length
		'p': [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75], # Maximum adapter difference percentage
		'l': [18], # Minimum remaining sequence length
	}
}
delimiter = "\t"
test_done_header = [] # list of the headers of the tests
delta_to_analyze = range(0,int(args.deltaMax)+1) # values > k in delta become false IS, else true:: example of delta=1, for each tool param config, if a value is >1 (strinctly greather than 1!!!!) then it is false (=not returned as correct value), else true.
stat_measures = ['precision', 'specificity', 'accuracy', 'fdr', 'sensitivity', 'mcc', 'f1score', 'diagnosticOR', 'positiveLR', 'gainAcc', 'gainPpv', 'FPR']
roc_measures = ['FPR', 'sensitivity']
stat_matrix_components = ["TP", "FP", "FN", "TN"]

#######################################################################
####### FUNCTIONS
########################################################################

def runTrimmingFromDictionary(seq_dictionary, inputfasta, adapterfasta, tmpdir, geno_str, inputfastq):
	"""
	"""
	trim_dictionary = {}
	for app, program_args in application_dictionary.iteritems():
		myexec = program_args['exec_name'] # get executable name in local target machine
		# now use the prototype to reconstruct how to run the program
		if app == 'flexbar':
			for ai in program_args['ai']:
				for ao in program_args['ao']:
					for at in program_args['at']:
						# RUN
						params = 'flexbar_' + '_'.join(x for x in [str(ai), str(ao), str(at)])
						test_done_header.append(params)
						tmpfile = os.path.join(tmpdir, params)
						try:
							runthis = myexec + " --reads %s --target %s -a %s --threads 10 -ae LEFT -m 18 -q 10 --max-uncalled 8 -ai %s -ao %s -at %s" %(inputfasta, tmpfile, adapterfasta, str(ai), str(ao), str(at))
							print runthis
							retcode = call(runthis, shell=True)
							if retcode < 0:
								print >>sys.stderr, "Child was terminated by signal", -retcode
							else:
								print >>sys.stderr, "Child returned", retcode
						except OSError as e:
							print >>sys.stderr, "Execution failed:", e
						# read fastq, compute delta, acquire data
						handle = open("%s.fasta" %(tmpfile), "rU")
						record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
						handle.close()
						for res_seq, seq_vals in record_dict.iteritems():
							seq_delta = getJunctionDelta(seq_vals.seq, geno_str)
							if not trim_dictionary.has_key(res_seq):
								trim_dictionary[res_seq] = {params: seq_delta}
							else:
								trim_dictionary[res_seq][params] = seq_delta
						os.remove("%s.fasta" %(tmpfile))
		elif app == 'eautils':
			for m in program_args['m']:
				for p in program_args['p']:
					for l in program_args['l']:
						# RUN
						params = 'eautils_' + '_'.join(x for x in [str(m), str(p), str(l)])
						test_done_header.append(params)
						tmpfile = os.path.join(tmpdir, params + '.fasta')
						try:
							# fastq-mcf /opt/applications/scripts/isatk/elements/sequences/LTR.32bp.fa 32_LTR_trimming.fastq -m 30 -p 15 -l 20 -q 10 -P 33 -R -k 0 -x 0 -S | fastq_to_fasta > eautils/32_LTR_trimming.fa
							runthis = myexec + " %s %s -m %s -p %s -l %s -q 10 -P 33 -R -k 0 -x 0 -S | fastq_to_fasta > %s " %(adapterfasta, inputfastq, str(m), str(p), str(l), tmpfile)
							print runthis
							retcode = call(runthis, shell=True)
							if retcode < 0:
								print >>sys.stderr, "Child was terminated by signal", -retcode
							else:
								print >>sys.stderr, "Child returned", retcode
						except OSError as e:
							print >>sys.stderr, "Execution failed:", e
						# read fastq, compute delta, acquire data
						handle = open("%s" %(tmpfile), "rU")
						record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
						handle.close()
						for res_seq, seq_vals in record_dict.iteritems():
							seq_delta = getJunctionDelta(seq_vals.seq, geno_str)
							if not trim_dictionary.has_key(res_seq):
								trim_dictionary[res_seq] = {params: seq_delta}
							else:
								trim_dictionary[res_seq][params] = seq_delta
						os.remove("%s" %(tmpfile))
		# else:
	return trim_dictionary

def checkArgs(args):
	"""
	Check file path
	"""
	if not os.path.isfile(args.groundTruth) or not os.path.isfile(args.fastaSeqs) or not os.path.isdir(args.tmpfolder) or not os.path.isfile(args.adapterFasta) or not os.path.isfile(args.fastqSeqs):
		print "\n[AP]\tError while reading files/folders: no valid paths.\n\tExit\n"
		sys.exit()
	if args.doStats and not args.doTests:
		# the matrix file MUST be available!!!
		if not os.path.isfile(args.outMatrix):
			print "[AP]\tUsing the option --doStats without the option --doTests, the results file (the matrix file) MUST be already available and set. Check it and use the option --outMatrix to link the file in the app.\n\tExit\n"
			sys.exit()


def getJunctionDelta(fasta_str, geno_str):
	"""
	Input:
		1. fasta string
		2. genomic reference string
	Output:
		numeric delta (positive or negative {-n, 0, +n}):
		negative:: if expected trimming position is at base X whereas obtained trimmed position is at base x-k => delta = k
	Logics:
		string split and parse returned list
		- if returns a list of 2 elements, then the delta can be only or exact (thus both values are null) or negative (first element not null, second null)
		- otherwise it returns a list of 1 element. -> process this list recursively by removing the first base for each loop
	"""
	# returning value
	delta = 0 # default is 0 -> exact match
	# split fasta string with genomic string
	split_list = fasta_str.split(geno_str)
	# if list has 2 elements -> aquire values (exact or negative only)
	if len(split_list) == 2:
		# if both are empty then return the exact match (delta = 0)
		if not split_list[0] and not split_list[1]:
			delta = 0
		elif not split_list[1]:
			delta = -len(split_list[0])
	# (else) if list split has only 1 element -> positive value, of a number to be identified recursively
	elif len(split_list) == 1:
		shift = 0 # bp of shift (into the genomic string), positive delta
		while len(split_list) == 1:
			shift += 1
			split_list = fasta_str.split(str(geno_str[shift:]))
		# now we have the positive delta, do a check on list values
		if len(split_list) == 2:
			# if both are empty then return the shift match (delta = shift)
			if not split_list[0] and not split_list[1]:
				delta = shift
			else:
				print "Error, you went over the junction... roll back please."
	else:
		print "[AP]\tERROR: String split is weird:", split_list
	return delta

def saveTestResults(seq_dictionary, trimmed_seq_dict, outfile_tsv):
	"""
	Input: 
	 - trim_dictionary
	 - ground truth dictionary
	Output:
	 - CSV file of results, from the input ground truth extend dataset with the list of tests
	"""
	# first find columns and rows
	row_headers = trimmed_seq_dict.keys() # sequence headers
	column_headers = test_done_header
	# build outfile
	with open(outfile_tsv, 'w') as outf:
		writer = csv.writer(outf, delimiter = delimiter)
		writer.writerow(['sequence_header', 'full_sequence', 'ltr_sequence', 'label'] + column_headers)
		# now parse rows
		for seq_header, seq_values in trimmed_seq_dict.iteritems():
			list_to_write_prefix = [seq_header, seq_dictionary[seq_header]['full_sequence'], seq_dictionary[seq_header]['ltr_sequence'], seq_dictionary[seq_header]['label']]
			list_to_write_data = []
			for col in column_headers:
				if col in seq_values.keys():
					list_to_write_data.append(str(seq_values[col]))
				else:
					list_to_write_data.append(str(999))
			writer.writerow(list_to_write_prefix + list_to_write_data)
	return True

def saveStatMeasuresResults(delta_stat_dictionary, testLabelList, stat_measure, outfile_tsv):
	"""
	Input: 
	 - delta stat dictionary results
	 - stat measure as string [db like format]
	Output:
	 - CSV file of results, rows = delta, cols = tests
	"""
	# build outfile
	with open(outfile_tsv, 'w') as outf:
		writer = csv.writer(outf, delimiter = delimiter)
		writer.writerow(['%s_by_delta' %(stat_measure)] + testLabelList)
		for delta, stat_dictionary in delta_stat_dictionary.iteritems():
			what_to_write = [delta]
			for testLabel in testLabelList:
				what_to_write += [stat_dictionary[testLabel][stat_measure]]
			writer.writerow(what_to_write)
	return True

def saveStatResultsByTest(dictionary_byTest, header_cols, outfile_tsv):
	"""
	Input: 
	 - dictionary results by test
	 - cols headers
	Output:
	 - CSV file of results, rows = delta, cols = tests
	"""
	# build outfile
	with open(outfile_tsv, 'w') as outf:
		writer = csv.writer(outf, delimiter = delimiter)
		writer.writerow(header_cols)
		for test, vals in dictionary_byTest.iteritems():
			writer.writerow([test] + vals)
	return True

def reshapeDataByTest(testLabelList, dictionary_byDelta):
	"""
	From a dictionary of delta - values, to a dictionary of test - values.
	"""
	pair_delta_statvalue = [tuple([delta, sv]) for delta in delta_to_analyze for sv in stat_matrix_components + stat_measures] # pair of delta and stat values (that are both stat components TP.., and stat assessments Precision...)
	header_cols = ["testLabel"] + ["delta %d - %s" %(delta, sv) for delta, sv in pair_delta_statvalue]
	dictionary_byTest = {} # k = test, v = array of vals in the given order (pair_delta_statvalue)
	for test in testLabelList:
		# loop over delta and vals (dict keys)
		test_values = [] # to be next converted in numpy array based on the test index
		for delta, dkeys in pair_delta_statvalue:
			test_values.append(dictionary_byDelta[delta][test][dkeys])
		# now collect results
		dictionary_byTest[test] = test_values
	return dictionary_byTest, header_cols

def reshapeDataByTestROC(testLabelList, dictionary_byDelta):
	"""
	From a dictionary of delta - values, to a dictionary of test - values, only filtering ROC axes (FPR and TPR or sensitivity)
	"""
	pair_delta_statvalue = [tuple([delta, sv]) for delta in delta_to_analyze for sv in roc_measures] 
	header_cols = ["testLabel"] + ["delta %d - %s" %(delta, sv) for delta, sv in pair_delta_statvalue]
	dictionary_byTest = {} # k = test, v = array of vals in the given order (pair_delta_statvalue)
	for test in testLabelList:
		# loop over delta and vals (dict keys)
		test_values = [] # to be next converted in numpy array based on the test index
		for delta, dkeys in pair_delta_statvalue:
			test_values.append(dictionary_byDelta[delta][test][dkeys])
		# now collect results
		dictionary_byTest[test] = test_values
	return dictionary_byTest, header_cols

class JunctionDeltaStep:
	"""

	"""
	def __init__(self, dataframe, ):
		"""
		Init class
		"""
		# create a copy of the input matrix
		self.dataframe = dataframe
		self.innerDf = dataframe.copy()
		self.stat_dictionary = {} # k = test label, v = dict: k = {TP, FN, FP, TN} v = 4 counts -> then extended with measures keys and values
		self.deltaThreshold = None
	
	def reshape(self, deltaThreshold, groundTruthLabel):
		"""
		Input: dataframe of delta values
		Output: dataframe of binary values based on selected delta threshold
		Note: matrix values in reshape can be overlapping between source matrix and destination matrix (=> it is not a function!)
		"""
		print "[AP]\tReshape with deltaThreshold =", deltaThreshold
		# data formatting
		tmp_negative_value = -1
		# self.innerDf[self.innerDf < 0] = abs(self.innerDf) # first put all values as positive values
		self.innerDf[self.innerDf < 0] = -self.innerDf # first put all values as positive values
		self.innerDf[self.innerDf <= deltaThreshold] = tmp_negative_value # if delta is lt deltaThreshold then this IS is counted as TRUE (here first put these values to a number that will NEVER be overlapping to real data, a temp number that at the end will be reconverted)
		self.innerDf[self.innerDf > deltaThreshold] = 0 # if delta is gt deltaThreshold then this IS is counted as FALSE
		self.innerDf[self.innerDf == tmp_negative_value] = 1 # now re-change tmp values to 1
		self.innerDf[groundTruthLabel] = self.dataframe[groundTruthLabel] # now re-change tmp values to 1
		self.deltaThreshold = deltaThreshold
	
	def statMatrix(self, groundTruthLabel, testLabelList):
		"""
		Do the confusion matrix.
		Input: list of column headers of groundTruth (string) and test columns (list)
		Output: dictionary of dictionaries of the statistical measure: k = test label, v = dict: k = {TP, FN, FP, TN} v = 4 counts
		"""
		# print "Do stat matrix using groundTruth = ", groundTruthLabel, " and testing testLabelList = ", testLabelList, "\n"
		# data filtering
		for labeltest in testLabelList:
			print "[AP]\tDelta %d -> labeltest: [%s] vs groundTruthLabel: [%s]" %(self.deltaThreshold, labeltest, groundTruthLabel)
			# results with respect to expected positive values
			expectedPositive_labelResults = self.innerDf[self.innerDf[groundTruthLabel]==1][labeltest].value_counts() # series
			#print "expectedPositive_labelResults =", expectedPositive_labelResults
			TP = expectedPositive_labelResults.get(1, 0) # retrieve results 1 with expected positive values, default 0 (the second parameter of the get function)
			FN = expectedPositive_labelResults.get(0, 0)
			# results with respect to expected false values
			expectedNegative_labelResults = self.innerDf[self.innerDf[groundTruthLabel]==0][labeltest].value_counts() # series
			#print "expectedNegative_labelResults =", expectedNegative_labelResults
			FP = expectedNegative_labelResults.get(1, 0)
			TN = expectedNegative_labelResults.get(0, 0)
			# ## ---- tmp or backup start ----
			# # results with respect to expected positive values
			# expectedPositive_labelResults = self.innerDf[self.innerDf[groundTruthLabel]==1][labeltest].value_counts() # series
			# print "expectedPositive_labelResults =", expectedPositive_labelResults
			# TP = 0 # default value = 0
			# if expectedPositive_labelResults.get(1, 0) > 0: # if the series does have the key (here 1), take the value, else default 0
			# 	TP = expectedPositive_labelResults[1]
			# FN = 0 # default value = 0
			# if expectedPositive_labelResults.get(0, 0) > 0:
			# 	FN = expectedPositive_labelResults[0]
			# print "TP and FN:", TP, FN
			# # results with respect to expected false values
			# expectedNegative_labelResults = self.innerDf[self.innerDf[groundTruthLabel]==0][labeltest].value_counts() # series
			# print "expectedNegative_labelResults =", expectedNegative_labelResults
			# FP = 0 # default value = 0
			# if expectedNegative_labelResults.get(1, 0) > 0: # if the series does have the key (here 1), take the value, else default 0
			# 	print "label 1 exists"
			# 	FP = expectedNegative_labelResults[1]
			# TN = 0 # default value = 0
			# if expectedNegative_labelResults.get(0, 0) > 0: # if the series does have the key (here 1), take the value, else default 0
			# 	TN = expectedNegative_labelResults[0]
			# ## ---- tmp end ----
			# fill local dictionary
			self.stat_dictionary[labeltest] = {
					'TP': TP,
					'FN': FN,
					'FP': FP,
					'TN': TN,
				}
	
	def computeStatMeasures(self, testLabelList, ):
		"""
		Input: the list of labels for which compute stat assessments.
		Output: update dictionary of results by label self.stat_dictionary
		"""
		# print "Compute StatMeasures"
		for labeltest in testLabelList:
			smr = StatMatrixResults(self.stat_dictionary[labeltest]['TP'], self.stat_dictionary[labeltest]['FP'], self.stat_dictionary[labeltest]['FN'], self.stat_dictionary[labeltest]['TN']) # init object of stat results
			self.stat_dictionary[labeltest]['precision'] = smr.precision()
			self.stat_dictionary[labeltest]['sensitivity'] = smr.sensitivity()
			self.stat_dictionary[labeltest]['accuracy'] = smr.accuracy()
			self.stat_dictionary[labeltest]['specificity'] = smr.specificity()
			self.stat_dictionary[labeltest]['fdr'] = smr.fdr()
			self.stat_dictionary[labeltest]['mcc'] = smr.mcc()
			self.stat_dictionary[labeltest]['gainAcc'] = smr.gainAcc()
			self.stat_dictionary[labeltest]['gainPpv'] = smr.gainPpv()
			self.stat_dictionary[labeltest]['positiveLR'] = smr.positiveLR()
			self.stat_dictionary[labeltest]['negativeLR'] = smr.negativeLR()
			self.stat_dictionary[labeltest]['diagnosticOR'] = smr.diagnosticOR()
			self.stat_dictionary[labeltest]['f1score'] = smr.f1score()
			self.stat_dictionary[labeltest]['FPR'] = smr.FPR()
			del smr
	
	def free(self, ):
		"""
		Free memory from copy of dataframe.
		"""
		self.innerDf = None

class StatMatrixResults:
	"""
	Do the statistical measures from stat matrix.
	"""
	def __init__(self, TP, FP, FN, TN):
		"""
		Init class
		"""
		self.TP = TP
		self.FN = FN
		self.TN = TN
		self.FP = FP
	def precision (self, ):
		"""
		TP / (TP + FP)
		"""
		return float(self.TP) / (self.TP + self.FP)
	def sensitivity (self, ):
		"""
		TP / (TP + FN)
		"""
		return float(self.TP) / (self.TP + self.FN)
	def specificity (self, ):
		"""
		TN / (FP + TN)
		"""
		return float(self.TN) / (self.FP + self.TN)
	def accuracy (self, ):
		"""
		TP + TN / Total population
		"""
		return ( float(self.TP) + float(self.TN) )/ (self.FN + self.TP + self.FP + self.TN)
	def fdr (self, ):
		"""
		FP / (TP + FP)
		"""
		return float(self.FP) / (self.TP + self.FP)
	def mcc (self, ):
		"""
		Mattew correlation coefficient, related to chi2 (http://en.wikipedia.org/wiki/Matthews_correlation_coefficient)
		( (TP*TN)-(FP*FN) ) / math.sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) )
		"""
		return float( (self.TP*self.TN)-(self.FP*self.FN) ) / math.sqrt( (self.TP+self.FP)*(self.TP+self.FN)*(self.TN+self.FP)*(self.TN+self.FN) )
	def gainAcc (self, ):
		"""
		Gain in accuracy http://en.wikipedia.org/wiki/Gain_%28information_retrieval%29
		"""
		total = float(self.FN + self.TP + self.FP + self.TN)
		acc = ( float(self.TP) + float(self.TN) )/ (self.FN + self.TP + self.FP + self.TN)
		r_pos = ( self.TP + self.FN ) / total
		r_neg = ( self.TN + self.FP ) / total
		r = math.pow(r_pos, 2) + math.pow(r_neg, 2)
		gain = acc/float(r)
		return gain
	def gainPpv (self, ):
		"""
		Gain in precision http://en.wikipedia.org/wiki/Gain_%28information_retrieval%29
		"""
		total = float(self.FN + self.TP + self.FP + self.TN)
		ppv = float(self.TP) / (self.TP + self.FP)
		r_pos = ( self.TP + self.FN ) / total
		gain = ppv/float(r_pos)
		return gain
	def positiveLR (self, ):
		"""
		Positive likelihood ratio
		"""
		TPR = float(self.TP) / (self.TP + self.FN)
		FPR = float(self.FP) / (self.FP + self.TN)
		return TPR/FPR
	def negativeLR (self, ):
		"""
		Negative likelihood ratio
		"""
		return (1-(float(self.TP) / (self.TP + self.FN)))/ (float(self.TN) / (self.FP + self.TN))
	def diagnosticOR (self, ):
		"""
		Diagnostic odds ration
		"""
		TPR = float(self.TP) / (self.TP + self.FN)
		FPR = float(self.FP) / (self.FP + self.TN)
		lr_pos = TPR/FPR
		lr_neg = (1-(float(self.TP) / (self.TP + self.FN)))/ (float(self.TN) / (self.FP + self.TN))
		return lr_pos / lr_neg
	def f1score (self, ):
		"""
		F1-score
		"""
		return float(2*self.TP)/(2*self.TP + self.FP + self.FN)
	def FPR (self, ):
		"""
		False positive rate: 1 − specificity
		"""
		return 1 - self.specificity()

def plotROC(file_dataframe, fileoutput, zoombox_shape = [-0.0, 0.1, 0.9, 1.0]):
	"""
	Input:
	 - file_dataframe: file written with all statistical resuts
	 fileoutput: pdf file
	 - zoombox_shape= array of xy values to zoom (zoombox_shape[xmin, xmax, ymin, ymax] = [-0.04, 0.14, 0.86, 1.04])
	"""
	# some predefined parameters
	plt.rcParams['font.size'] = 14
	#plt.rcParams['font.family'] = 'Tahoma'
	plt.rcParams['axes.labelsize'] = 14
	plt.rcParams['xtick.labelsize'] = 14
	plt.rcParams['ytick.labelsize'] = 14
	plt.rcParams['legend.fontsize'] = 'medium'
	plt.rcParams['legend.fancybox'] = True
	# plt.rcParams['legend.numpoints'] = 1
	plt.rcParams['legend.shadow'] = True
	plt.rcParams['figure.figsize'] = 16,9
	plt.rcParams['figure.dpi'] = 150
	#
	pp = PdfPages(fileoutput)
	df = pd.read_csv(file_dataframe, sep=delimiter, header=0, index_col=[0]) 
	for index in delta_to_analyze:
		x='delta %d - FPR' %(index)
		y='delta %d - sensitivity' %(index)
		color = my_cnames[index%len(my_cnames)]
		plt.scatter(df[x], df[y], c=color, marker=my_markers[index%len(my_markers)], alpha=0.5, label='Delta %d' %(index), s=80)
	##
	plt.plot([0, 1], [0, 1], ls="--", c='red', label="Random", linewidth=3)
	plt.xlim(xmin=-0.04, xmax=1.04)
	plt.ylim(ymin=-0.04, ymax=1.04)
	plt.title("ROC curve - Trimming optimization" )
	plt.xlabel('FPR (1 - Specificity)')
	plt.ylabel('TPR (Sensitivity)')
	plt.grid(True)
	# # plt.grid()
	plt.legend(loc='lower right');
	# plt.show()
	plt.savefig(pp, format='pdf', dpi=150, facecolor='w', edgecolor='w', orientation='landscape', papertype='a3', transparent=True, bbox_inches=None, pad_inches=0.1)
	plt.close()
	
	# details of best solutions
	graphical_threshold = 0.002
	for index in delta_to_analyze:
		x='delta %d - FPR' %(index)
		y='delta %d - sensitivity' %(index)
		color = my_cnames[index%len(my_cnames)]
		plt.scatter(df[x], df[y], c=color, marker=my_markers[index%len(my_markers)], alpha=0.5, label='Delta %d' %(index), s=80)
	
	## plt.xticks(0, 1, 0.05)
	xmin = zoombox_shape[0] - graphical_threshold
	xmax = zoombox_shape[1] + graphical_threshold
	ymin = zoombox_shape[2] - graphical_threshold
	ymax = zoombox_shape[3] + graphical_threshold
	plt.xlim(xmin=xmin, xmax=xmax)
	plt.ylim(ymin=ymin, ymax=ymax)
	plt.title("ROC curve - Trimming optimization" )
	plt.xlabel('FPR (1 - Specificity)')
	plt.ylabel('TPR (Sensitivity)')
	plt.grid(True)
	# # plt.grid()
	# plt.legend(loc='upper left');
	plt.legend(loc='upper right');
	# plt.show()
	plt.savefig(pp, format='pdf', dpi=150, facecolor='w', edgecolor='w', orientation='landscape', papertype='a3', transparent=True, bbox_inches=None)
	plt.close()

	# details of best solutions up to the forth delta 
	for index in delta_to_analyze[:4]:
		x='delta %d - FPR' %(index)
		y='delta %d - sensitivity' %(index)
		color = my_cnames[index%len(my_cnames)]
		plt.scatter(df[x], df[y], c=color, marker=my_markers[index%len(my_markers)], alpha=0.5, label='Delta %d' %(index), s=100)
	## plt.xticks(0, 1, 0.05)
	plt.xlim(xmin=-0.001, xmax=0.021)
	plt.ylim(ymin=0.979, ymax=1.001)
	plt.title("ROC curve - Trimming optimization" )
	plt.xlabel('FPR (1 - Specificity)')
	plt.ylabel('TPR (Sensitivity)')
	plt.grid(True)
	# # plt.grid()
	# plt.legend(loc='upper left');
	plt.legend(loc='upper right');
	# plt.show()
	plt.savefig(pp, format='pdf', dpi=150, facecolor='w', edgecolor='w', orientation='landscape', papertype='a3', transparent=True, bbox_inches=None)
	plt.close()

	## 
	pp.close()
		
# def findBestROCZoom(file_dataframe):
# 	"""
# 	Input:
# 	 - file_dataframe: file name to be processed with pandas
# 	Output:
# 	 - array of zoombox_shape as [xmin, xmax, ymin, ymax]
# 	"""
# 	df = pd.read_csv(file_dataframe, sep=delimiter, header=0, index_col=[0]) 
# 	df[df[x1]<=0.1][df[y1]>=0.8].index




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
	# strftime("%Y%m%d %H:%M:%S", gmtime())
	checkArgs(args)
	### --------- DO TESTS COMBINATORIAL --------------
	if args.doTests:
		# acquire ground truth file as CSV, tab delimited
		seq_gt_dict = {} # sequence ground truth dictionary: k=header, v=dict:: k=col name (full_sequence, ltr_sequence, label), v=value
		with open(args.groundTruth, 'r') as inf:
			reader = csv.reader(inf, delimiter = delimiter)
			for row in reader:
				if not row[0].startswith("#"): # skip comments with line starting #
					seq_gt_dict[row[0]] = {
						'full_sequence': row[1], 
						'ltr_sequence': row[2], 
						'label': row[3],
						}
		# do trimming with all available software
		trimmed_seq_dict = runTrimmingFromDictionary(seq_gt_dict, args.fastaSeqs, args.adapterFasta, args.tmpfolder, args.genomicSeq, args.fastqSeqs)
		# save results in a CSV file
		saveTestResults(seq_gt_dict, trimmed_seq_dict, args.outMatrix)
	
	### ------------------ STATS ----------------------
	if args.doStats:
		# compute statistical assessments
		# create a dataframe from written file
		df = pd.read_csv(args.outMatrix, sep=delimiter, header=0, index_col=[0,1,2]) # acquire df as dataframe, considering the first row as header and the first 3 columns (0,1,2) as row keys, where these fiels have been defined here above as full_sequence, ltr_sequence, label (STATIC!!!) 
		gtLabel = df.columns[0] # keep the columns ground truth, the first one
		testLabels = [x for x in df.columns[1:]] # keep the columns to test, just exclude the column label
		print "[AP]\tReference (ground truth) and comparing columns (tests):\n", gtLabel, "\n", ' - '.join(str(x) for x in testLabels)
		# NOTE: rows -> df.index; cols -> df.columns [without row name, data only]
		results_dictionary_byDelta = {} # k = delta, v = JunctionDeltaStep.stat_dictionary
		for delta in delta_to_analyze:
			print "[AP]\tDelta:", delta
			current_time = time.time() # get read of starting time
			print "[AP]\tObject instantiation"
			jds = JunctionDeltaStep(df) # acquire df and create tmp copy
			print "[AP]\tReshape DataFrame"
			jds.reshape(delta, gtLabel) # reshape data
			print "[AP]\tCompute statistical matrixes"
			jds.statMatrix(gtLabel, testLabels) # now you have a dictionary of stat matrix for each test label in jds.stat_dictionary
			print "[AP]\tComputing statistical assessments"
			jds.computeStatMeasures(testLabels)
			print "[AP]\tElapsed time for this delta at %d has been %.2f [seconds]\n[AP]" %(delta, time.time() - current_time)
			results_dictionary_byDelta[delta] = jds.stat_dictionary
			jds = None
		# save global results in different files
		print "[AP]\tWrite results of the statistical measures."
		for sm in stat_measures:
			print "[AP]\t", sm
			saveStatMeasuresResults(results_dictionary_byDelta, testLabels, sm, os.path.splitext(args.outMatrix)[0] + ".%s.csv" %(sm) )
		# build a new easier data structure
		results_dictionary_byTest, column_headers = reshapeDataByTest(testLabels, results_dictionary_byDelta)
		# write output of the last data structure
		completeResultFile = os.path.splitext(args.outMatrix)[0] + ".allDataByTest.tsv"
		print "[AP]\tSave statistical results in a file (.allDataByTest.tsv):", completeResultFile
		saveStatResultsByTest(results_dictionary_byTest, column_headers, completeResultFile)
		## ROC
		# build a new easier data structure for ROC
		print "[AP]\tBuild a ROC dataset"
		ROC_results_dictionary_byTest, ROC_column_headers = reshapeDataByTestROC(testLabels, results_dictionary_byDelta)
		# write output of the last data structure
		ROC_data_outfile = os.path.splitext(args.outMatrix)[0] + ".ROC.tsv"
		print "[AP]\tSave ROC file:", ROC_data_outfile
		saveStatResultsByTest(ROC_results_dictionary_byTest, ROC_column_headers, ROC_data_outfile)
		# plot results
		ROC_plot_outfile = os.path.splitext(args.outMatrix)[0] + ".ROC.pdf"
		print "[AP]\tPlot ROC data in", ROC_plot_outfile
		plotROC(ROC_data_outfile, ROC_plot_outfile)

	elapsed_time = time.time() - start_time # get read of finishing time
	print "\n[AP]\tTask Finished, closing.\n\tElapsed time: %.2f [seconds]\n" %(elapsed_time)


# sentinel
if __name__ == "__main__":
    main()
