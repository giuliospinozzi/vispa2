#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys, csv
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
import random
from operator import itemgetter, attrgetter

"""

"""
header = """


"""

#########################################################################################
####### GLOBAL VARS
#########################################################################################

nucleotides = list("ACGTACGTACGTACGT") # uniform distribution of the 4 bp (equal probability) + only 1 N (1/4 on groups of nucleotides, that corresponds to 1/17 so far)
nucleotides_set = set(nucleotides)
genomic_seq = "AGAGATCAAGTCTCACTATGTTGCCCAGGTTGGTCTCGAACTCTTGGGCTCAAGCAATCCTCTCACCTCAGCCTCCCAAAGTGCTGGGATTACAGACATGAGCCACCCTGCTCGGCTAGAATT"
ltr_63 = "TGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA"
ltr_32 = "ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA"
sequence_to_scrumble = ltr_32
sequence_to_append = genomic_seq
max_search_random_sequences = 200 # 200
outfile_prefix = "test_outfile_v4"
max_mutation_span = 4 # number of max consecutive nucleotides to mutate
max_insertion_span = 4 # number of max consecutive nucleotides to insert
changes_max_numb_in_sequence = 5 # number of max random mutations in the sequence
changes_max_span = 8 # number of maximum span per random mutation change
changes_types = ['mut', 'ins', 'del'] # DO NOT CHANGE THESE NAMES!!! THESE NAMES ARE ALSO LINKED IN THE DICTIONARY OF CLASSIFICATION RULES 
changes_max_random_trials = 2000 # 20000
groundtruth_file_header = ['header', 'full sequence', 'ltr sequence', 'label'] # header file, sorted fields: what and how to write in the ground truth file 

# concetto di fondo: usare 3 criteri di regole qui sotto, tali per cui queste sono il CONFINE  per cui riportare una etichetta come buona o no (quindi <= per ogni numero mentre per le stringhe e' ad occorrenza, ovvero ne basta solo uno di evento per etichettare quella come da scartare): 
# 1) sequenza in se' dove si posizionano le X nei punti dove NON si devono mai applicare variations; 
# 2) max number of allowed variations
# 3) max len for eah variation (if applicable)
## Mutation section, rules
label_mut_denidedbases = "ACCCTTTTAGTCAGTGTGGAAAATCTCTXXXX" # given the same LTR input sequence, recreate the sequence composed by Y (or the input base) and X, where Y are accepted and X are denided bases.
label_mut_maxnumb = 9 # labeling definition, variation type: mutation, max number of variations along the overall sequence. rule number 2
label_mut_maxspan = 3 # rule number 3 for mutation
## Insertions
label_ins_denidedbases = "ACCCTTTTAGTCAGTGTGGAAAATCTCTAGXX" # given the same LTR input sequence, recreate the sequence composed by Y (or the input base) and X, where Y are accepted and X are denided bases.
label_ins_maxnumb = 2 # labeling definition, variation type: mutation, max number of variations along the overall sequence. rule number 2
label_ins_maxspan = 4 # rule number 3 for mutation
## Deletions
label_del_denidedbases = "ACCCTTTTAGTCAGTGTGGAAAATCTCTXXCA" # given the same LTR input sequence, recreate the sequence composed by Y (or the input base) and X, where Y are accepted and X are denided bases.
label_del_maxnumb = 3 # labeling definition, variation type: mutation, max number of variations along the overall sequence. rule number 2
label_del_maxspan = 15 # rule number 3 for mutation

# now collect the rules in a dictionary
classification_rules = {
	'mut': {
		'denidedbases': [bp_index for bp_index in range(0,len(label_mut_denidedbases)) if list(label_mut_denidedbases)[bp_index]=='X'], # 0-based array (look at the range interval to confirm this statement) of bp where this variation type MUST BE avoided
		'maxnumb': label_mut_maxnumb,
		'maxspan': label_mut_maxspan,
		},
	'ins': {
		'denidedbases': [bp_index for bp_index in range(0,len(label_ins_denidedbases)) if list(label_ins_denidedbases)[bp_index]=='X'], # 0-based array (look at the range interval to confirm this statement) of bp where this variation type MUST BE avoided
		'maxnumb': label_ins_maxnumb,
		'maxspan': label_ins_maxspan,
		},
	'del': {
		'denidedbases': [bp_index for bp_index in range(0,len(label_del_denidedbases)) if list(label_del_denidedbases)[bp_index]=='X'], # 0-based array (look at the range interval to confirm this statement) of bp where this variation type MUST BE avoided
		'maxnumb': label_del_maxnumb,
		'maxspan': label_del_maxspan,
		},
	}


#########################################################################################
### MY FUNCTIONS
#########################################################################################

def doDeletion(input_string, del_type, del_parameters):
	"""
	deletions: <del_type> | <del_parameters>:
		from_last | len >0
		from_first | len >0
		Nbp | [N size {1..seqlen}, starting position 0-based >0 && <len(input_string)]
	"""
	instring_list = list(input_string)
	output_string = ""
	if del_type == "from_last":
		output_string = ''.join( instring_list[:len(instring_list)-del_parameters] )
	elif del_type == "from_first":
		output_string = ''.join( instring_list[del_parameters:] )
	elif del_type == "Nbp":
		output_string = ''.join( instring_list[:del_parameters[1]] + instring_list[del_parameters[1]+del_parameters[0]:] )
	else:
		print "[AP]\tError, deletion type not defined."
		sys.exit()
	return output_string

def doInsertion(input_string, ins_type, ins_parameters, random_seq = True):
	"""
	if random(ACGTN) do:
		insertions: <ins_type> | <ins_parameters>:
			from_last | len >0 
			from_first | len >0
			Nbp | [N size {1..seqlen}, starting position 0-based >0 && <len(input_string)]
	if defined dstring do:
		insertions: <ins_type> | <ins_parameters>:
			from_Nbp (not used because so far it is the only one) | [user defined string, starting position 0-based >0 && <len(input_string)]
	"""
	instring_list = list(input_string)
	output_string = ""
	if random_seq:
		if ins_type == "from_last":
			output_string = instring_list + ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters))
		elif ins_type == "from_first":
			output_string = ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters)) + instring_list
		elif ins_type == "Nbp":
			output_string = ''.join( instring_list[:ins_parameters[1]] ) + ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters[0])) + ''.join(instring_list[ins_parameters[1]:] )
		else:
			print "[AP]\tError, insertion type not defined."
			sys.exit()
	else: # if not random
		output_string = ''.join( instring_list[:ins_parameters[1]] ) + ins_parameters[0] + ''.join(instring_list[ins_parameters[1]:])
	return output_string

def doMutation(input_string, mut_type, mut_parameters):
	"""
	(random) mutation excluding same string to replace: <mut_type> | <mut_parameters>:
		Nbp | [start bp, end bp 0-based] where this is a closed interval thus you are mutating from the staring the the ending bp included.

	NB: each base in the interval MUST be different, else you could obtain an unexpected mutation in 2 positions in the same span: from AAAA, span 2 from start, you may obtain CCAA but also ACAA -> the first one is ok but the second one corresponds to the case of SINGLE mutation at base 2.
	"""
	instring_list = list(input_string)
	output_string = ''.join( instring_list[:mut_parameters[0]] )
	if mut_type == "Nbp":
		if mut_parameters[1]-mut_parameters[0]<=0:
			print "[AP]\tWarning: your mutated string has an input interval NON positive! mut_parameters[1] - mut_parameters[0] <= 0::", mut_parameters[1], mut_parameters[0]
		else:
			for bp_index in range(mut_parameters[0],mut_parameters[1]):
				char_tomutate = str(instring_list[bp_index])
				random_string = str(random.choice(nucleotides))
				# random_is_different_fromsource = False
				# print "Before:: ", char_tomutate, random_string
				while char_tomutate == random_string:
					random_string = str(random.choice(nucleotides))
				# print "\tAfter random is ", random_string
				output_string += random_string
			output_string += ''.join(instring_list[mut_parameters[1]:])
	else:
		print "[AP]\tError, mutation type not defined."
		sys.exit()
	return output_string

def classifyRead(classification_rules, variation_type, variation_starting_bp0b, variation_span, variation_count):
	"""
	Input:
	 - classification rules
	 - read(s)
	Output:
	 - for the single input read and variation type, return 0 (not to keep) or 1 (to keep), that is the classification label
	"""
	read_classification_label = 1 # this will be 0 or 1
	# for each variation type apply a specific classification use case decision rule
	# #### --- versione sempre funzionante ------
	# # variated bp in the region of the denied bp?
	# for variated_bp in range(variation_starting_bp0b, variation_starting_bp0b+variation_span):
	# 	if variated_bp in classification_rules[variation_type]['denidedbases']:
	# 		read_classification_label = 0
	#### --- versione piu elegante
	variated_bp_not_allowed = [ variated_bp for variated_bp in range(variation_starting_bp0b, variation_starting_bp0b+variation_span) if variated_bp in classification_rules[variation_type]['denidedbases'] ] # this array contains OLNY overlapped reads with at least ONE base of the NOT ALLOWED array (deniedbp)
	if len(variated_bp_not_allowed) > 0:
		read_classification_label = 0
	if variation_span > classification_rules[variation_type]['maxspan']:
		read_classification_label = 0
	if variation_count > classification_rules[variation_type]['maxnumb']:
		read_classification_label = 0
	## return label
	return read_classification_label


#########################################################################################
### MAIN
#########################################################################################

def main():
	dict_sequences = {} # k=header, v=string of sequence + details (tsv)
	dict_labeled_sequences = {} # k=header, v=dict:: k={'full sequence', 'ltr sequence', 'label'}, v = values where label in {0,1} as label for classification (0=not IS, 1=IS)
	## header composition: <global_type>:<local_type>:<paramteres>(-<parameters>)*:<len>
	# -------- DELETIONS -------------
	# delete N bp {1..11} from last
	for k in range(1,16):
		header = "del:from_last:%dbp" %(k)
		scrumbled_seq = doDeletion(sequence_to_scrumble, "from_last", k)
		details = "\t%s\tdel\tfrom_last\t%d\tltrlen\t%d" %(scrumbled_seq, k, len(scrumbled_seq))
		dict_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = scrumbled_seq + sequence_to_append + details
		dict_labeled_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = {
			'full sequence': scrumbled_seq + sequence_to_append.strip(), 
			'ltr sequence': scrumbled_seq, 
			'label': classifyRead(classification_rules, 'del', len(sequence_to_scrumble) - k, k, 1),
			} # add label as classification rules
	# delete N bp {1..16} from first
	for k in range(1,16):
		header = "del:from_first:%dbp" %(k)
		scrumbled_seq = doDeletion(sequence_to_scrumble, "from_first", k)
		details = "\t%s\tdel\tfrom_first\t%d\tltrlen\t%d" %(scrumbled_seq, k, len(scrumbled_seq))
		dict_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = scrumbled_seq + sequence_to_append + details
		dict_labeled_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = {
			'full sequence': scrumbled_seq + sequence_to_append.strip(), 
			'ltr sequence': scrumbled_seq, 
			'label': classifyRead(classification_rules, 'del', 0, k, 1),
			} # add label as classification rules
	# delete N bp within the sequence
	for bpsize in range(1,20): # length of consecutive bp deletion
		for starting_position in range(1,len(sequence_to_scrumble)-1): # starting position [0-based]
			# do the operation only if you will not delete a number of non LTR bases to avoid duplicated events (deletion of 5bp of a seq with len 32 from bp 30 and 28 will produce the same results of removing 4 bp)
			if len(sequence_to_scrumble) - starting_position - bpsize > 0:
				header = "del:Nbp:len_%dbp-startingat_%d" %(bpsize, starting_position+1)
				scrumbled_seq = doDeletion(sequence_to_scrumble, "Nbp", [bpsize, starting_position])
				details = "\t%s\tdel\tNbp\tlen\t%d\tstartingat\t%d\tltrlen\t%d" %(scrumbled_seq, bpsize, starting_position+1, len(scrumbled_seq))
				dict_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = scrumbled_seq + sequence_to_append + details
				dict_labeled_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = {
					'full sequence': scrumbled_seq + sequence_to_append.strip(), 
					'ltr sequence': scrumbled_seq, 
					'label': classifyRead(classification_rules, 'del', starting_position, bpsize, 1),
					} # add label as classification rules
	
	# -------- INSERTIONS -------------
	# generate a dict of sets of random bases -> this dataset will be used as predetermined sequences
	frombp = 20
	rand_dict = {} # k=len, v=array of seqs
	for seqlen in range(1,max_insertion_span): # max range sequences 10 from 1
		seq_set = set()
		for k in range(0,max_search_random_sequences): # on a wide set of cases, change this INT to change choice search
			seq_string = ''.join(str(random.choice(nucleotides)) for x in range(0,seqlen)) 
			seq_set.add(seq_string)
		rand_dict[seqlen] = list(seq_set)
	# insert all sequences predefined of N bp within the sequence in sliding positions
	for bpsize, bpseq_list in rand_dict.iteritems(): # bp insertions
		for bpseq in bpseq_list:
			for starting_position in range(frombp,len(sequence_to_scrumble)-1): # starting position [0-based] # avoid simulating insertions as perfect attachment to the sequence ## OLD BACKUP::for starting_position in range(len(sequence_to_scrumble)-1-frombp,len(sequence_to_scrumble)+1):
				header = "ins:from_Nbp:len_%dbp-%s-startingat_%d" %(bpsize, bpseq, starting_position+1)
				scrumbled_seq = doInsertion(sequence_to_scrumble, "from_Nbp", [bpseq, starting_position], random_seq=False)
				details = "\t%s\tins\tfrom_Nbp\tlen\t%d\tinserted_seq\t%s\tstartingat\t%d\tltrlen\t%d" %(scrumbled_seq, bpsize, bpseq, starting_position+1, len(scrumbled_seq))
				dict_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = scrumbled_seq + sequence_to_append + details
				dict_labeled_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = {
					'full sequence': scrumbled_seq + sequence_to_append.strip(), 
					'ltr sequence': scrumbled_seq, 
					'label': classifyRead(classification_rules, 'ins', starting_position, bpsize, 1),
					} # add label as classification rules

	# --------- MUTATIONS --------------
	# mutations random
	frombp = 8
	for ending_position in range(frombp+max_mutation_span,len(sequence_to_scrumble)-1): # ending position [0-based] ## OLD BACKUP::for ending_position in range(len(sequence_to_scrumble)-max_mutation_span-frombp,len(sequence_to_scrumble)):
		for starting_position in range(frombp,len(sequence_to_scrumble)): # starting position [0-based]
			if (ending_position - starting_position <= max_mutation_span) and (ending_position > starting_position):
				header = "mut:Nbp:starting_%dbp-ending_%d-span_%dbp" %(starting_position+1, ending_position, ending_position-starting_position )
				scrumbled_seq = doMutation(sequence_to_scrumble, "Nbp", [starting_position, ending_position])
				details = "\t%s\tmut\tNbp\tstarting\t%d\tending\t%d\tspan\t%d\tltrlen\t%d" %(scrumbled_seq, starting_position+1, ending_position, ending_position-starting_position, len(scrumbled_seq) )
				dict_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = scrumbled_seq + sequence_to_append + details
				dict_labeled_sequences[header + ":ltrlen_%d" %(len(scrumbled_seq))] = {
					'full sequence': scrumbled_seq + sequence_to_append.strip(), 
					'ltr sequence': scrumbled_seq, 
					'label': classifyRead(classification_rules, 'mut', starting_position, abs(ending_position-starting_position), 1),
					} # add label as classification rules

	# --------- N VARIATIONS (MUT, DEL, INS) RANDOM -------------
	# do replicates of each test
	for replica in range(0,changes_max_random_trials): # do N tests (changes_max_random_trials)
		# print "----- replica ", replica, "--------\nvar_type, var_nbp, var_startingbase, shift_by_var_nbp, reference_var_startingbase"
		# first acquire the variations to do
		random_variations = [] # list of sets, composed by the triplet <mutation type>, <number of bp span>, <starting base of the reference sequence>
		seq_positions = list(xrange(3,len(sequence_to_scrumble)-1))
		for num_changes in range(2,random.randint(3,changes_max_numb_in_sequence)): # vary the number of changes in a range between 2 and max (2 because 1 will be one of the base cases mut, del or ins)
			random_base = random.choice(seq_positions)
			random_variations.append( (random.choice(changes_types), random.randint(1,changes_max_span), random_base ) ) # <mutation type>, <number of bp span>, <starting base of the reference sequence>. NB: it is prevented/avoided that a single base is randomply selected for 2 mutations
			seq_positions.remove(random_base)
		# sort the variation array -> so that variations will be applied sorted in silico
		random_variations_sortedbyposition = sorted(random_variations, key=itemgetter(2, 1), reverse = True)
		# print random_variations_sortedbyposition
		# now given all the acquired variations, apply them all to the original sequence (doing a variation per time)
		sequence_variant_content = 'N'*len(sequence_to_scrumble) # sequence representation content in a string of N with variations with X
		header = "var:"
		details = "\t"
		shift_by_var_nbp = 0 # given the mutation type and span, for subsequent mutations you have to do a shift in the reference number
		scrumbled_seq = sequence_to_scrumble
		variation_count = 0 # how many variations are we considering? if > 2 then write the read, else discard it. this option is required for internal checks (see next nested if in the variation types)
		deletion_count = 0 # questa var contiene le delezioni!!!
		mutation_count = 0 # nuber of mutations
		insertion_count = 0 # number of insertions
		classification_label_multiple_variation_array = [] # this array will contain a label for each variation. once you have a complete array, than just multiplying values will return you the final binary label (without so far accounting for the max number of mutations allowed)
		for (var_type, var_nbp, var_startingbase) in random_variations_sortedbyposition:
			# reference_var_startingbase = var_startingbase-shift_by_var_nbp
			reference_var_startingbase = var_startingbase
			#print var_type, var_nbp, var_startingbase+1, reference_var_startingbase
			if var_type == 'mut': # this variation type does not contribute to any sequence len modification -> no shift in starting bp
				if reference_var_startingbase+var_nbp < len(scrumbled_seq):
					scrumbled_seq = doMutation(scrumbled_seq, "Nbp", [reference_var_startingbase, reference_var_startingbase+var_nbp])
					# print scrumbled_seq
					header += "%s:Nbp:span_%dbp-startingat_%d|" %(var_type, var_nbp, var_startingbase+1)
					details += "%s\tNbp\tspan\t%d\tstartingat\t%d\t" %(var_type, var_nbp, var_startingbase+1)
					sequence_variant_content = sequence_variant_content[0:reference_var_startingbase] + 'X'*var_nbp + sequence_variant_content[reference_var_startingbase+var_nbp:]
					variation_count += 1 # if you arrived here then this variation is done => increment the relative count
					classification_label_multiple_variation_array.append(classifyRead(classification_rules, 'mut', reference_var_startingbase, var_nbp, 1)) # add for this mutation the corresponding putative classification label
					mutation_count += 1
			elif var_type == 'del':
				if reference_var_startingbase >= 5:
					# print scrumbled_seq, var_type, var_nbp, var_startingbase+1, reference_var_startingbase
					scrumbled_seq = doDeletion(scrumbled_seq, "Nbp", [var_nbp, reference_var_startingbase])
					# print scrumbled_seq
					header += "%s:Nbp:span_%dbp-startingat_%d|" %(var_type, var_nbp, var_startingbase+1)
					details += "%s\tNbp\tspan\t%d\tstartingat\t%d\t" %(var_type, var_nbp, var_startingbase+1)
					deletion_count += sequence_variant_content[reference_var_startingbase:reference_var_startingbase+var_nbp].count('N')
					sequence_variant_content = sequence_variant_content[0:reference_var_startingbase] + sequence_variant_content[reference_var_startingbase+var_nbp:]
					variation_count += 1 # if you arrived here then this variation is done => increment the relative count
					# print sequence_variant_content, variation_count
					classification_label_multiple_variation_array.append(classifyRead(classification_rules, 'del', reference_var_startingbase, var_nbp, 1)) # add for this mutation the corresponding putative classification label
					deletion_count += 1
			elif var_type == 'ins':
				if reference_var_startingbase <= len(scrumbled_seq)-1:
					scrumbled_seq = doInsertion(scrumbled_seq, "Nbp", [var_nbp, reference_var_startingbase], random_seq=True)
					# print scrumbled_seq
					header += "%s:Nbp:span_%dbp-startingat_%d|" %(var_type, var_nbp, var_startingbase+1)
					details += "%s\tNbp\tspan\t%d\tstartingat\t%d\t" %(var_type, var_nbp, var_startingbase+1)
					sequence_variant_content = sequence_variant_content[0:reference_var_startingbase] + 'X'*var_nbp + sequence_variant_content[reference_var_startingbase:]
					variation_count += 1 # if you arrived here then this variation is done => increment the relative count
					classification_label_multiple_variation_array.append(classifyRead(classification_rules, 'ins', reference_var_startingbase, var_nbp, 1)) # add for this mutation the corresponding putative classification label
					insertion_count += 1
			else:
				print "[AP]\tError on variation type!!!!!! ->", var_type, " not recognized, exiting...\n"
				sys.exit()
		
		# evaluate the labels for classification
		sequence_classification_label = 1
		# 1. how many mutations/ins/dels did you have? are they overcoming the max allowed values?
		for number_observed_variations, number_accepted_variations in [(mutation_count, classification_rules['mut']['maxnumb']), (deletion_count, classification_rules['del']['maxnumb']), (insertion_count, classification_rules['ins']['maxnumb'])]:
			if number_observed_variations > number_accepted_variations:
				sequence_classification_label *= 0
		
		# report read details
		if variation_count > 1:
			dict_sequences[header.rstrip('|') + ":ltrlen_%d" %(len(scrumbled_seq))] = scrumbled_seq + sequence_to_append + "\t%s\tvar\tvariation_count_noDel\t%d\tvariation_count_withDel\t%d" %(scrumbled_seq, sequence_variant_content.count('X'), sequence_variant_content.count('X') + deletion_count) + details + "\tltrlen\t%d" %(len(scrumbled_seq))
			# print scrumbled_seq
			# 2. the final classification label derives from the multiplication of all labels in classification_label_multiple_variation_array
			sequence_classification_label *= reduce(lambda x, y: x*y, classification_label_multiple_variation_array)
			dict_labeled_sequences[header.rstrip('|') + ":ltrlen_%d" %(len(scrumbled_seq))] = {
					'full sequence': scrumbled_seq + sequence_to_append.strip(), 
					'ltr sequence': scrumbled_seq, 
					'label': sequence_classification_label,
					} # add label as classification rules

	# --------- WRITE OUTPUT ------------
	# 1. file of sequences with details
	outf = open(outfile_prefix + ".csv", 'w')
	outf.write("header\tauto label\tfull sequence\tltr sequence\tvariation type\tdetails\n")
	for k, v in dict_sequences.iteritems():
		outf.write(k + "\t%s\t" %(dict_labeled_sequences[k]['label']) + v + "\n")
	outf.close()
	# 2. file of sequences with classification labels (file format as acquired by trimming optimization)
	# outf = open(outfile_prefix + ".labeled.groundtruth.csv", 'w')
	with open(outfile_prefix + ".labeled.groundtruth.csv", 'wb') as csvfile:
		writer = csv.writer(csvfile, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
		writer.writerow(["#" + groundtruth_file_header[0]] + groundtruth_file_header[1:])
		for k, v in dict_labeled_sequences.iteritems():
			what_to_write = [k] + [ v[i] for i in groundtruth_file_header[1:] ] 
			writer.writerow(what_to_write)
	# 3. write fasta file of sequences
	outf = open(outfile_prefix + ".fa", 'w')
	for k, v in dict_labeled_sequences.iteritems():
		outf.write(">%s\n%s\n" %(k, v['full sequence']))
	outf.close()
	# 4. write fastq file of sequences
	outf = open(outfile_prefix + ".fq", 'w')
	for k, v in dict_labeled_sequences.iteritems():
		outf.write("@%s\n%s\n+\n%s\n" %(k, v['full sequence'], str('9'*len(v['full sequence'])) ))
	outf.close()

# sentinel
if __name__ == "__main__":
    main()







