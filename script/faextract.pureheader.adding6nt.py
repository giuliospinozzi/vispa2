#!/usr/bin/python
# -*- coding: utf-8 -*-
# da usare quando vuoi estrarre FASTQ avendo i dati dal BED.

######## WARNING!!!!! ##########
### FASTA SEQUENCE WITH NEW LINES WILL CAUSE PROBLEMS!!! ######
### NO NEW LINES ARE ALLOWED!!! #######

import sys, random
# acquire ID from file as list, adding the prefix @
ids, lineno, flag = set('>' + x.strip() for x in open(sys.argv[1])), 0, False
#print ids
for l in sys.stdin:
	lineno += 1
	if lineno%2 == 1: 
		flag = (l.split(' ')[0].strip() in ids)
		if flag: 
			header = l.strip()
	if lineno%2 == 0: 
		added = random.sample(['A', 'C', 'G', 'T']*10,6)
		if flag: 
			#print added, l, flag
			print "%s\n%s%s" %(header, ''.join(added), l.strip())

