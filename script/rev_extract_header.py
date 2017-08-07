#!/usr/bin/python
# -*- coding: utf-8 -*-
# da usare quando vuoi estrarre linee da un file con headers presenti in un altro.

import sys 
# acquire ID from file as list
ids, flag = set(x.strip().replace("\n","") for x in open(sys.argv[1])), False
#print ids
for l in sys.stdin:
	flag = (l.split('\t')[0].strip().replace("\n","") not in ids)
	if flag: print l,
