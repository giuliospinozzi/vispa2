#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys 
# acquire ID from file as list, adding the prefix @
ids, lineno, flag = set('@' + x for x in open(sys.argv[1])), 0, False
for l in sys.stdin:
    lineno += 1
    if lineno%4 == 1: flag = (l in ids)
    if flag: print l,

