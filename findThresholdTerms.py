# reads in gene ontology terms and true labels to find interesting terms 

# inputs: gene ontology - parent to child relations, gene to term mapping,true labels,gene pairs tested

import sys
from numpy import *

# read in gene term relations to build term gene sets
#with open('/Users/Admin/Documents/ontology_randomforrest/GO_biological_process.gene_term_fixed.txt') as f:
with open(sys.argv[2]) as f:
	geneTermData = f.readlines()
termGeneArrayHash = {} # key: term, value: array of genes corresponding to term
for line in geneTermData:
	tokens = line.strip().split()	
	term = tokens[1]
	gene = tokens[0]
	if term in termGeneArrayHash:
		termGeneArrayHash[term].append(gene)
	else:
		termGeneArrayHash[term] = [gene]

termGeneNumHash = {}
for key in termGeneArrayHash.keys()











