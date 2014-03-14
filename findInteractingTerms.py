# reads in gene ontology terms and true labels to find interesting terms 

# inputs: gene ontology - parent to child relations, gene to term mapping,true labels,gene pairs tested

import sys
from igraph import *
from numpy import *

# build ontology tree
#with open('/Users/Admin/Documents/ontology_randomforrest/GO_biological_process.term_term.txt') as f:
with open(sys.argv[1]) as f:
	ontologyData = f.readlines()
# map each child parent term to a numeric value
vertexCount = 0
edges = []
edgeRelations = []
termIndexHash = {}
for line in ontologyData:
	tokens = line.strip().split("\t")
	child = tokens[0]
	parent = tokens[1]
	relation = tokens[2]

	parentIndex = -1
	if parent in termIndexHash:
		parentIndex = termIndexHash[parent]		
	else:
		parentIndex = vertexCount
		vertexCount += 1
		termIndexHash[parent] = parentIndex

	childIndex = -1
	if child in termIndexHash:
		childIndex = termIndexHash[child]		
	else:
		childIndex = vertexCount
		vertexCount += 1
		termIndexHash[child] = childIndex
	edge = (parentIndex, childIndex)
	edges.append(edge)
	edgeRelations.append(relation)

del ontologyData
ontologyTree = Graph(edges,directed=True)
ontologyTree.es["relation"] = edgeRelations
# find root node
degrees = ontologyTree.degree(type="in")
root = degrees.index(0)

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
del geneTermData
# read in tested gene pairs and Labels and build genepairLabelHash
genepairLabelHash = {} # key tested gene pair in the form of a tuple, value: label - 1 if there is a negative interaction and 0 otherwise
#with open('/Users/Admin/Documents/ontology_randomforrest/costanzo/true_labels.txt') as f:
with open(sys.argv[3]) as f:
	trueLabelData = f.readlines()
labels = []
for line in trueLabelData:
	labels.append(int(line.strip()))
del trueLabelData
#with open('/Users/Admin/Documents/ontology_randomforrest/costanzo/gene_sets.txt') as f:
with open(sys.argv[4]) as f:
	genePairData = f.readlines()
for i in range(len(genePairData)):
	tokens = genePairData[i].strip().split("$")
	genePair = tuple(sorted(tokens)) # gene pairs are always sorted lexographically
	genepairLabelHash[genePair] = labels[i]
del genePairData

# calculate background rate
backgroundRate = float(sum(labels))/float(len(labels)) # exactly 0.1 <- this is kind of suspicious

# find interesting gene pairs
terms = sorted(termGeneArrayHash.keys()) # term pairs are always examined in lexographical order
likelihoodArray = zeros(shape = (len(terms),len(terms)))
for i in range(len(terms)-1):
	for j in range(i+1,len(terms)):
		term1 = terms[i]
		term2 = terms[j]	
		# compute sets for each term, the set of genes that appear under term1, but not term2 as well as the set of genes that appear under term2, but not term1
		geneSet1 = set(termGeneArrayHash[term1])
		geneSet2 = set(termGeneArrayHash[term2])
		intersection = geneSet1.intersection(geneSet2)
		geneSet1 = sorted(list(geneSet1 - intersection))
		geneSet2 = sorted(list(geneSet2 - intersection))
		numTestedPairs = 0
		numInteractingPairs = 0
		for gene1 in geneSet1:
			for gene2 in geneSet2:
				currentGenePair = (gene1,gene2)
				if currentGenePair in genepairLabelHash:
					numTestedPairs += 1
					if genepairLabelHash[currentGenePair] == 1:
						numInteractingPairs += 1
		if numInteractingPairs > 0:
			rate = float(numInteractingPairs)/float(numTestedPairs)
		else:
			rate = 0.0
		likelihoodEnrichment = rate/backgroundRate
#		if likelihoodEnrichment > 1.0:
#			print term1, term2, likelihoodEnrichment
		likelihoodArray[i][j] = rate
		likelihoodArray[j][i] = rate
		
from tempfile import TemporaryFile
save("likelihoodArray",likelihoodArray)












