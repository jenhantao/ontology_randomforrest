# reads in any number of random forrests, oncotypes, and gene ontology terms to characterize the random forrest
# computes the following: 

# inputs: path to trees, oncotypes, predicted labels, true labels, gene terms, gene ontology

import sys
from igraph import *
from os import listdir
from os.path import isfile, join

# read in trees
pathToTrees = sys.argv[1] # path to directory containing tree files
forrest = []
treeFilePaths= [ f for f in listdir(pathToTrees) if isfile(join(pathToTrees,f)) ]
for treePath in treeFilePaths:
	treeFile = open(pathToTrees + "/" +treePath)
	treeData = treeFile.readlines()
	treeFile.close()
	currentTree = Graph()
	edges = [] # array of 2 member tuples denoting 
	splitVar = []
	splitPoint = []
	status = [] # 1 if a node has children, -1 otherwise
	prediction = [] # NA for non terminal nodes, Neg/Pos otherwise
	vertexNumber = 0
	for line in treeData[1:]: # skip first line which is the header
		tokens = line.strip().split()
		if not tokens[4] == "-1":
			leftEdge = (vertexNumber, int(tokens[0])-1) # edge leading from the current vertex to the left child
			rightEdge = (vertexNumber, int(tokens[1])-1) # edge leading from the current vertex to the right child
			edges.append(leftEdge)
			edges.append(rightEdge)
		splitVar.append(tokens[2])
		splitPoint.append(tokens[3])
		status.append(int(tokens[4]))
		prediction.append(tokens[5])
		vertexNumber += 1
	currentTree.add_vertices(vertexNumber)
	currentTree.add_edges(edges)
	currentTree.vs["splitVar"] = splitVar
	currentTree.vs["status"] = status
	currentTree.vs["prediction"] = prediction
	forrest.append(currentTree)


for tree in forrest:
	print(tree)
