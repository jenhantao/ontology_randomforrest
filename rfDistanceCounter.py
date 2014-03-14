# reads in any number of random forrests, oncotypes, and gene ontology terms to characterize the random forrest
# computes the following: 

# inputs: path to trees, gene ontology,  outputDirectory

import sys
from igraph import *
from os import listdir
from os.path import isfile, join
from numpy import *

with open(sys.argv[2]) as f:
        ontologyData = f.readlines()
# map each child parent term to a numeric value
termIndexHash = {}
terms = []
for line in ontologyData:
        tokens = line.strip().split("\t")
        child = tokens[0]
        parent = tokens[1]
	terms.append(child)
	terms.append(parent)
terms = sorted(list(set(terms)))
for i in range(len(terms)):
	termIndexHash[terms[i]] = i
save(sys.argv[1]+"terms",terms)
	
treeSize = []
levelHash = {} # key - gene, value - array of levels
pairDistanceHash = {} # key - gene1_gene2, value array of distances
cooccurenceHash = {} # key - gene1_gene2, value array of 1s and 0s indicating cooccurrence in a given tree
# read in trees
pathToTrees = sys.argv[1] # path to directory containing tree files

forrest = []
treeFilePaths= [ f for f in listdir(pathToTrees) if isfile(join(pathToTrees,f)) ]
for treePath in treeFilePaths:
	if not "npy" in treePath:
		treeFile = open(pathToTrees + "/" +treePath)
		treeData = treeFile.readlines()
		treeFile.close()
		currentTree = Graph(directed=True)
		edges = [] # array of 2 member tuples denoting 
		splitVar = []
		splitPoint = []
		status = [] # 1 if a node has children, -1 otherwise
		prediction = [] # NA for non terminal nodes, Neg/Pos otherwise
		vertexNumber = 0
		for line in treeData[1:]: # skip first line which is the header
			# operations needed to build tree
			tokens = line.strip().split()
			if not tokens[4] == "-1":
				leftEdge = (vertexNumber, int(tokens[0])-1) # edge leading from the current vertex to the left child
				rightEdge = (vertexNumber, int(tokens[1])-1) # edge leading from the current vertex to the right child
				edges.append(leftEdge)
				edges.append(rightEdge)
			if not tokens[2] == "NA":
				splitVar.append(tokens[2][:tokens[2].index("_",3)])
			else:
				splitVar.append(tokens[2])
			splitPoint.append(tokens[3])
			status.append(int(tokens[4]))
			prediction.append(tokens[5])
			vertexNumber += 1
			# operations to build other objects
			levelHash[tokens[2]] = []	
		currentTree.add_vertices(vertexNumber)
		currentTree.add_edges(edges)
		currentTree.vs["splitVar"] = splitVar
		currentTree.vs["status"] = status
		currentTree.vs["prediction"] = prediction
		forrest.append(currentTree)

# traverse tree and compute properties 
distanceArray = zeros((len(terms),len(terms)))
for i in range(len(forrest)):
	tree = forrest[i]
	# compute pairwise distance
	distances = tree.shortest_paths()
	# compute co-ocurrence
	for i in range(len(tree.vs)-1):
		for j in range(i+1,len(tree.vs)):
			if not tree.vs[i]["splitVar"] == "NA" and not tree.vs[j]["splitVar"] == "NA": 
				index1 = termIndexHash[tree.vs[i]["splitVar"]]
				index2 = termIndexHash[tree.vs[j]["splitVar"]]
				if not str(distances[i][j]) == "inf":
					distanceArray[index1][index2] += distances[i][j]
					distanceArray[index2][index1] += distances[i][j]
for i in range(len(terms)):
	for j in range(len(terms)):
		if distanceArray[i][j] > 0.0:
			distanceArray[i][j] = float(len(forrest))/float(distanceArray[i][j])
save(pathToTrees+"treeDistanceArray",distanceArray)
	
	









