# reads in any number of random forrests, oncotypes, and gene ontology terms to characterize the random forrest
# computes the following: 

# inputs: path to tree, gene ontology

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

tokens = sys.argv[1].split("/")
tokens[-1] = "terms"

savez_compressed("/".join(tokens),terms)
	
treeSize = []
# read in tree
treePath = sys.argv[1] # path to tree
treeFile = open(treePath)
treeData = treeFile.readlines()
treeFile.close()
tree = Graph(directed=True)
edges = [] # array of 2 member tuples denoting 
splitVar = []
splitPoint = []
status = [] # 1 if a node has children, -1 otherwise
prediction = [] # NA for non terminal nodes, Neg/Pos otherwise
vertexNumber = 0
# create graph from tree data
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
tree.add_vertices(vertexNumber)
tree.add_edges(edges)
tree.vs["splitVar"] = splitVar
tree.vs["status"] = status
tree.vs["prediction"] = prediction

# traverse tree and compute properties 
distanceArray = zeros((len(terms),len(terms)))
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
savez_compressed(treePath.replace(".txt",""), distanceArray)
	
	









