# reads in any number of random forrests, oncotypes, and gene ontology terms to characterize the random forrest
# computes the following: 

# inputs: path to trees, oncotypes, predicted labels, true labels, gene terms, gene ontology

import sys
from igraph import *
from os import listdir
from os.path import isfile, join

treeSize = []
levelHash = {} # key - gene, value - array of levels
pairDistanceHash = {} # key - gene1_gene2, value array of distances
cooccurenceHash = {} # key - gene1_gene2, value array of 1s and 0s indicating cooccurrence in a given tree
# read in trees
pathToTrees = sys.argv[1] # path to directory containing tree files
outputDir = sys.argv[2]
forrest = []
treeFilePaths= [ f for f in listdir(pathToTrees) if isfile(join(pathToTrees,f)) ]
for treePath in treeFilePaths:
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

# assign levels
for tree in forrest:
	queue = [0]
	tree.vs[0]["level"]=1
	while queue:
		current = queue[0]
		queue = queue[1:]
		edges = tree.es.select(_source_in=[current])
		for edge in edges:
			tree.vs[edge.tuple[1]]["level"] = tree.vs[current]["level"] + 1
			queue.append(edge.tuple[1])
# traverse tree and compute properties 
for tree in forrest:
	# operations for traversal
	queue = [0]
	tree.vs[0]["level"]=1
	# calculations on current tree
	# compute size of tree
	treeSize.append(len(tree.vs))	
	# compute pairwise distance
	distances = tree.shortest_paths()

	# compute co-ocurrence
	for i in range(len(tree.vs)-1):
		for j in range(i+1,len(tree.vs)):
			if not tree.vs[i]["splitVar"] == "NA" and not tree.vs[j]["splitVar"] == "NA": 
				key = None
				key1 = tree.vs[i]["splitVar"] + "/" + tree.vs[j]["splitVar"]
				key2 = tree.vs[j]["splitVar"] + "/" + tree.vs[i]["splitVar"]
	#			key1 = str(i)+ "/" + str(j)
	#			key2 = str(j)+ "/" + str(i)
				if key1 in cooccurenceHash:
					key = key1
				if key2 in cooccurenceHash:
					key = key2
				if key == None:
					key = key1
					cooccurenceHash[key] = []
				if not str(distances[i][j]) == "inf" or not str(distances[j][i]) == "inf": # threshold by distance?
					cooccurenceHash[key].append(1)
				else:
					cooccurenceHash[key].append(0)
	while queue:
		current = queue[0]
		queue = queue[1:]
		# calculations on current vertex
		levelHash[tree.vs[current]["splitVar"]].append(tree.vs[current]["level"])
		edges = tree.es.select(_source_in=[current])
		for edge in edges:
			queue.append(edge.tuple[1])
avgLevelFile = open(outputDir+"/averageLevel.tsv", "w")
avgLevelFile.write("gene\taverage level\n")
for gene in levelHash.keys():
	avgLevelFile.write( gene+ "\t" + str(mean(levelHash[gene]))+"\n")
cooccurenceFile = open(outputDir+"/cooccurence.tsv","w")
cooccurenceFile.write("gene1\tgene2\taverage cooccurence\n")
for key in cooccurenceHash:
	tokens = key.split("/")
	cooccurenceFile.write(tokens[0] + "\t" + tokens[1] + "\t" + str(mean(cooccurenceHash[key])) + "\n")
cooccurenceFile.close()
#	print tree.vs["status"]
#	print tree.vs["prediction"]
#	print tree.vs.select(status=-1)["prediction"]
avgLevelFile.close()
