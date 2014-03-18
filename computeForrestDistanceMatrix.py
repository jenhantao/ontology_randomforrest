# computes the "average" average distance matrix for all term pairs in a random forrest
# avg dist = (d1+d2+...+dk)/(i1+i2+...+ik) where d is the pairwise distance between two terms. i is 1 if d is nonzero and 0 otherwise.

# inputs: path to tree distance matrices directory

import sys
from igraph import *
from os import listdir
from os.path import isfile, join
from numpy import *

treeFilePaths= [ f for f in listdir(pathToTrees) if isfile(join(pathToTrees,f)) ]
