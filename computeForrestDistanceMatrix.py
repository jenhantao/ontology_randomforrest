# computes the "average" average distance matrix for all term pairs in a random forrest
# avg dist = (d1+d2+...+dk)/(i1+i2+...+ik) where d is the pairwise distance between two terms. i is 1 if d is nonzero and 0 otherwise.

# inputs: path to tree distance matrices directory

import sys
from igraph import *
from os import listdir
from os.path import isfile, join
from numpy import *

path = sys.argv[1]
outName = sys.argv[1].split("/")[-2].strip()+"_forrestArray"
filePaths= [ f for f in listdir(path) if isfile(join(path,f)) and "npz" in f and not "term" in f]
forrestDistanceArray = None # stores the sum of all distances seen at every position
nonZeroPositionsArray = None # stores the number of ocurrences of nonzero distance values at each position


# sum all the arrays array together

counter=0
for fp in filePaths:
	counter +=1
	currentArray = load(path+fp)
	currentArray = currentArray[currentArray.files[0]]
	currentNonZeros= (currentArray > 0) * 1
	print counter
	if forrestDistanceArray ==None:
		forrestDistanceArray = currentArray
		nonZeroPositionArray = currentNonZeros
	else:
		forrestDistanceArray += currentArray
		nonZeroPositionArray += currentNonZeros

# normalize arrays 
print forrestDistanceArray
print nonZeroPositionArray
outputArray = nonZeroPositionArray/forrestDistanceArray
outputArray[isnan(outputArray)] =0
print outputArray
savez_compressed(outName, outputArray)
