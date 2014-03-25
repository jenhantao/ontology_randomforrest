import sys
from numpy import *
import matplotlib.pyplot as plt

# plots interaction likelihoods for terms VS term to term distance scores for random forrest
# inputs: path to likelihood array file, path distance array file

likelihoodArray = load(sys.argv[1])
likelihoodArray = likelihoodArray[likelihoodArray.files[0]]

distanceArray = load(sys.argv[2])
distanceArray = distanceArray[distanceArray.files[0]]

# reshape arrays
likelihoodArray = asarray(reshape(likelihoodArray, (1,likelihoodArray.size))[0])
distanceArray = asarray(reshape(distanceArray, (1,distanceArray.size))[0])
if True:
	distanceArray = -log(asarray(reshape(distanceArray, (1,distanceArray.size))[0]))
	distanceArray[isinf(distanceArray)] = 0
	print distanceArray.min()
	print distanceArray.max()

heatmap, xedges, yedges = histogram2d(likelihoodArray, distanceArray, bins=4)
plt.figure()
# grey scale
plt.imshow(heatmap)
axis = asarray(plt.axis())
axis[0]=0
#axis[2]=0
plt.axis(axis)
plt.set_cmap('gray')
plt.xlabel("Interaction Likelihood Score")
plt.ylabel("Forrest Distance Score")

plt.savefig(sys.argv[3]+'.png')
