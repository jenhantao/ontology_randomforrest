# reads in npy file and converts to to npz
import sys
from numpy import *

array = load(sys.argv[1])
savez_compressed(sys.argv[1].replace(".npy",""), array)
