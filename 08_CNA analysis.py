# Load the packages
import os
import dynetan as dna
from dynetan.toolkit import *
from dynetan.viz import *
from dynetan.proctraj import *
from dynetan.gencor import *
from dynetan.contact import *

import warnings

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from statistics import mean, stdev
from matplotlib import colors
from IPython import get_ipython

dnap = DNAproc(notebookMode=False)

workDir = "./"
psfFile = os.path.join(workDir, "cas.psf")
dcdFiles = os.path.join(workDir, "cas.dcd") 

pathToData = "./casResults/"
fileNameRoot = "dnaData_Cas"

segIDs = ["SYS"]
customResNodes = {}

numWinds = 1 
cutoffDist = 4.5
contactPersistence = 0.75

dnap.setNumWinds(numWinds)
dnap.setCutoffDist(cutoffDist)
dnap.setContactPersistence(contactPersistence)
dnap.setSegIDs(segIDs)
dnap.setNodeGroups(customResNodes)
fullPathRoot = os.path.join(pathToData, fileNameRoot)

print(f"Loading topology file {psfFile} and trajectory file(s) {dcdFiles}.")
dnap.loadSystem(psfFile, dcdFiles)
print("System loaded.")
print("MDAnalysis universe:", dnap.getU().trajectory)

dnap.checkSystem()
dnap.prepareNetwork(verbose=True)

print("Aligning trajectory...")
dnap.alignTraj(verbose=2)

print("Finding contacts...")
dnap.findContacts(stride=10, verbose=2, n_cores=8)

print("Recursion limit:", sys.getrecursionlimit())
sys.setrecursionlimit(3000)
print("New recursion limit:", sys.getrecursionlimit())

print("Filtering contacts...")
dnap.filterContacts(notSameRes=True, notConsecutiveRes=False, removeIsolatedNodes=True)

dnap.calcCor(ncores=4, verbose=2)

numSampledFrames = 50000
dnap.setNumSampledFrames(numSampledFrames)
dnap.calcCartesian(backend="openmp")

dnap.calcGraphInfo()

print("Graph with {} nodes and {} edges".format(len(dnap.nxGraphs[0].nodes),len(dnap.nxGraphs[0].edges)))

print("Calculating optimal paths...")
dnap.calcOptPaths(ncores=4)

print("Calculating edge betweeness...")
dnap.calcBetween(ncores=4)

dnap.calcEigenCentral()

dnap.calcCommunities()

dnap.saveData(fullPathRoot)

dcdstride = 10

print("We will save {} heavy atoms and {} frames.".format(dnap.workU.atoms.n_atoms,len(dnap.workU.trajectory[::dcdstride]) ))

warnings.filterwarnings("ignore", category=UserWarning)

dnap.saveReducedTraj(fullPathRoot, stride = dcdstride)

print("\n\n END OF Calculate \n\n")



























