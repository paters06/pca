#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:26:33 2020

@author: hernando
@collaborator: paternina
"""

# Python libraries
import numpy as np

#######################################################################
# DO NOT REMOVE THIS SEGMENT
import os
import sys

# Get current working directory
dir1 = os.getcwd()
# Insert .. command to go to the upper directory
dir2 = dir1 + '/..'
# Change directory
os.chdir(dir2)
# Get the new current working directory
dir3 = os.getcwd()
# Setting the package directory path for the modules execution
sys.path.append(dir3)
#######################################################################

# Local project
import src.nurbs as rbs
import src.plottingScripts as plts
import src.surfaceRefinements as srfn


# Complete list of control points
fullP = np.array([[0,0,0],[0.2,0,0],[0,0.6,0],[0.2,0.4,0],
                  [0.6,0.4,0],[0.6,0.6,0]])

fullw = np.array([[1],[1],[1],[1],[1],[1]])

idcontrolpoints = [[0,1,2,3],[3,4,2,5]]
localcontrolpoints = [fullP[idcontrolpoints[0],:],fullP[idcontrolpoints[1],:]]

multip = [1,1]
multiq = [1,1]

multiU = [np.array([0,0,1,1]),np.array([1,1,2,2])]
multiV = [np.array([0,0,1,1]),np.array([0,0,1,1])]

# First Patch
#P = np.array([[0,0,0],
#              [0.2,0,0],
#              [0,0.6,0],
#              [0.2,0.4,0]])

#w = np.array([[1],[1],[1],[1]])

#p = 1
#q = 1

#U = np.array([0,0,1,1])
#V = np.array([0,0,1,1])

# Second Patch
#P = np.array([[0.2,0.4,0],
#              [0.6,0.4,0],
#              [0,0.6,0],
#              [0.6,0.6,0]])

#w = np.array([[1],[1],[1],[1]])

#p = 1
#q = 1

#U = np.array([0,0,1,1])
#V = np.array([0,0,1,1])

localRefinement = 'Y'
patchesToRefine = [0,1]
reflist = [['h','h','h','h'],['h','h']]
dirlist = [['U','V','U','V'],['U','V']]

if localRefinement == 'Y':
    multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints = \
    srfn.localPatchRefinement(patchesToRefine,reflist,dirlist,multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints)


#Computing the surface
fullcpts = rbs.multipatchNurbsSurface(multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints)
fullcpu,fullcpv = rbs.multipatchNurbsSurfaceTangent(multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints)
plts.plotmultipatchSurface(fullcpts,localcontrolpoints)
#plts.plotmultipatchTangentSurface(fullcpts,fullcpu)
#plts.plotmultipatchTangentSurface(fullcpts,fullcpv)
