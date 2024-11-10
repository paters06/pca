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
# Command found in https://note.nkmk.me/en/python-script-file-path/
dir1 = os.path.dirname(os.path.abspath(__file__))
# Insert .. command to go to the upper directory
dir2 = dir1 + '/..'
# Setting the package directory path for the modules execution
sys.path.append(dir2)
#######################################################################

# Local project
import src.nurbs as rbs
import src.plotting.plottingScripts as plts
import src.spline_functions.surfaceRefinements as srfn


# Complete list of control points
controlPointsPatch1 = np.array([[0.0,0.0,0.0],[0.2,0,0],[0.0,0.6,0],[0.2,0.4,0.0]])
weightsPatch1 = np.ones((controlPointsPatch1.shape[0],1))

controlPointsPatch2 = np.array([[0.2,0.4,0.0],[0.6,0.4,0.0],[0.0,0.6,0.0],[0.6,0.6,0.0]])
weightsPatch2 = np.ones((controlPointsPatch2.shape[0],1))

multiP = [controlPointsPatch1,controlPointsPatch2]
multiw = [weightsPatch1,weightsPatch2]

multip = [1,1]
multiq = [1,1]

multiU = [np.array([0,0,1,1]),np.array([1,1,2,2])]
multiV = [np.array([0,0,1,1]),np.array([0,0,1,1])]

surface = rbs.MultiPatchNURBSSurface(multiU,multiV,multip,multiq,\
                                     multiP,multiw)

localRefinement = 'Y'
patchesToRefine = [0,1]
numreflist = [2,1]
reflist = [['h'],['h']]
dirlist = [['U','V'],['U','V']]

if localRefinement == 'Y':
    srfn.localPatchRefinement(surface,patchesToRefine,numreflist,reflist,dirlist)

fullcpts = surface.createMultipatchSurface()
fullcpu,fullcpv = surface.createMultipatchTangentSurface()

# surface.plotMultipatchSurface()
# surface.plotMultipatchTangentSurface("u")
