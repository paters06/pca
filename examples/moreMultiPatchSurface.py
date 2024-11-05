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

Ra = 1.0
Rb = 4.0
# ka = cos(pi/2 - pi/8)
ka = 0.5*np.sqrt(2 - np.sqrt(2))
# kb = cos(pi/4)
kb = 0.5*np.sqrt(2)
# kw = cos(pi/8)
kw = 0.5*np.sqrt(2 + np.sqrt(2))

# Complete list of control points
fullP = np.array([[0,Ra,0],[Ra*ka,Ra,0],[Ra*kb,Ra*kb,0],
                  [0,Rb,0],[Rb*ka,Rb,0],[Rb,Rb,0],
                  [Ra,Ra*ka,0],[Ra,0,0],[Rb,Rb*ka,0],
                  [Rb,0,0]])

fullw = np.array([[1],[kw],[1],[1],[kw],[1],[kw],[1],[kw],[1]])

idcontrolpoints = [[0,1,2,3,4,5],[2,6,7,5,8,9]]
localcontrolpoints = [fullP[idcontrolpoints[0],:],fullP[idcontrolpoints[1],:]]

multip = [2,2]
multiq = [1,1]

multiU = [np.array([0,0,0,1,1,1]),np.array([0,0,0,1,1,1])]
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

localRefinement = 'N'
patchesToRefine = [0,1]
reflist = [['h','h','h','h'],['h','h']]
dirlist = [['U','V','U','V'],['U','V']]

if localRefinement == 'Y':
    multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints = \
    srfn.localPatchRefinement(patchesToRefine,reflist,dirlist,multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints)


#Computing the surface
fullcpts = rbs.MultiPatchNURBSSurface(multiU,multiV,multip,multiq,fullP,fullw)
# fullcpu,fullcpv = rbs.MultipatchNurbsSurfaceTangent(multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints)
fullcpts.plotMultipatchSurface()
# plts.plotMultipatchSurface(fullcpts,localcontrolpoints)
#plts.plotmultipatchTangentSurface(fullcpts,fullcpu)
#plts.plotmultipatchTangentSurface(fullcpts,fullcpv)
