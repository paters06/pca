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

#P = np.array([[0,0,0],[3,0,3],[6,0,3],[9,0,0],[0,2,2],[3,2,5],[6,2,5],[9,2,2],
#               [0,4,0],[3,4,3],[6,4,3],[9,4,0]])

#w = np.array([[1],[3],[3],[1],[1],[5],[5],[1],[1],[4],[4],[1]])

#P = np.array([[1,0,0],[1,np.sqrt(2)-1,0],[np.sqrt(2)-1,1,0],[0,1,0],[2.5,0,0],[2.5,0.75,0],
#              [0.75,2.5,0],[0,2.5,0],[4,0,0],[4,4,0],[4,4,0],[0,4,0]])

#w = np.array([[1],[0.5*(1+(1/np.sqrt(2)))],[0.5*(1+(1/np.sqrt(2)))],[1],[1],[1],
#              [1],[1],[1],[1],[1],[1]])

#p = 2
#q = 2

#U = np.array([0,0,0,0.5,1,1,1])
#V = np.array([0,0,0,1,1,1])

P = np.array([[0,0,0],
              [0.2,0,0],
              [0,0.6,0],
              [0.2,0.4,0]])

w = np.array([[1],[1],[1],[1]])

p = 1
q = 1

U = np.array([0,0,1,1])
V = np.array([0,0,1,1])

#cpts = rbs.nurbsSurface(U,V,p,q,P,w)
cpts = rbs.nurbsSurface(U,V,p,q,P,w)
#cpu,cpv = rbs.nurbsSurfaceTangent(U,V,p,q,P,w)
plts.plottingSurface(cpts[0,:,:],cpts[1,:,:],cpts[2,:,:])
#plts.plotTangentSurface(cpts[0,:,:],cpts[1,:,:],cpts[2,:,:],cpu[0,:,:],cpu[1,:,:],cpu[2,:,:])
#plts.plotTangentSurface(cpts[0,:,:],cpts[1,:,:],cpts[2,:,:],cpv[0,:,:],cpv[1,:,:],cpv[2,:,:])
