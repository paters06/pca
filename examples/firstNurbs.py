#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 18:16:56 2020
Modified on Sun Aug 02 17:59:20 2020

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

P = np.array([[0,0],[1,1],[3,2],[4,1],[5,-1]])
U = np.array([0,0,0,1,2,3,3,3])
w = np.array([[1],[1],[1],[1],[1]])
p = 2

cpts = rbs.nurbsCurve(U,p,P,w)
cppts = rbs.nurbsCurveTangent(U,p,P,w)
cppts_norm = np.sqrt(cppts[:,0]**2 + cppts[:,1]**2)
cppts_norm = np.reshape(cppts_norm,(len(cppts_norm),1))
cppts_unit = cppts/cppts_norm

#print(np.hstack((cpts,cppts)))

#plts.plotCurve2d(cpts,P,"no")
plts.plotTangentCurve2d(cpts,cppts_unit,P,"no")
#plts.plotTangentCurve2d(cpts,cppts,P,"no")
