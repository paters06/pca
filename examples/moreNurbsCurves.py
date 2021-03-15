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

# Local project
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

import src.nurbs as rbs
import src.plottingScripts as plts


#P = np.array([[0,0],[1,1],[3,2],[4,1],[5,-1]])
#U = np.array([0,0,0,1,2,3,3,3])
#w = np.array([[1],[1],[1],[1],[1]])
#p = 2

#P = np.array([[1,0],
#              [1,1],
#              [0,1],
#              [-1,1],
#              [-1,0],
#              [-1,-1],
#              [0,-1],
#              [1,-1],
#              [1,0]])

#w = np.array([[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1]])

#U = np.array([0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1])
#p = 2

P = np.array([[1,0],
              [1,1],
              [0,1]])

w = np.array([[1],[1/np.sqrt(2)],[1]])

U = np.array([0,0,0,1,1,1])
p = 2

curve1 = rbs.NURBSCurve(U,p,P,w)
cpts = curve1.createCurve()
cppts = curve1.createTangentCurve()
#curve1.plotCurve()
curve1.plotTangentCurve()
