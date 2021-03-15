#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:26:33 2020

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

# Rectangle with central circular hole
L = 1.0
H = 0.2

# Simple rectangle
#P = np.array([[0,0,0],
#                  [0.5*L,0,0],
#                  [L,0,0],
#                  [0,0.5*H,0],
#                  [0.5*L,0.5*H,0],
#                  [L,0.5*H,0],
#                  [0,H,0],
#                  [0.5*L,H,0],
#                  [L,H,0]])

#w = np.array([[1],[1],[1],[1],[1],[1],[1],[1],[1]])

#p = 2
#q = 2

#U = np.array([0,0,0,1,1,1])
#V = np.array([0,0,0,1,1,1])

# Annular disk
#P = np.array([[1,0,4],
#              [1,1,4],
#              [0,1,4],
#              [-1,1,4],
#              [-1,0,4],
#              [-1,-1,4],
#              [0,-1,4],
#              [1,-1,4],
#              [1,0,4],
#              [2,0,4],
#              [2,2,4],
#              [0,2,4],
#              [-2,2,4],
#              [-2,0,4],
#              [-2,-2,4],
#              [0,-2,4],
#              [2,-2,4],
#              [2,0,4]])

#w = np.array([[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1],
#              [1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1],[1/np.sqrt(2)],[1]])

#p = 2
#q = 1

#U = np.array([0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1])
#V = np.array([0,0,1,1])

# Square with central circular hole (eight part)
#P = np.array([[1,0,4],
#              [1,1,4],
#              [0,1,4],
#              [-1,1,4],
#              [-1,0,4],
#              [-1,-1,4],
#              [0,-1,4],
#              [1,-1,4],
#              [1,0,4],
#              [2,0,4],
#              [2,2,4],
#              [0,2,4],
#              [-2,2,4],
#              [-2,0,4],
#              [-2,-2,4],
#              [0,-2,4],
#              [2,-2,4],
#              [2,0,4]])

Ra = 1.0
Rb = 4.0
# ka = cos(pi/2 - pi/8)
ka = 0.5*np.sqrt(2 - np.sqrt(2))
# kb = cos(pi/4)
kb = 0.5*np.sqrt(2)
# kw = cos(pi/8)
kw = 0.5*np.sqrt(2 + np.sqrt(2))

P = np.array([[0,Ra,4],
              [Ra*ka,Ra,4],
              [Ra*kb,Ra*kb,4],
              [0,Rb,4],
              [Rb*ka,Rb,4],
              [Rb,Rb,4]])

#P = np.array([[Ra*kb,Ra*kb,4],
#              [Ra,Ra*ka,4],
#              [Ra,0,4],
#              [Rb,Rb,4],
#              [Rb,Rb*ka,4],
#              [Rb,0,4]])

w = np.array([[1],[kw],[1],[1],[kw],[1]])

p = 2
q = 1

U = np.array([0,0,0,1,1,1])
V = np.array([0,0,1,1])


surface1 = rbs.NURBSSurface(U,V,p,q,P,w)
cpts = surface1.createSurface()
cpu,cpv = surface1.createTangentSurface()
#surface1.plotSurface()
surface1.plotTangentSurface("v")
