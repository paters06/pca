#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 19:23:10 2020

@author: hernando
@collaborator: paternina
This code creates B-splines surfaces
"""

import numpy as np
import pca_01 as pca
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

Px = np.array([[0,3,6,9],[0,3,6,9],[0,3,6,9]])
Py = np.array([[0,0,0,0],[2,2,2,2],[4,4,4,4]])
Pz = np.array([[0,3,3,0],[2,5,5,2],[0,3,3,0]])

U = np.array([0,0,0,0.5,1,1,1])
V = np.array([0,0,0,1,1,1])
p = 2
q = 2

# U = np.array([0,0,0,0.4,0.6,1,1,1])
# V = np.array([0,0,0,0.2,0.5,0.8,1,1,1])
# p = 2
# q = 2

# U = [0,0,0,0,0.25,0.5,0.75,1,1,1,1]
# V = [0,0,0,0.2,0.4,0.6,0.6,0.8,1,1,1]
# p = 3
# q = 2

# u = bs.uGenerator(0.0,1.0,0.01)
# v = bs.uGenerator(0.0,1.0,0.01)

cx,cy,cz = pca.bSplineSurface(U,V,p,q,Px,Py,Pz)
pca.plottingSurface(cx,cy,cz,Px,Py,Pz)