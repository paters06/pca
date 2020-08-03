#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:26:33 2020

@author: hernando
@collaborator: paternina
"""
import pca_01 as pca
import numpy as np

Px = np.array([[0,3,6,9],[0,3,6,9],[0,3,6,9]])
Py = np.array([[0,0,0,0],[2,2,2,2],[4,4,4,4]])
Pz = np.array([[0,3,3,0],[2,5,5,2],[0,3,3,0]])
w = np.array([[1,3,3,1],[1,5,5,1],[1,4,4,1]])
p = 2
q = 2

U = np.array([0,0,0,0.5,1,1,1])
V = np.array([0,0,0,1,1,1])

cx,cy,cz = pca.nurbsSurface(U,V,p,q,Px,Py,Pz,w)
pca.plottingSurface(cx,cy,cz,Px,Py,Pz)
