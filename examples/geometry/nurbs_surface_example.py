#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:26:33 2020

@author: hernando
@collaborator: paternina
"""

# Python libraries
import numpy as np

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# Local project
from src.spline_functions.nurbs_surface import NURBSSurface
from src.plotting.nurbs_plotter import NURBSPlotter

P = np.array([[0,0,0],
              [0.2,0,0],
              [0,0.6,0],
              [0.2,0.4,0]])

w = np.array([[1],[1],[1],[1]])

p = 1
q = 1

U = np.array([0,0,1,1])
V = np.array([0,0,1,1])

surface1 = NURBSSurface(P,w,p,q,U,V)
plotter = NURBSPlotter()
plotter.nurbs_surface = surface1
# plotter.plot_surface()
plotter.plot_tangent_surface("u")
