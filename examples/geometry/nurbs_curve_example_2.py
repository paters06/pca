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

import os
import sys

# print(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

# Local project
from src.spline_functions.nurbs_curve import NURBSCurve
from src.plotting.nurbs_plotter import NURBSPlotter
from src.spline_functions.curveRefinements import curveRefinement

Rmax = 1.0

P = np.array([[Rmax,0],[Rmax,Rmax],[0,Rmax]])
w = np.array([[1],[0.5*np.sqrt(2)],[1]])
U = np.array([0,0,0,1,1,1])
p = 2


curve1 = NURBSCurve(P,w,p,U)
plotter = NURBSPlotter()
plotter.nurbs_curve = curve1
# plotter.plot_curve()
# plotter.plot_tangent_curve()

reflist = ['p']
Uref,pref,Pref,wref = curveRefinement(reflist,U,p,P,w)

curve2 = NURBSCurve(Pref,wref,pref,Uref)
plotter.nurbs_curve = curve2
print(Pref)
# plotter.plot_curve()