#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 18:16:56 2020
Modified on Sun Aug 02 17:59:20 2020
Modified on Sat Jul 01 16:56:20 2023

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
from src.multipatch_nurbs import MultiPatchNURBSCurve


def plot_curve1():
    P = np.array([[0,0],[0.5,1],[1,0]])
    w = np.array([[1],[1],[1]])
    U = np.array([0,0,0,1,1,1])
    p = 2

    curve1 = rbs.NURBSCurve(P,w,p,U)
    # cpts = curve1.createCurve()
    curve1.plotBasisFunctions()
    # curve1.plotCurve()


def plot_curve2():
    # P = np.array([[0,0],[1,1],[3,2],[4,1],[5,-1]])
    P = np.array([[0,0],[0.5,1],[1,0],[1,0],[1.5,-1],[2,0]])
    U = np.array([0,0,0,1,1,2,2,2])
    w = np.array([[1],[1],[1],[1],[1]])
    p = 2

    curve2 = rbs.NURBSCurve(P,w,p,U)
    cpts = curve2.createCurve()
    curve2.plotCurve()


def plot_basis_functions():
    Rmax = 1.0
    P = np.array([[Rmax,0],[Rmax,Rmax],[0,Rmax],[-Rmax,Rmax],[-Rmax,0]])
    w = np.array([[1],[0.5*np.sqrt(2)],[1],[0.5*np.sqrt(2)],[1]])
    U = np.array([0,0,0,1,2,3,4,4,4])
    p = 2

    curve1 = rbs.NURBSCurve(P,w,p,U)
    curve1.plotBasisFunctions()


def plot_patch_curves():
    P1 = np.array([[0,0],[0.5,1],[1,0]])
    w1 = np.array([[1],[1],[1]])
    U1 = np.array([0,0,0,0.5,1,1,1])
    p1 = 2
    
    P2 = np.array([[1,0],[1.5,-1],[2,0]])
    w2 = np.array([[1],[1],[1]])
    U2 = np.array([0,0,0,1,1,1])
    p2 = 2

    curve1 = rbs.NURBSCurve(P1,w1,p1,U1)
    curve2 = rbs.NURBSCurve(P2,w2,p2,U2)

    curve_patch = MultiPatchNURBSCurve(curve1, curve2)
    curve_patch.plot_multipatch_basis_functions()
    curve_patch.plot_multipatch_curve()


def plot_patch_curves_2():
    P1 = np.array([[0,0],[0.5,1],[1,0]])
    w1 = np.array([[1],[0.25],[0.25],[1],[1]])
    U1 = np.array([0,0,0,0,0.5,1,1,1,1])
    p1 = 3
    
    P2 = np.array([[1,0],[1.5,-1],[2,0]])
    w2 = np.array([[0.5],[1],[1],[0.75],[1]])
    U2 = np.array([0,0,0,0,0.25,0.75,1,1,1,1])
    p2 = 3

    curve1 = rbs.NURBSCurve(P1,w1,p1,U1)
    curve2 = rbs.NURBSCurve(P2,w2,p2,U2)

    curve_patch = MultiPatchNURBSCurve(curve1, curve2)
    curve_patch.plot_multipatch_basis_functions()
    # curve_patch.plot_multipatch_curve()


if __name__ == '__main__':
    id_curve = 4
    if id_curve == 1:
        plot_curve1()
    elif id_curve == 2:
        plot_curve2()
    elif id_curve == 3:
        plot_basis_functions()
    elif id_curve == 4:
        plot_patch_curves()
    else:
        plot_patch_curves_2()
