import numpy as np
import pca_01 as pca
import nurbs as rbs
import plottingScripts as plts
import matplotlib.pyplot as plt
import refinements as rfn

Pinp = np.array([[0.5,0],[0.5,0.5],[0,0.5],[1,0],[1,1],[0,1]])
winp = np.array([[1],[0.5*np.sqrt(2)],[1],[1],[0.5*np.sqrt(2)],[1]])
pinp = 2
qinp = 1
Uinp = np.array([0,0,0,1,1,1])
Vinp = np.array([0,0,1,1])

cx,cy = rbs.nurbs2DField(Uinp,Vinp,pinp,qinp,Pinp,winp)
plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])))
