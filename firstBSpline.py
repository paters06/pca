#First B-spline curve
import numpy as np
import matplotlib.pyplot as plt

import pca_01 as pca

# Px = np.array([-1,-1,1,1])
# Py = np.array([0,1,1,0])
# U = np.array([0,0,0,1/2,1,1,1])
# p = 2

Px = np.array([0,1,3,5,5,8,9])
Py = np.array([0,2,4,2,0,0,3])
U = np.array([0,0,0,0,2,3,3,5,5,5,5])
p = 3

cx,cy = pca.bSplineCurve(U,p,Px,Py)
pca.plotCurve2d(cx,cy,Px,Py)