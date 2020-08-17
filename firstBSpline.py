#First B-spline curve
import numpy as np
import matplotlib.pyplot as plt

import pca_01 as pca

# Examples with curves

# Px = np.array([-1,-1,1,1])
# Py = np.array([0,1,1,0])
# U = np.array([0,0,0,1/2,1,1,1])
# p = 2

# Px = np.array([0,1,3,5,5,8,9])
# Py = np.array([0,2,4,2,0,0,3])
# U = np.array([0,0,0,0,2,3,3,5,5,5,5])
# p = 3
#
# #B-Spline curve
# cx,cy = pca.bSplineCurve(U,p,Px,Py)
# pca.plotCurve2d(cx,cy,Px,Py)

# Examples with derivatives

#Derivative B-Spline
# U = np.array([0,0,0,1,2,3,4,4,5,5,5])
# p = 2
# u = 5.0
#
# dnb = pca.derNFunction(U,p,u)
# print(dnb)

Px = np.array([1,2,6,7,9])
Py = np.array([1,5,5,1,3])
U = np.array([0,0,0,0.4,0.6,1,1,1])
p = 2

cX,cY = pca.bSplineCurve(U,p,Px,Py)
# pca.plotCurve2d(cX,cY,Px,Py,"yes","test1")
# pca.plotCurve2d(cX,cY,Px,Py,"no")

cPrimeX,cPrimeY = pca.bSplineCurveDerivative(U,p,Px,Py)
# pca.plotTangentCurve2d(cX,cY,cPrimeX,cPrimeY,Px,Py)
