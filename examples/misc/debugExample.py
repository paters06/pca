# Python libraries
import numpy as np

import sys
sys.path.append('/home/paters/pca/')

# Local project
import src.preprocessor2D as pre2D
import src.debug.debugScripts as dbg_scrpt

Rmax = 1.0
Rmin = 0.5
tv = 10

P = np.array([[Rmin,0],[Rmin,Rmin],[0,Rmin],[Rmax,0],[Rmax,Rmax],[0,Rmax]])

w = np.array([[1],[0.5*np.sqrt(2)],[1],[1],[0.5*np.sqrt(2)],[1]])

#Isogeometric routines
U = np.array([0,0,0,1,1,1])
V = np.array([0,0,1,1])

p = 2
q = 1

numGaussPoints = 4
# gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)
gaussLegendreQuadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)
# print(gaussLegendreQuadrature1D)

displacementConditions = [[[0.0,Rmin],[0.0,Rmax],"S",0.0],[[Rmin,0.0],[Rmax,0.0],"S",0.0]]
neumannConditions = [[[0.0,0.0],[1.0,0.0],"normal",tv]]

parametricNodes,nodesInElement = pre2D.parametricGrid(U,V)
loadElements,loadFaces = pre2D.loadPreprocessing(parametricNodes,nodesInElement,neumannConditions)
dirichletBCList = pre2D.dirichletBCPreprocessingOnFaces(P,displacementConditions)

dbg_scrpt.calculateAreaAndLength(U,V,w,p,q,P,parametricNodes,nodesInElement,gaussLegendreQuadrature,loadElements,loadFaces)
