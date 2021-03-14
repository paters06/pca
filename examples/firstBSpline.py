#First B-spline curve
import numpy as np
import matplotlib.pyplot as plt
import basisFunctions as bfunc
import bSplines as bs
import plottingScripts as plts

# Examples with curves

#P = np.array([[-1,0],[-1,1],[1,1],[1,0]])
#U = np.array([0,0,0,1/2,1,1,1])
#p = 2

# P = np.array([[0,1,3,5,5,8,9],[0,2,4,2,0,0,3]])
# U = np.array([0,0,0,0,2,3,3,5,5,5,5])
# p = 3

# P = np.array([[1,2,3,4,5,6,7,8],[1,1,1,1,0,0,0,0]])
# p = 2
# U = pca.knotGeneratorUniform(len(P[0])-1,p)

# num = 5
# P = np.array([np.linspace(0,1,num),np.zeros(num)])
# p = 2
# U = np.array([0,0,0,1/3,2/3,1.0,1.0,1.0])

#U = np.array([0,0,0,1,2,3,4,5,5,5])
#V = np.array([0,0,0.5,1,1])
#p = 2
#q = 1

#u = 2.5
#v = 0.0563508326896

#na = bfunc.basisFunction(u,len(U)-1,p,U)
#dna = bfunc.derBasisFunction(u,len(U)-1,p,U,1)
#nb = pca.nFunction(V,q,v)
#dnb = pca.derNFunction(V,q,v)

#dnadu = nb.T@dna

#print(na)
#print(dna)
#print(nb)
# print(dnb)

#print(dnadu)
#print(np.reshape(dnadu,(1,dnadu.shape[0]*dnadu.shape[1])))

#B-Spline curve
#cpts = bs.bSplineCurve(U,p,P)
#plts.plotCurve2d(cpts,P)

# Examples with derivatives

#Derivative B-Spline

P = np.array([[1,1],[2,5],[6,5],[7,1],[9,3]])
U = np.array([0,0,0,0.4,0.6,1,1,1])
p = 2

cpts = bs.bSplineCurve(U,p,P)
cppts = bs.bSplineCurveDerivative(U,p,P)
#print(cppts)

#plts.plotCurve2d(cpts,P,"no")
plts.plotTangentCurve2d(cpts,cppts_norm,P,"no")

# plts.plotCurve2d(cX,cY,Px,Py,"yes","test1")
#plts.plotTangentCurve2d(cX,cY,cPrimeX,cPrimeY,P,"yes","test2")
