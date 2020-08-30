import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
import pca_01 as pca
import curveFitting as cfit

def scalingSegment(ua,ub):
    localpts = np.zeros(len(gaussLegendrePoints))
    for i in range(len(localpts)):
        localpts[i] = 0.5*(ub - ua)*gaussLegendrePoints[i] + 0.5*(ub + ua)
    return localpts

def localStiffnessMatrix(intgpt):
    derNB = pca.derNFunction(U,p,intgpt)
    lke = derNB.transpose()*E@derNB
    return lke

def localForceVector(intgpt):
    nB = pca.nFunction(U,p,intgpt)
    if intgpt < 5.0:#Be careful with this condition
        bx = 10.0
    else:
        bx = 0.0
    # print(nB.T*bx)
    return nB.transpose()*bx

####################################################
###################MAIN PROBLEM#####################
####################################################

L = 10.0
numNodes = 21
p = 1

U = pca.knotGeneratorUniform(numNodes-1,p)*L
Ured = U[p:-p] #this guarantees that u for evaluation goes from [0,1]
# print(Ured) #Uncomment this line to check that knot vector is fine
numPatches = len(Ured) - 1

E = 1000
u0 = 0.0
tx = 25.0

numGaussPoints = 3
gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)
gaussLegendrePoints = gaussLegendreQuadrature[0]
gaussLegendreWeights = gaussLegendreQuadrature[1]

K = np.zeros((numNodes,numNodes))
F = np.zeros((numNodes,1))

for i in range(0,numPatches):
    intgPoints = scalingSegment(Ured[i],Ured[i+1])
    lenSegment = Ured[i+1] - Ured[i]
    Ke = np.zeros((numNodes,numNodes))
    Fe = np.zeros((numNodes,1))
    for j in range(len(intgPoints)):
        Ke += 0.5*lenSegment*localStiffnessMatrix(intgPoints[j])*gaussLegendreWeights[j]
        Fe += 0.5*lenSegment*localForceVector(intgPoints[j])*gaussLegendreWeights[j]

    K += Ke
    F += Fe

#Applying Neumann BC
F[-1] += tx

# print(K)
# print(F)

#Applying Dirichlet BC
dirNB = pca.nFunction(U,p,0.0)
#Replacing basis funcion evaluted in u=0
K[0,:] = dirNB

choice = 2

if choice == 1:
    #Removing the row and column with the Dirichlet BC
    Kred = K[1:,1:]
    #Modifying F column for suitable operation
    Kcol = np.reshape(K[:,0],(len(K[:,0]),1))
    #Updating F vector
    F -= Kcol*u0
    Fred = F[1:,:]

    #Solving the linear system
    d = la.solve(Kred,Fred)

    #Assembly the total solution
    dtotal = np.zeros((numNodes,1))
    dtotal[0] = u0
    dtotal[1:] = d
else:
    F[0] = u0
    d = la.solve(K,F)
    dtotal = d

#Construction of geometry
x = np.linspace(0,L,numNodes)
qPoints = np.zeros((numNodes,1))
qPoints[:,0] = x

pgeom = 1
method = "chord"
ctrlptsgeom,Ugeom = cfit.interpolateCurve(qPoints,pgeom,method)
cx = pca.bSpline1DCurve(Ugeom,pgeom,ctrlptsgeom.T)

#Building the solution field
numpoints = 101
urank = np.linspace(U.min(),U.max(),numpoints)

ux = np.zeros(len(urank))
sx = np.zeros(len(urank))

for i in range(len(urank)):
    nVec = pca.nFunction(U,p,urank[i])
    dnVec = pca.derNFunction(U,p,urank[i])
    ux[i] = dtotal.T@nVec.transpose()
    sx[i] = E*dtotal.T@dnVec.T

# plt.figure()
# plt.plot(cx,ux)
# plt.show()

# plt.figure()
# plt.plot(cx,sx)
# plt.show()
