import numpy as np
from numpy.linalg import inv,det,solve
# from scipy import linalg as la
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import pca_01 as pca
import bSplines as bs
# import curveFitting as cfit
import plottingScripts as plts
# import refinements as rfn

def checkingSymmetricMatrix(A):
    check = np.allclose(A, A.T, rtol=1e-2, atol=1e-3)
    if check:
        print("The given matrix is symmetric")
    else:
        print("The given matrix is not symmetric")

def checkRankMatrix(A):
    mRank = np.linalg.matrix_rank(Kred)
    mRows = Kred.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        print("The matrix has full rank. It is invertible")
    else:
        print("The matrix hast not full rank. It is not invertible")

def plotSparsity(A):
    fig = plt.figure()
    plt.spy(A,markersize=5)
    # plt.imshow(A,cmap=cm.viridis)
    # plt.colorbar()
    plt.show()

def coordinateParametrization(ua,ub,va,vb,gausspta,gaussptb):
    localpts = np.zeros((1,2))
    localpts[0][0] = 0.5*(ub - ua)*gausspta + 0.5*(ub + ua)
    localpts[0][1] = 0.5*(vb - va)*gaussptb + 0.5*(vb + va)
    return localpts

def jacobian(U,V,p,q,pta,ptb,px,py,parentgrad):
    # print(pta)
    # print(ptb)
    n2 = bs.bivariateShapeFunction(U,V,p,q,pta,ptb)
    dn2u = bs.bivariateDerivativeU(U,V,p,q,pta,ptb)
    dn2v = bs.bivariateDerivativeV(U,V,p,q,pta,ptb)

    dXdu = dn2u@px
    dXdv = dn2v@px
    dYdu = dn2u@py
    dYdv = dn2v@py

    jacob = np.zeros((2,2))

    jacob[0][0] = dXdu
    jacob[0][1] = dXdv
    jacob[1][0] = dYdu
    jacob[1][1] = dYdv

    # print("N2")
    # print(n2@px)
    # print("dXdu")
    # print(dn2u)
    jacob = jacob@parentgrad

    return jacob

def weightedJacobian(jac,gwpts,ipta,iptb):
    return abs(det(jac))*gwpts[ipta]*gwpts[iptb]

def strainDisplacementMatrix(U,V,p,q,pta,ptb,jacob):
    dN2u = bs.bivariateDerivativeU(U,V,p,q,pta,ptb)
    dN2v = bs.bivariateDerivativeV(U,V,p,q,pta,ptb)

    # print(jacob)
    invJac = inv(jacob)
    # print(invJac)
    dN2 = np.vstack((dN2u,dN2v))
    dN2dxi = invJac@dN2

    numpts = dN2dxi.shape[1]
    bMat = np.zeros((3,2*numpts))
    #dNx
    bMat[0,0::2] = dN2dxi[0]
    bMat[2,0::2] = dN2dxi[0]
    #dNy
    bMat[1,1::2] = dN2dxi[1]
    bMat[2,1::2] = dN2dxi[1]
    return bMat

def elasticMatrix(E,nu):
    dmat = np.zeros((3,3))
    dmat[0][0] = 1 - nu
    dmat[1][1] = 1 - nu
    dmat[2][2] = (1 - 2*nu)/2
    dmat[0][1] = nu
    dmat[1][0] = nu
    dmat *= E/((1+nu)*(1-2*nu))
    return dmat

def shapeFunctionMatrix(U,V,p,q,pta,ptb):
    N2 = bs.bivariateShapeFunction(U,V,p,q,pta,ptb)
    nMat = np.zeros((2,2*N2.shape[1]))

    nMat[0,0::2] = N2
    nMat[1,1::2] = N2
    return nMat

################ WEAK FORM INTEGRALS ####################

def localStiffnessMatrix(U,V,useg,vseg,p,q,px,py,gausspoints,gaussweights,parentgrad,dmat):
    lke = np.zeros((2*px.shape[0],2*px.shape[0]))
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = coordinateParametrization(useg[0],useg[1],vseg[0],vseg[1],gausspoints[qi],gausspoints[qj])
            # print(coor)
            jac = jacobian(U,V,p,q,coor[0][0],coor[0][1],px,py,parentgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            bMat = strainDisplacementMatrix(U,V,p,q,coor[0][0],coor[0][1],jac)
            lke += (bMat.T@dmat@bMat)*wJac

    # print(lke)
    return lke

def localBodyVector(U,V,useg,vseg,p,q,px,py,gausspoints,gaussweights,parentgrad,rho):
    lbe = np.zeros((2*px.shape[0],1))
    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8

    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = coordinateParametrization(useg[0],useg[1],vseg[0],vseg[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,p,q,coor[0][0],coor[0][1],px,py,parentgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            nMat = shapeFunctionMatrix(U,V,p,q,coor[0][0],coor[0][1])
            lbe += (nMat.T@bvec)*wJac

    return lbe

def appliedLoadVector(U,V,uval,vseg,p,q,px,py,gausspoints,gaussweights,load):
    lle = np.zeros((2*px.shape[0],1))
    tvec = np.zeros((2,1))
    tvec[0][0] = -load

    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        coor = coordinateParametrization(uval,uval,vseg[0],vseg[1],gausspoints[qj],gausspoints[qj])
        jac = 0.5*(vseg[1] - vseg[0])
        nMat = shapeFunctionMatrix(U,V,p,q,coor[0][0],coor[0][1])
        lle += (nMat.T@tvec)*jac*gaussweights[qj]

    return lle

def elementArea(U,V,useg,vseg,p,q,px,py,gausspoints,gaussweights,parentgrad):
    elemA = 0
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = coordinateParametrization(useg[0],useg[1],vseg[0],vseg[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,p,q,coor[0][0],coor[0][1],px,py,parentgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            elemA += 1.0*wJac

    return elemA

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(U,V,Ured,Vred,p,q,px,py,gaussquad,dmat,rho,uneu,load):
    K = np.zeros((2*px.shape[0],2*px.shape[0]))
    F = np.zeros((2*px.shape[0],1))
    Fb = np.zeros((2*px.shape[0],1))
    Fl = np.zeros((2*px.shape[0],1))
    totalArea = 0
    gaussLegendrePoints = gaussquad[0]
    gaussLegendreWeights = gaussquad[1]

    parentElemGrad = np.zeros((2,2))
    numElems = 0

    for ji in range(len(Vred)-1):
        for ii in range(len(Ured)-1):
            #Inside each element
            numElems += 1

            parentElemGrad[0][0] = 0.5*(Ured[ii+1] - Ured[ii])
            parentElemGrad[1][1] = 0.5*(Vred[ji+1] - Vred[ji])

            print("Element #",numElems)
            print("U coor:",np.array([[Ured[ii],Ured[ii+1]]]))
            print("V coor:",np.array([[Vred[ji],Vred[ji+1]]]))
            print("---")
            usegment = np.array([Ured[ii],Ured[ii+1]])
            vsegment = np.array([Vred[ji],Vred[ji+1]])
            K += localStiffnessMatrix(U,V,usegment,vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,parentElemGrad,dmat)
            Fb += localBodyVector(U,V,usegment,vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,parentElemGrad,rho)
            totalArea += elementArea(U,V,usegment,vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,parentElemGrad)
            if abs(Ured[ii+1] - uneu) < 1e-5:
                Fl += appliedLoadVector(U,V,Ured[ii+1],vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,load)

            print("---")

    F = Fb + Fl
    return K,F,totalArea

def boundaryConditionsEnforcement(K,F,udisp,ucond):
    dofs1 = 2*(np.array(udisp) - 1)
    dofs2 = dofs1 + 1

    remdofs = np.hstack((dofs1,dofs2))
    remdofs.sort()
    print(remdofs)

    print("First reduction")
    Fred = np.delete(F,remdofs,0)
    Kred = np.delete(K,remdofs,0)

    print("Modification of Freduced")
    for u in udisp:
        Kcol = Kred[:,u]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        Fred -= Kcol*ucond

    print("Second reduction")
    Kred = np.delete(Kred,remdofs,1)

    return Kred,Fred,remdofs

def solveMatrixEquations(Kred,Fred,totaldofs,remdofs):
    dred = solve(Kred,Fred)
    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    dtotal = np.zeros((totaldofs.shape[0],1))
    dtotal[reduceddofs,:] = dred
    return dtotal,dred

################ POSTPROCESSING ####################

def postProcessing(U,V,p,q,px,py,dx,dy):
    cx,cy = bs.bSpline2DField(U,V,p,q,px,py)
    ux,uy = bs.bSpline2DField(U,V,p,q,dx,dy)
    # plts.plotting2DField(cx,cy,ux)
    # plts.plotting2DField(cx,cy,uy)

####################################################
###################MAIN PROBLEM#####################
####################################################

#Data
E = 2e11 #Pa
nu = 0.3
rho = 0 #kg/m3
u0 = 0.0
tv = 1e6 #Pa
uDirichlet = [1,2]
uNeumann = 1

numGaussPoints = 4
gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

P = np.array([[-1,0],[-2,0],[-1,0.5],[-2,1],[-0.5,1],[-1,2],[0,1],[0,2]])
# P = np.array([[0,0],[-1,0],[-2,0],[0,1],[-1,2],[-2,2],[1,1.5],[1,4],[1,5],[3,1.5],[3,4],[3,5]])

#Isogeometric routines
Uinit = np.array([0,0,1,1])
Vinit = np.array([0,0,0,0.5,1,1,1])

# Vinit = np.array([0,0,0,0.5,1,1,1])
# Uinit = np.array([0,0,0,1,1,1])

p = 1
# p = 2
q = 2

Ured = Uinit[p:-p]
Vred = Vinit[q:-q]

dMat = elasticMatrix(E,nu)
K,F,totalArea = assemblyWeakForm(Uinit,Vinit,Ured,Vred,p,q,P[:,0],P[:,1],gaussLegendreQuadrature,dMat,rho,uNeumann,tv)

Kred,Fred,removedDofs = boundaryConditionsEnforcement(K,F,uDirichlet,u0)

totaldofs = np.arange(2*P.shape[0])
dtotal,dred = solveMatrixEquations(Kred,Fred,totaldofs,removedDofs)

dx = dtotal[0::2]
dy = dtotal[1::2]

postProcessing(Uinit,Vinit,p,q,P[:,0],P[:,1],dx,dy)
