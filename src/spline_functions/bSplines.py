import numpy as np

import basisFunctions as bfunc

####################################################
####################B-SPLINES#######################
####################################################

################B-SPLINES CURVES####################

def bSpline1DCurve(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))

    for i in range(len(urank)):
        nVec = bfunc.basisFunction(U,p,urank[i])
        cx[i] = nVec@px.T
    return cx

def bSplineCurve(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cpts = np.zeros((numpoints,2))

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))
    
    mu = len(U) - 1
    nu = mu - p - 1

    for i in range(len(urank)):
        uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
        idxU = uspan + np.arange(0,p+1) - p
        nbas = bfunc.basisFunction(uspan,urank[i],mu,p,U)
        
        nVec = np.zeros(nu+1)
        nVec[idxU] = nbas
        nVec = np.reshape(nVec,(1,len(nVec)))
        
        cpts[i][0] = nVec@px
        cpts[i][1] = nVec@py
    return cpts

def bSplineCurveDerivative(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cppts = np.zeros((numpoints,2))

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))
    
    mu = len(U) - 1
    nu = mu - p - 1
    
    # Derivative order
    d = 1

    for i in range(len(urank)):
        uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
        idxU = uspan + np.arange(0,p+1) - p
        dnbas = bfunc.derBasisFunction(uspan,urank[i],mu,p,U,d)
        
        dnVec = np.zeros(nu+1)
        dnVec[idxU] = dnbas[d,:]
        dnVec = np.reshape(dnVec,(1,len(dnVec)))
        
        cppts[i][0] = dnVec@px
        cppts[i][1] = dnVec@py
    return cppts

################B-SPLINES SURFACES###################

def bivariateShapeFunction(U,V,p,q,u,v):
    nVecU = bfunc.basisFunction(U,p,u)
    nVecV = bfunc.basisFunction(V,q,v)

    n2 = nVecV.T@nVecU
    n2 = np.reshape(n2,(1,n2.shape[0]*n2.shape[1]))
    return n2

def bivariateDerivativeU(U,V,p,q,u,v):
    dnVecU = bfunc.derNFunction(U,p,u)
    nVecV = bfunc.basisFunction(V,q,v)

    dn2u = nVecV.T@dnVecU
    dn2u = np.reshape(dn2u,(1,dn2u.shape[0]*dn2u.shape[1]))
    return dn2u

def bivariateDerivativeV(U,V,p,q,u,v):
    nVecU = bfunc.basisFunction(U,p,u)
    dnVecV = bfunc.derNFunction(V,q,v)

    dn2v = dnVecV.T@nVecU
    dn2v = np.reshape(dn2v,(1,dn2v.shape[0]*dn2v.shape[1]))
    return dn2v

def bSpline2DField(U,V,p,q,px,py):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            n2 = bivariateShapeFunction(U,V,p,q,urank[i],vrank[j])

            cx[i,j] = n2@px
            cy[i,j] = n2@py
    return cx,cy

def bSplineSurface2(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            n2 = bivariateShapeFunction(U,V,p,q,urank[i],vrank[j])

            cx[i,j] = n2@px
            cy[i,j] = n2@py
            cz[i,j] = n2@pz
    return cx,cy,cz


