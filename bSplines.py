import numpy as np
import pca_01 as pca

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
        nVec = pca.nFunction(U,p,urank[i])
        cx[i] = nVec@px.T
    return cx

def bSplineCurve(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        nVec = pca.nFunction(U,p,urank[i])
        cx[i] = nVec@px.T
        cy[i] = nVec@py.T
    return cx,cy

def bSplineCurveDerivative(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cprimex = np.zeros(len(urank))
    cprimey = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        dnVec = pca.derNFunction(U,p,urank[i])
        cprimex[i] = dnVec@px.T
        cprimey[i] = dnVec@py.T
    return cprimex,cprimey

################B-SPLINES SURFACES###################

def bSplineSurface(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = pca.nFunction(U,p,urank[i])
            nVecj = pca.nFunction(V,q,vrank[j])

            cx[i,j] = nVecj.T@px@nVeci
            cy[i,j] = nVecj.T@py@nVeci
            cz[i,j] = nVecj.T@pz@nVeci
    return cx,cy,cz

def bSplineSurfaceDerivativeU(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            dnVeci = pca.derNFunction(U,p,urank[i])
            nVecj = pca.nFunction(V,q,vrank[j])

            cprimex[i,j] = nVecj.T@px@dnVeci
            cprimey[i,j] = nVecj.T@py@dnVeci
            cprimez[i,j] = nVecj.T@pz@dnVeci
    return cprimex,cprimey,cprimez

def bSplineSurfaceDerivativeV(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = pca.nFunction(U,p,urank[i])
            dnVecj = pca.derNFunction(V,q,vrank[j])

            cprimex[i,j] = dnVecj.T@px@nVeci
            cprimey[i,j] = dnVecj.T@py@nVeci
            cprimez[i,j] = dnVecj.T@pz@nVeci
    return cprimex,cprimey,cprimez

def bSplineSurfaceDerivativeUV(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            dnVeci = pca.derNFunction(U,p,urank[j])
            dnVecj = pca.derNFunction(V,q,vrank[j])

            cprimex[i,j] = dnVecj.T@px@dnVeci
            cprimey[i,j] = dnVecj.T@py@dnVeci
            cprimez[i,j] = dnVecj.T@pz@dnVeci
    return cprimex,cprimey,cprimez
