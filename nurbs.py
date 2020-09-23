import numpy as np
import pca_01 as pca

####################################################
######################NURBS#########################
####################################################

##################NURBS CURVES######################

def nurbsCurve(U,p,P,w):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        nVec = pca.nFunction(U,p,urank[i])
        cx[i] = (px@(nVec*w).transpose())/(nVec@w.transpose())
        cy[i] = (py@(nVec*w).transpose())/(nVec@w.transpose())
    return cx,cy

#################NURBS SURFACES#####################

def nurbsSurface(U,V,p,q,px,py,pz,w):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = pca.nFunction(U,p,urank[i])
            nVecj = pca.nFunction(V,q,vrank[j])

            cx[i,j] = (nVecj@(px*w)@nVeci.transpose())/(nVecj@(w)@nVeci.transpose())
            cy[i,j] = (nVecj@(py*w)@nVeci.transpose())/(nVecj@(w)@nVeci.transpose())
            cz[i,j] = (nVecj@(pz*w)@nVeci.transpose())/(nVecj@(w)@nVeci.transpose())
    return cx,cy,cz
