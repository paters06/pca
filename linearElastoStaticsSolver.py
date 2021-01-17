import numpy as np
import numpy.linalg
import nurbs as rbs

def parametricCoordinate(ua,ub,va,vb,gausspta,gaussptb):
    localpts = np.zeros((1,2))
    localpts[0][0] = 0.5*(ub - ua)*gausspta + 0.5*(ub + ua)
    localpts[0][1] = 0.5*(vb - va)*gaussptb + 0.5*(vb + va)
    return localpts

def geometricCoordinate(paramcoor,U,V,w,p,q,px,py):
    ratFunc = rbs.ratFunction(U,V,w,p,q,paramcoor[0][0],paramcoor[0][1])
    geomcoor = np.zeros((1,2))
    geomcoor[0][0] = ratFunc@px
    geomcoor[0][1] = ratFunc@py
    return geomcoor

def jacobian(U,V,w,p,q,pta,ptb,px,py):
    n2 = rbs.ratFunction(U,V,w,p,q,pta,ptb)
    dn2u = rbs.dRatdU(U,V,w,p,q,pta,ptb)
    dn2v = rbs.dRatdV(U,V,w,p,q,pta,ptb)

    dXdu = dn2u@px
    dXdv = dn2v@px
    dYdu = dn2u@py
    dYdv = dn2v@py

    jacob = np.zeros((2,2))

    jacob[0][0] = dXdu
    jacob[0][1] = dXdv
    jacob[1][0] = dYdu
    jacob[1][1] = dYdv

    return jacob

def weightedJacobian(jac,paramgrad,gwpts,ipta,iptb):
    return abs(np.linalg.det(jac))*abs(np.linalg.det(paramgrad))*gwpts[ipta]*gwpts[iptb]

def strainDisplacementMatrix(U,V,w,p,q,pta,ptb,jacob):
    dN2u = rbs.dRatdU(U,V,w,p,q,pta,ptb)
    dN2v = rbs.dRatdV(U,V,w,p,q,pta,ptb)

    invJac = np.linalg.inv(jacob)
    dN2 = np.vstack((dN2u,dN2v))
    dN2dxi = invJac.T@dN2

    numpts = dN2dxi.shape[1]
    bMat = np.zeros((3,2*numpts))
    #dNx
    bMat[0,0::2] = dN2dxi[0,:]
    bMat[2,0::2] = dN2dxi[1,:]
    #dNy
    bMat[1,1::2] = dN2dxi[1,:]
    bMat[2,1::2] = dN2dxi[0,:]
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

def shapeFunctionMatrix(U,V,w,p,q,pta,ptb):
    N2 = rbs.ratFunction(U,V,w,p,q,pta,ptb)
    nMat = np.zeros((2,2*N2.shape[1]))

    nMat[0,0::2] = N2
    nMat[1,1::2] = N2
    return nMat

################ WEAK FORM INTEGRALS ####################

def localStiffnessMatrix(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt,dmat):
    lke = np.zeros((2*px.shape[0],2*px.shape[0]))
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
            wJac = weightedJacobian(jac,paramgrad,gaussweights,qi,qj)
            bMat = strainDisplacementMatrix(U,V,w,p,q,coor[0][0],coor[0][1],jac)
            lke += (bMat.T@dmat@bMat)*wJac

    return lke

def localBodyVector(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt,rho):
    lbe = np.zeros((2*px.shape[0],1))
    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8

    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
            wJac = weightedJacobian(jac,paramgrad,gaussweights,qi,qj)
            nMat = shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
            lbe += (nMat.T@bvec)*wJac

    return lbe
