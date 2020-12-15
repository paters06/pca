import numpy as np
from numpy.linalg import inv,det,solve
# from scipy import linalg as la
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import pca_01 as pca
import nurbs as rbs
# import curveFitting as cfit
import plottingScripts as plts
import surfaceRefinements as srfn

def checkingSymmetricMatrix(A):
    check = np.allclose(A, A.T, rtol=1e-2, atol=1e-3)
    if check:
        print("The given matrix is symmetric")
    else:
        print("The given matrix is not symmetric")

def checkRankMatrix(A):
    mRank = np.linalg.matrix_rank(A)
    mRows = A.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        print("The matrix has full rank. It is invertible")
    else:
        print("The matrix hast not full rank. It is not invertible")

def plotSparsity(A):
    fig = plt.figure()
    # plt.spy(A,markersize=5)
    plt.imshow(A,cmap=cm.viridis)
    plt.colorbar()
    plt.show()

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

def jacobian(U,V,w,p,q,pta,ptb,px,py,paramgrad):
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

    jacob = jacob@paramgrad

    return jacob

def weightedJacobian(jac,gwpts,ipta,iptb):
    return abs(det(jac))*gwpts[ipta]*gwpts[iptb]

def strainDisplacementMatrix(U,V,w,p,q,pta,ptb,jacob):
    dN2u = rbs.dRatdU(U,V,w,p,q,pta,ptb)
    dN2v = rbs.dRatdV(U,V,w,p,q,pta,ptb)

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
            # print(coor)
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py,paramgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            bMat = strainDisplacementMatrix(U,V,w,p,q,coor[0][0],coor[0][1],jac)
            lke += (bMat.T@dmat@bMat)*wJac

    # print(lke)
    return lke

def localBodyVector(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt,rho):
    lbe = np.zeros((2*px.shape[0],1))
    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8

    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py,paramgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            nMat = shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
            lbe += (nMat.T@bvec)*wJac

    return lbe

def appliedLoadVector(U,V,w,p,q,px,py,gausspoints,gaussweights,apt,bpt,load):
    lle = np.zeros((2*px.shape[0],1))
    tvec = np.zeros((2,1))
    tvec[0][0] = -load

    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        # coor = parametricCoordinate(uval,uval,vseg[0],vseg[1],gausspoints[qj],gausspoints[qj])
        coor = parametricCoordinate(apt[0],bpt[0],apt[1],bpt[1],gausspoints[qj],gausspoints[qj])
        gcoor = geometricCoordinate(coor,U,V,w,p,q,px,py)
        # print("Geometric coor")
        # print(gcoor)
        # Extracting the unique non-zero value
        jac2 = 0.5*float(np.extract((bpt-apt) > 1e-5, (bpt-apt)))
        # print(jac2)
        du = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])
        # print(du)
        dxdu = du@px
        dydu = du@py
        jac1 = np.sqrt(dxdu**2 + dydu**2)
        # print(jac1)
        nMat = shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
        # lle += (nMat.T@tvec)*4.0*jac2*gaussweights[qj]
        lle += (nMat.T@tvec)*jac1*jac2*gaussweights[qj]

    return lle

def elementArea(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt):
    elemA = 0
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py,paramgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            elemA += 1.0*wJac

    return elemA

def elementLength(U,V,w,p,q,px,py,gausspoints,gaussweights,apt,bpt):
    elemL = 0
    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        coor = parametricCoordinate(apt[0],bpt[0],apt[1],bpt[1],gausspoints[qj],gausspoints[qj])
        gcoor = geometricCoordinate(coor,U,V,w,p,q,px,py)
        # print("Geometric coor")
        # print(gcoor)
        # Extracting the unique non-zero value
        jac2 = 0.5*float(np.extract((bpt-apt) > 1e-5, (bpt-apt)))
        # print(jac2)
        dxdu = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])@px
        dydu = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])@py
        jac1 = np.sqrt(dxdu**2 + dydu**2)
        # print(jac1)
        elemL += 1.0*jac1*jac2*gaussweights[qj]

    return elemL

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(U,V,w,p,q,P,paramnodes,nodeselem,gaussquad,dmat,rho,loadnodes,loadelems,load):
    K = np.zeros((2*P.shape[0],2*P.shape[0]))
    F = np.zeros((2*P.shape[0],1))
    Fb = np.zeros((2*P.shape[0],1))
    Fl = np.zeros((2*P.shape[0],1))
    totalArea = 0
    totalLength = 0
    gaussLegendrePoints = gaussquad[0]
    gaussLegendreWeights = gaussquad[1]

    paramGrad = np.zeros((2,2))
    numElems = nodeselem.shape[0]

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

    for ielem in range(0,numElems):
        """
        - paramGrad is conformed by the difference between u values and v values
        - The difference in u values is determined by paramnodes[C][0] - paramnodes[A][0]
          where C and A stand for the corners of the parametric element ABCD and
          0 stands for the u component
        - C can be obtained as the 3rd element in a given row of the matrix nodeselem:
        ielem -> [A B C D]
                  0 1 2 3
        - Likewise, A is the 1st element in a ielem row of nodeselem
        - It means that uC = paramnodes[nodeselem[ielem][2]][0]
                        uA = paramnodes[nodeselem[ielem][0]][0]
                        vC = paramnodes[nodeselem[ielem][2]][1]
                        vA = paramnodes[nodeselem[ielem][0]][1]
        """
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        paramGrad[0][0] = 0.5*(uC - uA)
        paramGrad[1][1] = 0.5*(vC - vA)

        aPoint = np.array([uA,vA])
        cPoint = np.array([uC,vC])

        print("---")
        print("Element #",ielem)
        K += localStiffnessMatrix(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint,dmat)
        Fb += localBodyVector(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint,rho)
        totalArea += elementArea(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint)
        if ielem in loadelems:
            print('Included element')
            uB = paramnodes[loadnodes[1]][0]
            uA = paramnodes[loadnodes[0]][0]
            vB = paramnodes[loadnodes[1]][1]
            vA = paramnodes[loadnodes[0]][1]

            aPoint = np.array([uA,vA])
            bPoint = np.array([uB,vB])

            Fl += appliedLoadVector(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,aPoint,bPoint,load)
            totalLength += elementLength(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,aPoint,bPoint)

        print("---")

    # print("Total Length")
    # print(totalLength)
    F = Fb + Fl
    return K,F,totalArea

def boundaryConditionsEnforcement(K,F,udisp,axisrestrict,ucond):
    remdofs = 2*(np.array(udisp) - 1)
    restricteddofs = np.array(axisrestrict)
    remdofs = remdofs + restricteddofs

    # remdofs = np.hstack((dofs1,dofs2))
    remdofs.sort()
    # print(remdofs)

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

def nurbs2DField(numpoints,U,V,p,q,P,w,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

    cx = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    cy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):
                ratFunc = rbs.ratFunction(U,V,w,p,q,urank[i],vrank[j])

                cx[iu*numpoints + i,jv*numpoints + j] = ratFunc@px
                cy[iu*numpoints + i,jv*numpoints + j] = ratFunc@py

    return cx,cy

def displacementField(numpoints,U,V,p,q,D,w,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    dx = np.reshape(D[:,0],(D.shape[0],1))
    dy = np.reshape(D[:,1],(D.shape[0],1))

    ux = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    uy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):
                ratFunc = rbs.ratFunction(U,V,w,p,q,urank[i],vrank[j])

                ux[iu*numpoints + i,jv*numpoints + j] = ratFunc@dx
                uy[iu*numpoints + i,jv*numpoints + j] = ratFunc@dy

    return ux,uy

def stressField(numpoints,U,V,p,q,P,w,dtot,dmat,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    paramgrad = np.zeros((2,2))

    sx = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    sy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    sxy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):

                if abs(urank[i] - 0.5) > 1e-5:
                    xpcoor = urank[i]
                    ypcoor = vrank[j]
                else:
                    # print("Singularity")
                    # print(urank[i])
                    # print(vrank[j])
                    xpcoor = 0.51
                    ypcoor = vrank[j]

                paramgrad[0][0] = 0.5*(uC - uA)
                paramgrad[1][1] = 0.5*(vC - vA)

                jac = jacobian(U,V,w,p,q,xpcoor,ypcoor,P[:,0],P[:,1],paramgrad)
                bmat = strainDisplacementMatrix(U,V,w,p,q,xpcoor,ypcoor,jac)

                svec = dmat@(bmat@dtot)
                sx[iu*numpoints + i,jv*numpoints + j] = svec[0]
                sy[iu*numpoints + i,jv*numpoints + j] = svec[1]
                sxy[iu*numpoints + i,jv*numpoints + j] = svec[2]

    return sx,sy,sxy

def postProcessing(U,V,p,q,P,D,w,paramnodes,nodeselem,dtot,dmat):
    numpoints = 11
    cx,cy = nurbs2DField(numpoints,U,V,p,q,P,w,paramnodes,nodeselem)
    ux,uy = displacementField(numpoints,U,V,p,q,D,w,paramnodes,nodeselem)
    sx,sy,sxy = stressField(numpoints,U,V,p,q,P,w,dtot,dmat,paramnodes,nodeselem)
    # plts.plotting2DField(cx,cy,ux,P,["Ux Displacement Field","[m]"])
    plts.plotting2DField(cx,cy,sx,P,["Sx Component Stress Field","[Pa]"])
    # print(sx.max())
    # print(sx.min())

################ PREPROCESSING ####################

def parametricGrid(U,V):
    # Selecting the unique values of each array
    uniqueU = np.unique(U)
    uniqueV = np.unique(V)

    # Creating the 2D mesh
    ugrid,vgrid = np.meshgrid(uniqueU,uniqueV)

    # Resizing the grid component matrix to column vectors
    ucomp = np.reshape(ugrid,(ugrid.shape[0]*ugrid.shape[1],1))
    vcomp = np.reshape(vgrid,(vgrid.shape[0]*vgrid.shape[1],1))

    # Stacking the components of the parametric coordinates
    paramnodes = np.hstack((ucomp,vcomp))

    # Assembling the element matrix
    numU = len(uniqueU)
    numV = len(uniqueV)
    numelems = (numU - 1)*(numV - 1)

    elemmat = np.zeros((numelems,4),dtype=int)
    elemIndex = 0

    for j in range(0,numV-1):
        for i in range(0,numU-1):
            elemmat[elemIndex,0] = j*numU + i #A
            elemmat[elemIndex,1] = j*numU + (i+1) #B
            elemmat[elemIndex,2] = (j+1)*numU + (i+1) #C
            elemmat[elemIndex,3] = (j+1)*numU + i #D

            elemIndex += 1

    return paramnodes,elemmat

def loadPreprocessing(paramnodes,nodeselem,U,V,p,q,P,w,cload):
    loadnodes = []
    loadelements = []
    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

    for i in range(0,paramnodes.shape[0]):
        ratFunc = rbs.ratFunction(U,V,w,p,q,paramnodes[i][0],paramnodes[i][1])
        cx = ratFunc@px
        if abs(cx - cload) < 1e-4:
            loadnodes.append(i)

    for i in range(0,nodeselem.shape[0]):
        present = 0
        for ln in loadnodes:
            x = np.where(nodeselem[i,:] == ln)

            # x is a tuple. where does not have as output a list
            if len(x[0]) != 0:
                present += 1

        if present == 2:
            loadelements.append(i)

    return loadnodes,loadelements

def dirichletBCPreprocessing(P,cdirichlet):
    dirichletctrlpts = []
    axisrestrictions = []

    for i in range(0,P.shape[0]):
        if abs(P[i][0] - cdirichlet) < 1e-4:
            dirichletctrlpts.append(i+1)
            axisrestrictions.append(0)

        if abs(P[i][1] - cdirichlet) < 1e-4:
            dirichletctrlpts.append(i+1)
            axisrestrictions.append(1)

    return dirichletctrlpts,axisrestrictions

####################################################
################## MAIN PROBLEM ####################
####################################################

#Data
E = 1e5 #Pa
nu = 0.31
rho = 0.0 #kg/m3
u0 = 0.0
tv = 10 #Pa
# uDirichlet = [1,4,5,8,9,12]
# uAxis = [1,0,1,0,1,0]
# uNeumann = 0.5
# vNeumann = 1

numGaussPoints = 4
gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

Pinit = np.array([[-1,0],[-1,np.sqrt(2)-1],[1-np.sqrt(2),1],[0,1],[-2.5,0],[-2.5,0.75],
              [-0.75,2.5],[0,2.5],[-4,0],[-4,4],[-4,4],[0,4]])
winit = np.array([[1],[0.5*(1+(1/np.sqrt(2)))],[0.5*(1+(1/np.sqrt(2)))],[1],[1],[1],
              [1],[1],[1],[1],[1],[1]])

#Isogeometric routines
Uinit = np.array([0,0,0,0.5,1,1,1])
Vinit = np.array([0,0,0,1,1,1])

pinit = 2
qinit = 2

doRefinement = 'N'

if doRefinement == 'Y':
    reflist = ['k','k']
    dirlist = ['U','V']
    Uinp,Vinp,pinp,qinp,Pinp,winp = srfn.surfaceRefinement(reflist,dirlist,Uinit,Vinit,pinit,qinit,Pinit,winit)
else:
    Uinp = Uinit
    Vinp = Vinit
    pinp = pinit
    qinp = qinit
    Pinp = Pinit
    winp = winit

# cx,cy = rbs.nurbs2DField(Uinit,Vinit,p,q,P,w)
# plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])))

parametricNodes,nodesInElement = parametricGrid(Uinp,Vinp)
loadNodes,loadElements = loadPreprocessing(parametricNodes,nodesInElement,Uinp,Vinp,pinp,qinp,Pinp,winp,-4.0)
dirichletCtrlPts,axisRestrictions = dirichletBCPreprocessing(Pinp,0.0)
# print(dirichletCtrlPts)
# print(axisRestrictions)

dMat = elasticMatrix(E,nu)
K,F,totalArea = assemblyWeakForm(Uinp,Vinp,winp,pinp,qinp,Pinp,parametricNodes,nodesInElement,gaussLegendreQuadrature,dMat,rho,loadNodes,loadElements,tv)

Kred,Fred,removedDofs = boundaryConditionsEnforcement(K,F,dirichletCtrlPts,axisRestrictions,u0)

totaldofs = np.arange(2*Pinp.shape[0])
dtotal,dred = solveMatrixEquations(Kred,Fred,totaldofs,removedDofs)

dx = dtotal[0::2]
dy = dtotal[1::2]
D = np.hstack((dx,dy))

postProcessing(Uinp,Vinp,pinp,qinp,Pinp,D,winp,parametricNodes,nodesInElement,dtotal,dMat)
