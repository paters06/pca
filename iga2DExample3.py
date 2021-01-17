import numpy as np
from numpy.linalg import inv,det,solve
import numpy.linalg
import matplotlib.pyplot as plt

import pca_01 as pca
import nurbs as rbs
import plottingScripts as plts
import debugScripts as dbg_scrpt
import preprocessor2D as pre2D
import linearElastoStaticsSolver as linElastStat
import postprocessor2D as post2D
import surfaceRefinements as srfn

################ WEAK FORM INTEGRALS ####################

def appliedLoadVector(U,V,w,p,q,px,py,gausspoints,gaussweights,apt,bpt,load):
    lle = np.zeros((2*px.shape[0],1))
    tvec = np.zeros((2,1))
    tvec[0][0] = -load

    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        coor = linElastStat.parametricCoordinate(apt[0],bpt[0],apt[1],bpt[1],gausspoints[qj],gausspoints[qj])
        gcoor = linElastStat.geometricCoordinate(coor,U,V,w,p,q,px,py)
        # print("Geometric coor")
        # print(gcoor)
        # Extracting the unique non-zero value
        jac2 = 0.5*float(np.extract((bpt-apt) > 1e-5, (bpt-apt)))
        du = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])
        dxdu = du@px
        dydu = du@py
        # print('----')
        # print(dxdu)
        # print(dydu)
        # print('----')
        jac1 = np.sqrt(dxdu**2 + dydu**2)
        nMat = linElastStat.shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
        # print('----')
        # print(nMat.T@tvec)
        # print('----')
        lle += (nMat.T@tvec)*jac1*jac2*gaussweights[qj]

    return lle

def elementArea(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt):
    elemA = 0
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = linElastStat.parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = linElastStat.jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
            wJac = linElastStat.weightedJacobian(jac,paramgrad,gaussweights,qi,qj)
            elemA += 1.0*wJac

    return elemA

def elementLength(U,V,w,p,q,px,py,gausspoints,gaussweights,apt,bpt):
    elemL = 0
    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        coor = linElastStat.parametricCoordinate(apt[0],bpt[0],apt[1],bpt[1],gausspoints[qj],gausspoints[qj])
        gcoor = linElastStat.geometricCoordinate(coor,U,V,w,p,q,px,py)
        # print("Geometric coor")
        # print(gcoor)
        # Extracting the unique non-zero value
        jac2 = 0.5*float(np.extract((bpt-apt) > 1e-5, (bpt-apt)))
        du = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])
        dxdu = du@px
        dydu = du@py
        jac1 = np.sqrt(dxdu**2 + dydu**2)
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
        K += linElastStat.localStiffnessMatrix(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint,dmat)
        Fb += linElastStat.localBodyVector(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint,rho)
        totalArea += elementArea(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint)
        if ielem in loadelems:
            print('Loaded element')

            uB = paramnodes[nodeselem[ielem][2]][0]
            uA = paramnodes[nodeselem[ielem][3]][0]
            vB = paramnodes[nodeselem[ielem][2]][1]
            vA = paramnodes[nodeselem[ielem][3]][1]

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
    for rdof in remdofs:
        Kcol = Kred[:,rdof]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        Fred -= Kcol*ucond

    print("Second reduction")
    Kred = np.delete(Kred,remdofs,1)

    return Kred,Fred,remdofs

def solveMatrixEquations(Kred,Fred,totaldofs,remdofs):
    # Checking full rank in matrix
    mRank = np.linalg.matrix_rank(Kred)
    mRows = Kred.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        fullRank = True
    else:
        fullRank = False

    # fullRank = True
    if fullRank:
        print("The matrix has full rank. It is invertible")
        dred = np.linalg.solve(Kred,Fred)
    else:
        print("The matrix has not full rank. It is not invertible")
        dred = np.linalg.lstsq(Kred,Fred,rcond=None)[0]

    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    dtotal = np.zeros((totaldofs.shape[0],1))
    dtotal[reduceddofs,:] = dred
    return dtotal,dred

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

# Pinit = np.array([[1,0,0],[1,np.sqrt(2)-1,0],[np.sqrt(2)-1,1,0],[0,1,0],[2.5,0,0],[2.5,0.75,0],
#               [0.75,2.5,0],[0,2.5,0],[4,0,0],[4,4,0],[4,4,0],[0,4,0]])

winit = np.array([[1],[0.5*(1+(1/np.sqrt(2)))],[0.5*(1+(1/np.sqrt(2)))],[1],[1],[1],
              [1],[1],[1],[1],[1],[1]])

#Isogeometric routines
Uinit = np.array([0,0,0,0.5,1,1,1])
Vinit = np.array([0,0,0,1,1,1])

pinit = 2
qinit = 2

doRefinement = 'N'

if doRefinement == 'Y':
    reflist = ['h','h','h','h']
    dirlist = ['U','V','U','V']
    Uinp,Vinp,pinp,qinp,Pinp,winp = srfn.surfaceRefinement(reflist,dirlist,Uinit,Vinit,pinit,qinit,Pinit,winit)
else:
    Uinp = Uinit
    Vinp = Vinit
    pinp = pinit
    qinp = qinit
    Pinp = Pinit
    winp = winit

displacementConditions = [[0.0,0,"S"],[0.0,1,"S"]]
neumannConditions = [[[0.0,1.0],[0.5,1.0]]]

parametricNodes,nodesInElement = pre2D.parametricGrid(Uinp,Vinp)
loadNodes,loadElements = pre2D.loadPreprocessing(parametricNodes,nodesInElement,Uinp,Vinp,pinp,qinp,Pinp,winp,-4.0)
print(loadNodes)
print(loadElements)
loadElements,loadFaces = pre2D.loadPreprocessingv2(parametricNodes,nodesInElement,neumannConditions)
print(loadElements)
print(loadFaces)
dirichletCtrlPts,axisRestrictions = pre2D.dirichletBCPreprocessing(Pinp,displacementConditions)

# cx,cy = rbs.nurbs2DField(Uinit,Vinit,pinp,qinp,Pinp,winp)
# plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])),Pinp)

# dMat = linElastStat.elasticMatrix(E,nu)
# K,F,totalArea = assemblyWeakForm(Uinp,Vinp,winp,pinp,qinp,Pinp,parametricNodes,nodesInElement,gaussLegendreQuadrature,dMat,rho,loadNodes,loadElements,tv)

# Kred,Fred,removedDofs = boundaryConditionsEnforcement(K,F,dirichletCtrlPts,axisRestrictions,u0)

# totaldofs = np.arange(2*Pinp.shape[0])
# dtotal,dred = solveMatrixEquations(Kred,Fred,totaldofs,removedDofs)

# dx = dtotal[0::2]
# dy = dtotal[1::2]
# D = np.hstack((dx,dy))
# print(D)

# post2D.postProcessing(Uinp,Vinp,pinp,qinp,Pinp,D,winp,parametricNodes,nodesInElement,dtotal,dMat)
