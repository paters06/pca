# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Local project
import src.nurbs as rbs
import src.linearElastoStaticsSolver as linElastStat

def checkingSymmetricMatrix(A):
    check = np.allclose(A, A.T, rtol=1e-2, atol=1e-3)
    if check:
        print("The given matrix is symmetric")
    else:
        print("The given matrix is not symmetric")

def checkFullRankMatrix(A):
    mRank = np.linalg.matrix_rank(A)
    mRows = A.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        return True
    else:
        return False

def plotSparsity(A):
    fig = plt.figure()
    plt.spy(A,markersize=5)
    # plt.imshow(A,cmap=cm.viridis)
    # plt.colorbar()
    plt.show()

################ AREA INTEGRALS ####################

def calculateArea(surface,surfaceprep,numquad):
    totalArea = 0
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()

    numquad2d = numquad[0]
    numquad1d = numquad[1]

    # Extraction of surface preprocessing
    nonzeroctrlpts = surfaceprep[0]
    surfacespan = surfaceprep[1]
    elementcorners = surfaceprep[2]

    # Extraction of boundary preprocessing
    # nonzeroctrlptsload = boundaryprep[0]
    # boundaryspan = boundaryprep[1]
    # boundarycorners = boundaryprep[2]
    # axisselector = boundaryprep[3]

    paramGrad = np.zeros((2,2))
    numElems = len(elementcorners)
    # numLoadedElems = len(boundarycorners)

    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)

    # Precomputing info for the nurbs derivatives
    mU = len(U) - 1
    mV = len(V) - 1
    nU = mU - p - 1
    nV = mV - q - 1

    # Strain-energy and body force integrals
    print('Computing the strain-energy and body forces integrals')
    for ielem in range(0,numElems):
        # Extracting the indices of the non-zero control points
        idR = nonzeroctrlpts[ielem]
        # Extracting the indices for the location of the parametric element
        uspan = surfacespan[ielem][0]
        vspan = surfacespan[ielem][1]
        # Extracting the corners of the parametric element
        aPoint = elementcorners[ielem][0]
        cPoint = elementcorners[ielem][1]

        # Computing the parametric gradient and its determinant
        paramGrad[0][0] = 0.5*(cPoint[0] - aPoint[0])
        paramGrad[1][1] = 0.5*(cPoint[1] - aPoint[1])
        detJac2 = abs(np.linalg.det(paramGrad))

        # K stiffness matrix
        for iquad in range(numquad2d.shape[0]):
            coor = linElastStat.parametricCoordinate(aPoint[0],cPoint[0],aPoint[1],cPoint[1],numquad2d[iquad][0],numquad2d[iquad][1])

            # NURBS gradient
            biRatGrad = rbs.bivariateRationalGradient(mU,mV,p,q,uspan,vspan,coor[0][0],coor[0][1],U,V,Pw)

            # Jacobian
            jac = (biRatGrad[1:3,:]@P[idR,:]).T
            wJac = abs(np.linalg.det(jac))*detJac2*numquad2d[iquad][2]

            totalArea += wJac

        # End iquad loop
    # End ielem loop

    print("Total Area")
    print(totalArea)
    # return totalArea
