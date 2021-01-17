import numpy as np
import nurbs as rbs

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

# Improve this function
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

def loadPreprocessingv2(paramnodes,nodeselem,neumannconditions):
    loadelements = []
    loadfaces = []

    for cond in neumannconditions:
        apt = np.array(cond[0])
        bpt = np.array(cond[1])

        for ielem in range(0,nodeselem.shape[0]):
            nodecounter = 0

            print("Element #",ielem)

            for j in range(0,4):
                xpt = paramnodes[nodeselem[ielem][j]]

                print('----')
                print(apt)
                print(bpt)
                print(xpt)
                print('----')

                isAInSquare = ( (apt[0] - xpt[0])**2 + (apt[1] - xpt[1])**2 ) < 1e-4
                isBInSquare = ( (bpt[0] - xpt[0])**2 + (bpt[1] - xpt[1])**2 ) < 1e-4

                if isAInSquare or isBInSquare:
                    print("Node In")
                    nodecounter += 1

                if nodecounter == 2:
                    print("Element In")
                    loadelements.append(ielem)
                    loadfaces.append(j-1)

    return loadelements,loadfaces

def dirichletBCPreprocessing(P,dirichletconditions):
    dirichletctrlpts = []
    axisrestrictions = []

    for cond in dirichletconditions:
        cdirichlet = cond[0]
        axis = cond[1]
        restriction = cond[2]

        for i in range(0,P.shape[0]):
            if abs(P[i][axis] - cdirichlet) < 1e-4:
                dirichletctrlpts.append(i+1)

                if restriction == "C":
                    # Due to clamped condition
                    axisrestrictions.append(0)
                    axisrestrictions.append(1)
                else:
                    # Supported condition assumed
                    axisrestrictions.append(axis)

    return dirichletctrlpts,axisrestrictions
