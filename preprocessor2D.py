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

def checkColinearPoints(apt,bpt,cpt):
    abdist = np.sqrt( (apt[0] - bpt[0])**2 + (apt[1] - bpt[1])**2 )
    acdist = np.sqrt( (apt[0] - cpt[0])**2 + (apt[1] - cpt[1])**2 )
    cbdist = np.sqrt( (cpt[0] - bpt[0])**2 + (cpt[1] - bpt[1])**2 )

    # Check if AB = AC + CB
    if abs(abdist - acdist - cbdist) < 1e-5:
        return True
    else:
        return False

def loadPreprocessingv2(paramnodes,nodeselem,neumannconditions):
    loadednodes = []
    loadelements = []
    loadfaces = []

    for cond in neumannconditions:
        apt = np.array(cond[0])
        bpt = np.array(cond[1])

        # Check the other parametric nodes that belong to the
        # neumann boundary condition
        for inode in range(0,paramnodes.shape[0]):
            if checkColinearPoints(apt,bpt,paramnodes[inode,:]):
                loadednodes.append(inode)

        # print(loadednodes)
        for ielem in range(0,nodeselem.shape[0]):
            # print(nodeselem[ielem,:])
            commom_nodes = set.intersection(set(loadednodes),set(nodeselem[ielem,:]))
            # commom_nodes = list(commom_nodes)
            # print(len(commom_nodes))
            if len(commom_nodes) == 2:
                loadelements.append(ielem)

        # print(loadelements)
        for ldnelm in loadelements:
            for j in range(0,4):
                if j < 3:
                    side_nodes = [nodeselem[ldnelm][j],nodeselem[ldnelm][j+1]]
                else:
                    side_nodes = [nodeselem[ldnelm][j],nodeselem[ldnelm][0]]
                face = set.intersection(set(side_nodes),set(loadednodes))
                if len(face) == 2:
                    loadfaces.append(j)
                    # print(j)

        # print(loadfaces)

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
