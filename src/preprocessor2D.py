# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

# Local project
import src.nurbs as rbs
import src.plottingScripts as plts

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

def checkColinearPoints(apt,bpt,cpt):
    abdist = np.sqrt( (apt[0] - bpt[0])**2 + (apt[1] - bpt[1])**2 )
    acdist = np.sqrt( (apt[0] - cpt[0])**2 + (apt[1] - cpt[1])**2 )
    cbdist = np.sqrt( (cpt[0] - bpt[0])**2 + (cpt[1] - bpt[1])**2 )

    # Check if AB = AC + CB
    if abs(abdist - acdist - cbdist) < 1e-5:
        return True
    else:
        return False

def loadPreprocessing(paramnodes,nodeselem,neumannconditions):
    loadednodes = []
    loadelements = []
    loadfaces = []

    for cond in neumannconditions:
        apt = np.array(cond[0])
        bpt = np.array(cond[1])
        loadtype = cond[2]
        # print(loadtype)

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

def dirichletBCPreprocessingOnFaces(P,dirichletconditions):
    dirichletconds = []

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    for cond in dirichletconditions:
        apt = np.array(cond[0])
        bpt = np.array(cond[1])
        restriction = cond[2]
        value = cond[3]

        for i in range(0,P.shape[0]):
            nodecond = []

            # Check the control points that belong to the
            # Dirichlet boundary condition
            if checkColinearPoints(apt,bpt,P[i,:]):
                # Insertion of the node
                nodecond.append(i)
                # Insertion of the type of restriction
                nodecond.append(restriction)

                if restriction == "C":
                    # Clamped condition assumed
                    # Insertion of the axis where the condition are applied
                    nodecond.append([0,1])
                else:
                    # Supported condition assumed
                    tangetVec = bpt - apt
                    tangetVec = np.reshape(tangetVec,(len(tangetVec),1))
                    unitTangetVec = tangetVec/np.linalg.norm(tangetVec)

                    unitNormalVec = rotMat@unitTangetVec

                    for i in range(len(unitNormalVec)):
                        if abs(unitNormalVec[i]) > 1e-5:
                            axis = i

                    # Insertion of the axis where the condition is applied
                    nodecond.append([axis])

                # Insertion of the enforced value
                nodecond.append(value)

                # Insertion of the whole package of conditions
                dirichletconds.append(nodecond)

    return dirichletconds

def plotGeometry(U,V,p,q,P,w,dirichletconds,neumannconditions,paramnodes,nodeselem,loadelements,loadfaces):
    fig = plt.figure()
    ax = plt.axes()
    plt.axis('equal')
    ax.use_sticky_edges = False
    titlestring = "Geometry with boundary conditions"
    ax.axis("off")

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

    jmat = np.zeros((2,2))
    numpt = 5
    loadcoor = []
    loadfield = []

    loadtype = neumannconditions[0][2]
    loadvalue = neumannconditions[0][3]

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    # Boundary Geometry
    cxb,cyb = rbs.nurbs2DBoundary(U,V,p,q,P,w)
    fieldplot = ax.fill(cxb,cyb,facecolor='none',edgecolor='black',linewidth=1.5)

    # Control Points
    # ax.scatter(P[:,0],P[:,1])

    # Dirichlet Conditions
    dirctrlpts = []
    for drchcond in dirichletconds:
        inode = drchcond[0]
        dirctrlpts.append(inode)

    for i in range(0,len(dirichletconds)):
        if dirichletconds[i][1] == "C":
            dirichletplot = ax.scatter(P[dirctrlpts,0],P[dirctrlpts,1],c = "r",marker = "^")
        else:
            dirichletplot = ax.scatter(P[dirctrlpts,0],P[dirctrlpts,1],c = "g",marker = "o")

    # Neumann Conditions
    for ielem in range(0,len(loadelements)):

        paramside = loadfaces[ielem]
        if paramside == 0 or paramside == 1:
            startindex = paramside
            endindex = paramside + 1
        elif paramside == 2:
            startindex = paramside + 1
            endindex = paramside
        else:
            startindex = 3
            endindex = 0

        if paramside == 1 or paramside == 3:
            paramaxis = 1
        else:
            paramaxis = 0

        uB = paramnodes[nodeselem[loadelements[ielem]][endindex]][0]
        uA = paramnodes[nodeselem[loadelements[ielem]][startindex]][0]
        vB = paramnodes[nodeselem[loadelements[ielem]][endindex]][1]
        vA = paramnodes[nodeselem[loadelements[ielem]][startindex]][1]

        startpt = np.array([uA,vA])
        endpt = np.array([uB,vB])

        parampath = np.linspace(startpt,endpt,numpt,endpoint=True)
        geomcoor = np.zeros((parampath.shape[0],2))
        fieldpatch = np.zeros((parampath.shape[0],2))
        ipath = 0
        for ppath in parampath:
            ratFunc,dN2du,dN2dv = rbs.rationalFunctionAndGradient(U,V,w,p,q,ppath[0],ppath[1])

            dxdu = dN2du@px
            dxdv = dN2dv@px
            dydu = dN2du@py
            dydv = dN2dv@py

            jmat[0][0] = dxdu
            jmat[0][1] = dxdv
            jmat[1][0] = dydu
            jmat[1][1] = dydv

            jvec = jmat[:,paramaxis]
            normjvec = np.linalg.norm(jvec)

            if normjvec > 1e-6:
                unitTangetVec = jvec/normjvec
            else:
                unitTangetVec = np.zeros((2,1))

            if loadtype == "tangent":
                loadvec = (loadvalue/abs(loadvalue))*unitTangetVec
            elif loadtype == "normal":
                unitNormalVec = rotMat@unitTangetVec
                loadvec = (loadvalue/abs(loadvalue))*unitNormalVec
            else:
                print("Wrong load configuration")

            # ratFunc = rbs.ratFunction(U,V,w,p,q,ppath[0],ppath[1])
            geomcoor[ipath][0] = ratFunc@px
            geomcoor[ipath][1] = ratFunc@py
            fieldpatch[ipath,:] = loadvec.T
            ipath += 1

        loadcoor.append(geomcoor)
        loadfield.append(fieldpatch)

    for ld in range(len(loadcoor)):
        if ld == 0:
            loadcoor1 = loadcoor[ld]
            loadfield1 = loadfield[ld]
        else:
            loadcoor1 = np.vstack((loadcoor1,loadcoor[ld]))
            loadfield1 = np.vstack((loadfield1,loadfield[ld]))

    # neumannplot = ax.scatter(loadcoor1[:,0],loadcoor1[:,1],c = "b",marker = "s")
    neumannplot = ax.quiver(loadcoor1[:,0],loadcoor1[:,1],loadfield1[:,0],loadfield1[:,1],color=['b'])

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(titlestring)
    # Uncomment for example2
    plt.legend((dirichletplot,neumannplot),('Displacement restrictions','Load conditions'),loc='upper right',bbox_to_anchor=(1.2,1.0))
    # Uncomment for example3
    # plt.legend((dirichletplot,neumannplot),('Displacement restrictions','Load conditions'),loc='lower right',bbox_to_anchor=(1.2,0.0))
    plt.tight_layout()
    plt.show()