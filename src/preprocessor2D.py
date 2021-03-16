# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs
import src.plottingScripts as plts

################ PREPROCESSING ####################

def parametricGrid(U,V,p,q):
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

    """
    - The diagonals of a given ABCD parametric element are given by points C and A
      which values can be found as paramnodes[C][0] and paramnodes[A][0] respectively.
      Each corner is organized in the element as follows:
      ielem -> [A B C D]
                0 1 2 3
     """

    for j in range(0,numV-1):
        for i in range(0,numU-1):
            elemmat[elemIndex,0] = j*numU + i #A
            elemmat[elemIndex,1] = j*numU + (i+1) #B
            elemmat[elemIndex,2] = (j+1)*numU + (i+1) #C
            elemmat[elemIndex,3] = (j+1)*numU + i #D

            elemIndex += 1

    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1

    nonzeroctrlpts = []
    surfacespan = []
    elementcorners = []

    """
    - C can be obtained as the 3rd element in a given row of the matrix nodeselem
    - Likewise, A is the 1st element in a ielem row of nodeselem
    - It means that

      uC = paramnodes[C][0]
      uA = paramnodes[A][0]
      vC = paramnodes[C][1]
      vA = paramnodes[A][1]

      where
      C = elemmat[ielem][2]
      A = element[ielem][0]
      [0] stands for the u component
      [1] stands for the v component
    """

    for ielem in range(0,elemIndex):
        uC = paramnodes[elemmat[ielem][2]][0]
        uA = paramnodes[elemmat[ielem][0]][0]
        vC = paramnodes[elemmat[ielem][2]][1]
        vA = paramnodes[elemmat[ielem][0]][1]

        apt = np.array([uA,vA])
        cpt = np.array([uC,vC])

        uspan = bfunc.findKnotInterval(nu,p,0.5*(uC + uA),U)
        vspan = bfunc.findKnotInterval(nv,q,0.5*(vC + vA),V)

        idR = rbs.nonZeroIndicesSurface(uspan,vspan,p,q,nu)

        nonzeroctrlpts.append(idR)
        surfacespan.append([uspan,vspan])
        elementcorners.append([apt,cpt])

    surfaceprep = [nonzeroctrlpts,surfacespan,elementcorners]

    return paramnodes,elemmat,surfaceprep

def checkColinearPoints(apt,bpt,cpt):
    abdist = np.sqrt( (apt[0] - bpt[0])**2 + (apt[1] - bpt[1])**2 )
    acdist = np.sqrt( (apt[0] - cpt[0])**2 + (apt[1] - cpt[1])**2 )
    cbdist = np.sqrt( (cpt[0] - bpt[0])**2 + (cpt[1] - bpt[1])**2 )

    # Check if AB = AC + CB
    if abs(abdist - acdist - cbdist) < 1e-5:
        return True
    else:
        return False

def loadPreprocessing(paramnodes,nodeselem,neumannconditions,U,V,p,q):
    loadednodes = []
    loadfaces = []

    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1

    loadelements = []
    boundaryspan = []
    boundarycorners = []
    axisselector = []
    nonzeroctrlptsload = []

    for cond in neumannconditions:
        startpt = np.array(cond[0])
        endpt = np.array(cond[1])
        loadtype = cond[2]

        # Check the other parametric nodes that belong to the
        # neumann boundary condition
        for inode in range(0,paramnodes.shape[0]):
            if checkColinearPoints(startpt,endpt,paramnodes[inode,:]):
                loadednodes.append(inode)

        # Check the nodes that have applied load and match them
        # with the list of parametric elements
        for ielem in range(0,nodeselem.shape[0]):
            commom_nodes = set.intersection(set(loadednodes),set(nodeselem[ielem,:]))
            if len(commom_nodes) == 2:
                loadelements.append(ielem)

        # Check the sides of the loaded parametric elements that
        # belong t the neumann boundary condition
        for ldnelm in loadelements:
            for j in range(0,4):
                if j < 3:
                    side_nodes = [nodeselem[ldnelm][j],nodeselem[ldnelm][j+1]]
                else:
                    side_nodes = [nodeselem[ldnelm][j],nodeselem[ldnelm][0]]
                face = set.intersection(set(side_nodes),set(loadednodes))
                if len(face) == 2:
                    loadfaces.append(j)


    for iface in range(0,len(loadfaces)):
        elemface = loadfaces[iface]
        ielem = loadelements[iface]

        # Selecting the axis where the jacobian
        # is not a zero vector on the boundary
        if elemface == 1 or elemface == 3:
           paramaxis = 1
        else:
           paramaxis = 0

        if elemface == 0 or elemface == 1:
            startindex = elemface
            endindex = elemface + 1
        elif elemface == 2:
            startindex = elemface + 1
            endindex = elemface
        else:
            startindex = 3
            endindex = 0

        uB = paramnodes[nodeselem[ielem][endindex]][0]
        uA = paramnodes[nodeselem[ielem][startindex]][0]
        vB = paramnodes[nodeselem[ielem][endindex]][1]
        vA = paramnodes[nodeselem[ielem][startindex]][1]

        # Computing the corners of the segment
        apt = np.array([uA,vA])
        bpt = np.array([uB,vB])

        uspan = bfunc.findKnotInterval(nu,p,0.5*(uB + uA),U)
        vspan = bfunc.findKnotInterval(nv,q,0.5*(vB + vA),V)

        idR = rbs.nonZeroIndicesSurface(uspan,vspan,p,q,nu)

        boundarycorners.append([apt,bpt])
        boundaryspan.append([uspan,vspan])
        axisselector.append(paramaxis)
        nonzeroctrlptsload.append(idR)

    boundaryprep = [nonzeroctrlptsload,boundaryspan,boundarycorners,axisselector]

    return boundaryprep

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

                    for j in range(len(unitNormalVec)):
                        if abs(unitNormalVec[j]) > 1e-5:
                            axis = j

                    # Insertion of the axis where the condition is applied
                    nodecond.append([axis])

                # Insertion of the enforced value
                nodecond.append(value)

                # Insertion of the whole package of conditions
                dirichletconds.append(nodecond)

    return dirichletconds

def numericalIntegrationPreprocessing(numgauss):
    numericalquadrature = np.polynomial.legendre.leggauss(numgauss)

    quadraturepoints = numericalquadrature[0]
    quadratureweights = numericalquadrature[1]

    gausslegendre2d = np.zeros((numgauss*numgauss,3))
    gausslegendre1d = np.zeros((numgauss,2))

    igauss = 0
    for i in range(numgauss):
        gausslegendre1d[igauss][0] = quadraturepoints[i]
        gausslegendre1d[igauss][1] = quadratureweights[i]
        igauss += 1

    igauss = 0
    for j in range(numgauss):
        for i in range(numgauss):
            gausslegendre2d[igauss][0] = quadraturepoints[i]
            gausslegendre2d[igauss][1] = quadraturepoints[j]
            gausslegendre2d[igauss][2] = quadratureweights[i]*quadratureweights[j]
            igauss += 1

    numericalquad = [gausslegendre2d,gausslegendre1d]

    return numericalquad

def problemPreprocessing(surface,dirichletconditions,neumannconditions):
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()
    parametricNodes,nodesInElement,surfacePreprocessing = parametricGrid(U,V,p,q)
    boundaryPreprocessing = loadPreprocessing(parametricNodes,nodesInElement,neumannconditions,U,V,p,q)
    dirichletBCList = dirichletBCPreprocessingOnFaces(P,dirichletconditions)
    return surfacePreprocessing,boundaryPreprocessing,dirichletBCList

def plotGeometry(surface,dirichletconds,neumannconditions,boundaryprep):
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()

    fig = plt.figure()
    ax = plt.axes()
    plt.axis('equal')
    ax.use_sticky_edges = False
    titlestring = "Geometry with boundary conditions"
    ax.axis("off")

    jmat = np.zeros((2,2))
    numpt = 5
    loadcoor = []
    loadfield = []

    loadtype = neumannconditions[0][2]
    loadvalue = neumannconditions[0][3]

    # Extraction of boundary preprocessing
    nonzeroctrlptsload = boundaryprep[0]
    boundaryspan = boundaryprep[1]
    boundarycorners = boundaryprep[2]
    axisselector = boundaryprep[3]

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    # Boundary Geometry
    cbpts = rbs.nurbs2DBoundary(U,V,p,q,P,w)
    fieldplot = ax.fill(cbpts[:,0],cbpts[:,1],facecolor='none',edgecolor='black',linewidth=1.5)

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

    mu = len(U) - 1
    mv = len(V) - 1

    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)

    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)

    # Neumann Conditions
    for iload in range(0,len(boundarycorners)):
        # Extracting the indices of the non-zero control points of the loaded elements
        idR = nonzeroctrlptsload[iload]
        # Extracting the indices for the location of the parametric segment
        uspan = boundaryspan[iload][0]
        vspan = boundaryspan[iload][1]
        # Extracting the corners of the parametric segment
        startpt = boundarycorners[iload][0]
        endpt = boundarycorners[iload][1]
        # Extracting the non-zero column index of the boundary jacobian
        paramaxis = axisselector[iload]

        parampath = np.linspace(startpt,endpt,numpt,endpoint=True)
        geomcoor = np.zeros((parampath.shape[0],2))
        fieldpatch = np.zeros((parampath.shape[0],2))
        ipath = 0
        for ppath in parampath:
            biRatGrad = rbs.bivariateRationalGradient(mu,mv,p,q,uspan,vspan,ppath[0],ppath[1],U,V,Pw)

            jmat = (biRatGrad[1:3,:]@P[idR,:]).T

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

            geomcoor[ipath,:] = biRatGrad[0,:]@P[idR,:]
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
