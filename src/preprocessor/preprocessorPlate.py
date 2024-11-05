# Python libraries
import numpy as np
from numpy.lib.arraysetops import unique
import numpy.linalg
import matplotlib.pyplot as plt

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs
import src.plotting.plottingScripts as plts

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
    paramknots = np.hstack((ucomp,vcomp))

    # Assembling the element matrix
    numU = len(uniqueU)
    numV = len(uniqueV)
    numelems = (numU - 1)*(numV - 1)

    elemmat = np.zeros((numelems,4),dtype=int)
    elemIndex = 0

    """
    - The diagonals of a given ABCD parametric element are given by points C and A
      which values can be found as paramknots[C][0] and paramknots[A][0] respectively.
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
        # End for loop
    # End for loop

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

      uC = paramknots[C][0]
      uA = paramknots[A][0]
      vC = paramknots[C][1]
      vA = paramknots[A][1]

      where
      C = elemmat[ielem][2]
      A = element[ielem][0]
      [0] stands for the u component
      [1] stands for the v component
    """

    for ielem in range(0,elemIndex):
        uC = paramknots[elemmat[ielem][2]][0]
        uA = paramknots[elemmat[ielem][0]][0]
        vC = paramknots[elemmat[ielem][2]][1]
        vA = paramknots[elemmat[ielem][0]][1]

        apt = np.array([uA,vA])
        cpt = np.array([uC,vC])

        uspan = bfunc.findKnotInterval(nu,p,0.5*(uC + uA),U)
        vspan = bfunc.findKnotInterval(nv,q,0.5*(vC + vA),V)

        idR = rbs.nonZeroIndicesElement(uspan,vspan,p,q,nu)

        nonzeroctrlpts.append(idR)
        surfacespan.append([uspan,vspan])
        elementcorners.append([apt,cpt])
    # End for loop

    surfaceprep = [nonzeroctrlpts,surfacespan,elementcorners]

    return paramknots,elemmat,surfaceprep
# End function

def neumannBCPreprocessing_Plate(neumannconditionsdata,U,V,p,q):
    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1

    neumann_element = []
    neumann_localside = []
    neumann_axisselector = []
    neumann_type = []
    neumann_value = []

    uniqueU = np.unique(U)
    uniqueV = np.unique(V)

    numSegmentsU = len(uniqueU) - 1
    numSegmentsV = len(uniqueV) - 1

    for cond in neumannconditionsdata:
        startpt = np.array(cond[0])
        endpt = np.array(cond[1])

        segmentpt = endpt - startpt

        if abs(segmentpt[0]) < 1e-5:
            uspan = bfunc.findKnotInterval(nu,p,startpt[0],U)
            for iv in range(len(uniqueV)-1):
                vspan = bfunc.findKnotInterval(nv,q,uniqueV[iv],V)
                id_elem = numSegmentsU*(vspan - 1) + uspan
                neumann_element.append(id_elem)
                
                if abs(startpt[0]) < 1e-5:
                    neumann_localside.append(3)
                    neumann_axisselector.append(1)
                elif abs(startpt[0] - 1.0) < 1e-5:
                    neumann_localside.append(1)
                    neumann_axisselector.append(1)
                # End if

                neumann_type.append(cond[2])
                neumann_value.append(cond[3])
            # End iv for loop
        elif abs(segmentpt[1]) < 1e-5:
            vspan = bfunc.findKnotInterval(nv,q,startpt[1],V)
            for iu in range(len(uniqueU)-1):
                uspan = bfunc.findKnotInterval(nu,p,uniqueU[iu],U)
                id_elem = numSegmentsU*(vspan - 1) + uspan
                neumann_element.append(id_elem)
                
                if abs(startpt[1]) < 1e-5:
                    neumann_localside.append(0)
                    neumann_axisselector.append(0)
                elif abs(startpt[1] - 1.0) < 1e-5:
                    neumann_localside.append(2)
                    neumann_axisselector.append(0)
                # End if

                neumann_type.append(cond[2])
                neumann_value.append(cond[3])
            # End iu for loop
        else:
            print("Wrong dimensions")
        # End if

    # End for loop
    return [neumann_element,neumann_localside,neumann_axisselector,neumann_type,neumann_value]
# End function

def dirichletBCPreprocessing_Plate(dirichletconditions,U,V,p,q):
    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1

    enforcedctrlpts = []
    enforceddof = []
    enforcedvalues = []

    numdof = 3

    uniqueU = np.unique(U)
    uniqueV = np.unique(V)

    for cond in dirichletconditions:
        apt = np.array(cond[0])
        bpt = np.array(cond[1])
        value = cond[2]
        restriction = cond[3]

        cpt = bpt - apt
        id_ctrlpts = []

        # uspan = bfunc.findKnotInterval(nu,p,cpt[0],U)
        # vspan = bfunc.findKnotInterval(nv,q,cpt[1],V)

        if abs(cpt[0]) < 1e-5:
            uspan = bfunc.findKnotInterval(nu,p,apt[0],U)
            for iv in range(len(uniqueV)-1):
                vspan = bfunc.findKnotInterval(nv,q,uniqueV[iv],V)
                idr_bound = rbs.nonZeroIndiceBoundary(apt,bpt,uspan,vspan,p,q,nu)
                id_ctrlpts += idr_bound
            # End iv for loop
        elif abs(cpt[1]) < 1e-5:
            vspan = bfunc.findKnotInterval(nv,q,apt[1],V)
            for iu in range(len(uniqueU)-1):
                uspan = bfunc.findKnotInterval(nu,p,uniqueU[iu],U)
                idr_bound = rbs.nonZeroIndiceBoundary(apt,bpt,uspan,vspan,p,q,nu)
                id_ctrlpts += idr_bound
            # End iu for loop
        else:
            print("Wrong dimensions")
        # End if

        id_ctrlpts_arr = np.array(id_ctrlpts)
        id_ctrlpts = list(np.unique(id_ctrlpts_arr))

        for i in id_ctrlpts:
            enforcedctrlpts.append([i,restriction])
            if restriction == "C":
                for ii in range(0,numdof):
                    enforceddof.append(numdof*i + ii)
                    enforcedvalues.append(value)
                # End ii for loop
            elif restriction == "SS1":
                enforceddof.append(numdof*i)
                enforcedvalues.append(value)
            elif restriction == "SS2":
                for ii in range(0,numdof-1):
                    enforceddof.append(numdof*i + ii)
                    enforcedvalues.append(value)
                # End ii for loop
            # End if
        # End i for loop
    # End loop
    # print(enforcedctrlpts)
    # print(enforceddof)

    return enforcedctrlpts,enforceddof,enforcedvalues
# End function

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
    # End for loop

    igauss = 0
    for j in range(numgauss):
        for i in range(numgauss):
            gausslegendre2d[igauss][0] = quadraturepoints[i]
            gausslegendre2d[igauss][1] = quadraturepoints[j]
            gausslegendre2d[igauss][2] = quadratureweights[i]*quadratureweights[j]
            igauss += 1
        # End for loop
    # End for loop

    numericalquad = [gausslegendre2d,gausslegendre1d]

    return numericalquad
# End function

def problemPreprocessing(phenomenon,surface,dirichletconditionsdata,neumannconditionsdata):
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()
    parametricNodes,nodesInElement,surfacePreprocessing = parametricGrid(U,V,p,q)
    if neumannconditionsdata is not None:
        # boundaryPreprocessing = neumannBCPreprocessing(parametricNodes,nodesInElement,neumannconditionsdata,U,V,p,q)
        boundaryPreprocessing = neumannBCPreprocessing_Plate(neumannconditionsdata,U,V,p,q)
        print(boundaryPreprocessing[0])
    else:
        boundaryPreprocessing = None
    # End if

    if phenomenon == "Elasticity":
        # dirichletBCList,enforcedDOF,enforcedValues = dirichletBCPreprocessing_Elasticity(P,surface,dirichletconditionsdata,U,V,p,q)
        # print(enforcedDOF)
        enforcedCtrlPts,enforcedDOF,enforcedValues = dirichletBCPreprocessing_Plate(dirichletconditionsdata,U,V,p,q)
    # End if
    if phenomenon == "Heat" or phenomenon == "Schrodinger":
        enforcedCtrlPts,enforcedDOF,enforcedValues = dirichletBCPreprocessing_Plate(dirichletconditionsdata,U,V,p,q)
    # End if

    return surfacePreprocessing,boundaryPreprocessing,enforcedCtrlPts,enforcedDOF,enforcedValues
# End function

def plotGeometry(phenomenon,surface,dirichletconds,boundaryprep):
    if phenomenon == "Elasticity":
        first_string = "Displacement BC"
        second_string = "Load BC"
    elif phenomenon == "Heat":
        first_string = "Temperature BC"
        second_string = "Flux BC"
    elif phenomenon == "Schrodinger":
        first_string = "Wave Function BC"
        second_string = "Derivative Wave Function BC"
    else:
        first_string = "Null BC"
        second_string = "Null BC"

    U,V,p,q,P,w = surface.retrieveSurfaceInformation()

    fig = plt.figure()
    ax = plt.axes()
    plt.axis('equal')
    ax.use_sticky_edges = False
    titlestring = "Geometry with boundary conditions"
    ax.axis("off")
    # plt.xticks([])

    jmat = np.zeros((2,2))
    numpt = 5
    loadcoor = []
    loadfield = []

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    # Boundary Geometry
    cbpts = surface.createBoundary()
    fieldplot = ax.fill(cbpts[:,0],cbpts[:,1],facecolor='none',edgecolor='black',linewidth=1.5)

    # Control Points
    ax.scatter(P[:,0],P[:,1])

    # Dirichlet Conditions
    dirctrlpts = []
    for drchcond in dirichletconds:
        inode = drchcond[0]
        dirctrlpts.append(inode)

    for i in range(0,len(dirichletconds)):
        # print(dirichletconds[i])
        if len(dirichletconds[i]) == 2:
            dirichletplot = ax.scatter(P[dirctrlpts,0],P[dirctrlpts,1],c = "g",marker = "o")
        else:
            dofs = dirichletconds[i][2]
            if len(dofs) == 2:
                dirichletplot = ax.scatter(P[dirctrlpts,0],P[dirctrlpts,1],c = "r",marker = "^")
            else:
                dirichletplot = ax.scatter(P[dirctrlpts,0],P[dirctrlpts,1],c = "g",marker = "o")
            # End if
        # End if
    # End for loop

    mu = len(U) - 1
    mv = len(V) - 1

    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)

    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)

    if boundaryprep is not None:
        # Extraction of boundary preprocessing
        nonzeroctrlpts_boundary = boundaryprep[0]
        boundaryspan = boundaryprep[1]
        boundarycorners = boundaryprep[2]
        axisselector = boundaryprep[3]
        neumannbc_type = boundaryprep[4]
        neumannbc_value = boundaryprep[5]

        # Neumann Conditions
        for ineumann in range(0,len(boundarycorners)):
            # Extracting the indices of the non-zero control points of the loaded elements
            idR = nonzeroctrlpts_boundary[ineumann]
            # Extracting the indices for the location of the parametric segment
            uspan = boundaryspan[ineumann][0]
            vspan = boundaryspan[ineumann][1]
            # Extracting the corners of the parametric segment
            startpt = boundarycorners[ineumann][0]
            endpt = boundarycorners[ineumann][1]
            # Extracting the non-zero column index of the boundary jacobian
            paramaxis = axisselector[ineumann]
            # Extracting the type and value of the neumann conditions
            neumanntype = neumannbc_type[ineumann]
            neumannval = neumannbc_value[ineumann]

            if abs(neumannval) > 1e-5:
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
                        unitTangentVec = jvec/normjvec
                    else:
                        unitTangetVec = np.zeros((2,1))

                    if neumanntype == "tangent":
                        neumannvec = (neumannval/abs(neumannval))*unitTangentVec
                    elif neumanntype == "normal":
                        unitNormalVec = rotMat@unitTangentVec
                        neumannvec = (neumannval/abs(neumannval))*unitNormalVec
                    else:
                        print("Wrong Neumann condition configuration")

                    geomcoor[ipath,:] = biRatGrad[0,:]@P[idR,:]
                    fieldpatch[ipath,:] = neumannvec.T
                    ipath += 1
                # End ppath for loop

                if ineumann == 0:
                    loadcoor = geomcoor
                    loadfield = fieldpatch
                else:
                    # print("******")
                    # print(fieldpatch)
                    loadcoor = np.vstack((loadcoor,geomcoor))
                    loadfield = np.vstack((loadfield,fieldpatch))
                # End if
            # End if
        # End ineumann for loop
    # End if boundary is not None

    if len(loadcoor) != 0:
        # neumannplot = ax.scatter(loadcoor[:,0],loadcoor[:,1],c = "b",marker = "s")
        neumannplot = ax.quiver(loadcoor[:,0],loadcoor[:,1],loadfield[:,0],loadfield[:,1],color=['b'])

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(titlestring)
        # Uncomment for example2
        plt.legend((dirichletplot,neumannplot),(first_string,second_string),loc='upper right',bbox_to_anchor=(1.2,1.0))
        # Uncomment for example3
        # plt.legend((dirichletplot,neumannplot),('Displacement restrictions','Load conditions'),loc='lower right',bbox_to_anchor=(1.2,0.0))
    else:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(titlestring)
        # Uncomment for example2
        plt.legend([dirichletplot],[first_string],loc='upper right',bbox_to_anchor=(1.2,1.0))
        # Uncomment for example3
        # plt.legend((dirichletplot),('Displacement restrictions'),loc='lower right',bbox_to_anchor=(1.2,0.0))
    # End if

    plt.tight_layout()
    plt.show()
# End function
