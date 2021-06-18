# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs
import src.plottingScripts as plts
import src.preprocessor2D as pre2D

################ PREPROCESSING ####################

def multiPatchProblemPreprocessing(phenomenon,multisurface,dirichletconditionsdata,neumannconditionsdata):
    surfacePreprocessing = []
    boundaryPreprocessing = []
    dirichletBCList = []
    enforcedDOF = []
    enforcedValues = []

    multiU,multiV,multip,multiq,multiP,multiw,globalPatchIndices = multisurface.retrieveSurfaceInformation()
    numpatches = len(multiU)
    numdirichletzones = 0
    loadedpatches = []

    for ipatch in range(numpatches):
        U_i = multiU[ipatch]
        V_i = multiV[ipatch]
        p_i = multip[ipatch]
        q_i = multiq[ipatch]
        P_i = multiP[ipatch]
        w_i = multiw[ipatch]
        parametricNodes_i,nodesInElement_i,surfacePreprocessing_i = pre2D.parametricGrid(U_i,V_i,p_i,q_i)
        if dirichletconditionsdata[ipatch] is not None:
            print("There are Dirichlet BC on the patch #{0}".format(ipatch))
            numdirichletzones += 1
            if phenomenon == "Elasticity":
                surface_i = rbs.NURBSSurface(P_i,w_i,p_i,q_i,U=U_i,V=V_i)
                dirichletBCList_i,enforcedDOF_i,enforcedValues_i = \
                pre2D.dirichletBCPreprocessing_Elasticity(P_i,surface_i,[dirichletconditionsdata[ipatch]],U_i,V_i,p_i,q_i)
            # End if
            if phenomenon == "Heat":
                dirichletBCList_i,enforcedDOF_i,enforcedValues_i = \
                pre2D.dirichletBCPreprocessing([dirichletconditionsdata[ipatch]],P_i,U_i,V_i,p_i,q_i)
            # End if
        else:
            dirichletBCList_i = None
            enforcedDOF_i = None
            enforcedValues_i = None
        # End if

        if neumannconditionsdata[ipatch] is not None:
            print("There are Neumann BC on the patch #{0}".format(ipatch))
            boundaryPreprocessing_i = pre2D.neumannBCPreprocessing(parametricNodes_i,nodesInElement_i\
                                      ,neumannconditionsdata[ipatch],U_i,V_i,p_i,q_i)
            loadedpatches.append(ipatch)
            boundaryPreprocessing_i.insert(0,ipatch)
        else:
            boundaryPreprocessing_i = None
        # End if

        surfacePreprocessing.append(surfacePreprocessing_i)
        boundaryPreprocessing.append(boundaryPreprocessing_i)
        if dirichletBCList_i is not None:
            dirichletBCList_patch = [[ipatch] + dbci for dbci in dirichletBCList_i]
            dirichletBCList += dirichletBCList_patch
            enforcedDOF += enforcedDOF_i
            enforcedValues += enforcedValues_i
        # End if
    # End ipatch for loop
    # print(enforcedDOF)
    # print(enforcedValues)
    return surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues
# End function

def plotMultiPatchGeometry(phenomenon,multisurface,dirichletconds,boundaryprep):
    if phenomenon == "Elasticity":
        first_string = "Displacement BC"
        second_string = "Load BC"
    elif phenomenon == "Heat":
        first_string = "Temperature BC"
        second_string = "Flux BC"
    else:
        first_string = "Null BC"
        second_string = "Null BC"
    # End if

    fig = plt.figure()
    ax = plt.axes()
    plt.axis('equal')
    ax.use_sticky_edges = False
    titlestring = "Geometry with boundary conditions"
    ax.axis("off")

    multiU,multiV,multip,multiq,multiP,multiw,globalPatchIndices = multisurface.retrieveSurfaceInformation()

    # Boundary Geometry
    numpatches = len(multiU)

    for ipatch in range(0,numpatches):
        Ui = multiU[ipatch]
        Vi = multiV[ipatch]

        pi = multip[ipatch]
        qi = multiq[ipatch]

        Pi = multiP[ipatch]
        wi = multiw[ipatch]

        surface_i = rbs.NURBSSurface(Pi,wi,pi,qi,U=Ui,V=Vi)
        cbpts = surface_i.createBoundary()
        # cbpts = rbs.nurbs2DBoundary(Ui,Vi,pi,qi,Pi,wi)
        fieldplot = ax.fill(cbpts[:,0],cbpts[:,1],facecolor='none',edgecolor='black',linewidth=1.5)
    # End for loop

    # Control Points
    # ax.scatter(multisurface.fullP[:,0],multisurface.fullP[:,1])

    # Dirichlet Conditions
    dirctrlpts = []
    for drchcond in dirichletconds:
        # print(drchcond)
        idpatch = drchcond[0]
        idnode = drchcond[1]
        dirctrlpts.append(multiP[idpatch][idnode,:])
    # End for loop
    dirctrlpts = np.array(dirctrlpts)

    # Improve this part
    for i in range(0,len(dirichletconds)):
        dofs = dirichletconds[i][3]
        if len(dofs) == 2:
            dirichletplot = ax.scatter(dirctrlpts[:,0],dirctrlpts[:,1],c = "r",marker = "^")
        else:
            dirichletplot = ax.scatter(dirctrlpts[:,0],dirctrlpts[:,1],c = "g",marker = "o")
        # End if
    # End for loop

    jmat = np.zeros((2,2))
    numpt = 5
    loadcoor = []
    loadfield = []

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    loadedpatches = []
    # print(boundaryprep)

    for ipatch in range(0,numpatches):
        if boundaryprep[ipatch] is not None:
            # print(boundaryprep[ipatch])
            loadedpatches = boundaryprep[ipatch][0]
            nonzeroctrlptsload = boundaryprep[ipatch][1]
            boundaryspan = boundaryprep[ipatch][2]
            boundarycorners = boundaryprep[ipatch][3]
            axisselector = boundaryprep[ipatch][4]
            valuesload = boundaryprep[ipatch][6]
            loadtype = boundaryprep[ipatch][5]

            # Neumann Conditions
            for iload in range(0,len(boundarycorners)):
                # Extracting the indices of the patch
                # ipatch = loadedpatches[iload]
                ipatch = loadedpatches
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
                # Extracting the value of the load in the face
                load = valuesload[iload]

                Ui = multiU[ipatch]
                Vi = multiV[ipatch]

                pi = multip[ipatch]
                qi = multiq[ipatch]

                Pi = multiP[ipatch]
                wi = multiw[ipatch]

                mu = len(Ui) - 1
                mv = len(Vi) - 1

                Pwl = rbs.weightedControlPoints(Pi,wi)
                Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

                parampath = np.linspace(startpt,endpt,numpt,endpoint=True)
                geomcoor = np.zeros((parampath.shape[0],2))
                fieldpatch = np.zeros((parampath.shape[0],2))
                ipath = 0
                for ppath in parampath:
                    biRatGrad = rbs.bivariateRationalGradient(mu,mv,pi,qi,uspan,vspan,ppath[0],ppath[1],Ui,Vi,Pwi)

                    jmat = (biRatGrad[1:3,:]@Pi[idR,:]).T

                    jvec = jmat[:,paramaxis]
                    normjvec = np.linalg.norm(jvec)

                    if normjvec > 1e-6:
                        unitTangetVec = jvec/normjvec
                    else:
                        unitTangetVec = np.zeros((2,1))
                    # End if

                    if loadtype[iload] == "tangent":
                        loadvec = (load/abs(load))*unitTangetVec
                    elif loadtype[iload] == "normal":
                        unitNormalVec = rotMat@unitTangetVec
                        loadvec = (load/abs(load))*unitNormalVec
                    else:
                        print("Wrong load configuration")
                    # End if

                    geomcoor[ipath,:] = biRatGrad[0,:]@Pi[idR,:]
                    fieldpatch[ipath,:] = loadvec.T
                    ipath += 1
                # End for loop
                loadcoor.append(geomcoor)
                loadfield.append(fieldpatch)
            # End for loop
        # End if
    # End for loop

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
    # plt.legend((dirichletplot,neumannplot),('Displacement restrictions','Load conditions'),loc='upper right',bbox_to_anchor=(1.2,1.0))
    # Uncomment for example3
    plt.legend((dirichletplot,neumannplot),('Restrictions','Loads'),loc='right',bbox_to_anchor=(1.2,0.5))
    plt.tight_layout()
    plt.show()
# End function
