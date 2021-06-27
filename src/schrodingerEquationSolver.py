# Python libraries
import numpy as np
# import numpy.linalg

# Local project
import src.nurbs as rbs

def parametricCoordinate(ua,ub,va,vb,gausspta,gaussptb):
    localpts = np.zeros((1,2))
    localpts[0][0] = 0.5*(ub - ua)*gausspta + 0.5*(ub + ua)
    localpts[0][1] = 0.5*(vb - va)*gaussptb + 0.5*(vb + va)
    return localpts

def conductivityMatrix(kappa):
    kmat = np.zeros((2,2))
    kmat[0][0] = kappa
    kmat[1][1] = kappa
    return kmat

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(surface,surfaceprep,numquad,potentialfunction,boundaryprep):
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()

    KMat = np.zeros((P.shape[0],P.shape[0]))
    TMat = np.zeros((P.shape[0],P.shape[0]))
    VMat = np.zeros((P.shape[0],P.shape[0]))
    FVec = np.zeros((P.shape[0],1))
    Fb = np.zeros((P.shape[0],1))
    Fl = np.zeros((P.shape[0],1))
    MMat = np.zeros((P.shape[0],P.shape[0]))

    print("Size of K matrix: {}".format(KMat.shape))
    print("Size of M matrix: {}".format(MMat.shape))
    print("Size of F vector: {}".format(FVec.shape))

    numquad2d = numquad[0]
    numquad1d = numquad[1]

    # Extraction of surface preprocessing
    nonzeroctrlpts_surface = surfaceprep[0]
    surfacespan = surfaceprep[1]
    elementcorners = surfaceprep[2]

    paramGrad = np.zeros((2,2))
    numElems = len(elementcorners)

    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    # Definition of the material matrix
    # kappa = matprop[0]
    # rho = matprop[1]
    # source = matprop[2]
    # kMat = conductivityMatrix(kappa)
    rho = 1.0

    # Precomputing info for the nurbs derivatives
    mU = len(U) - 1
    mV = len(V) - 1
    nU = mU - p - 1
    nV = mV - q - 1

    # Strain-energy and body force integrals
    print('Computing the strain-energy and body forces integrals')
    for ielem in range(0,numElems):
        # Extracting the indices of the non-zero control points
        idR = nonzeroctrlpts_surface[ielem]
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

        # Global degrees of freedom
        globalDOF = idR
        globalDOFx,globalDOFy = np.meshgrid(globalDOF,globalDOF,indexing='xy')

        # K stiffness matrix
        for iquad in range(numquad2d.shape[0]):
            coor = parametricCoordinate(aPoint[0],cPoint[0],aPoint[1],cPoint[1],numquad2d[iquad][0],numquad2d[iquad][1])

            # NURBS gradient
            biRatGrad = rbs.bivariateRationalGradient(mU,mV,p,q,uspan,vspan,coor[0][0],coor[0][1],U,V,Pw)

            # Jacobian
            jac = (biRatGrad[1:3,:]@P[idR,:]).T
            wJac = abs(np.linalg.det(jac))*detJac2*numquad2d[iquad][2]

            # Strain displacement matrix
            invJac = np.linalg.inv(jac)
            dN2 = biRatGrad[1:3,:]
            dN2dxi = invJac.T@dN2

            # Global indexing
            TMat[globalDOFx,globalDOFy] += (dN2dxi.T@dN2dxi)*wJac

            # Body forces integral
            # if abs(source) > 1e-5:
                # nMat = biRatGrad[0,None,:]
                # Fb[globalDOF] += nMat.T*potentialfunction()*wJac
            # End if

            # Mass integral
            if abs(rho) > 1e-5:
                nMat = biRatGrad[0,:]
                nMat = np.reshape(nMat,(1,len(nMat)))
                geomcoor = (nMat@P[idR,:])[0]
                # print(np.sum(nMat))
                # print(geomcoor)
                # print(potentialfunction(geomcoor))
                VMat[globalDOFx,globalDOFy] += potentialfunction(geomcoor)*(nMat.T@nMat)*wJac
                MMat[globalDOFx,globalDOFy] += (nMat.T@nMat)*wJac
            # End if
        # End iquad loop
    # End ielem loop

    if boundaryprep is not None:
        # Extraction of boundary preprocessing
        nonzeroctrlpts_neumannbc = boundaryprep[0]
        boundaryspan = boundaryprep[1]
        boundarycorners = boundaryprep[2]
        axisselector = boundaryprep[3]
        neumannbc_type = boundaryprep[4]
        neumannbc_value = boundaryprep[5]

        numLoadedElems = len(boundarycorners)

        # Neumann boundary integrals
        print('Computing the Neumann boundary integrals')
        for ineumann in range(0,numLoadedElems):
            # Extracting the indices of the non-zero control points of the loaded elements
            idR = nonzeroctrlpts_neumannbc[ineumann]
            # Extracting the indices for the location of the parametric segment
            uspan = boundaryspan[ineumann][0]
            vspan = boundaryspan[ineumann][1]
            # Extracting the corners of the parametric segment
            aPoint = boundarycorners[ineumann][0]
            bPoint = boundarycorners[ineumann][1]
            # Extracting the non-zero column index of the boundary jacobian
            paramaxis = axisselector[ineumann]
            # Extracting the type and value of the neumann conditions
            neumanntype = neumannbc_type[ineumann]
            neumannval = neumannbc_value[ineumann]

            # Global degrees of freedom
            globalDOF = idR

            for iquad in range(numquad1d.shape[0]):
                coor = parametricCoordinate(aPoint[0],bPoint[0],aPoint[1],bPoint[1],numquad1d[iquad][0],numquad1d[iquad][0])

                biRatGrad = rbs.bivariateRationalGradient(mU,mV,p,q,uspan,vspan,coor[0][0],coor[0][1],U,V,Pw)
                Jac = (biRatGrad[1:3,:]@P[idR,:]).T
                jac1 = np.linalg.norm(Jac[:,paramaxis])
                jac2 = 0.5*np.sum(bPoint-aPoint)

                if jac1 > 1e-6:
                    unitTangentVec = Jac[:,paramaxis]/jac1
                else:
                    unitTangetVec = np.zeros((2,1))

                if neumanntype == "tangent":
                    tvec = neumannval*unitTangentVec
                elif neumanntype == "normal":
                    unitNormalVec = rotMat@unitTangetVec
                    tvec = neumannval*unitNormalVec
                else:
                    print("Wrong load configuration")

                tvec = 0.0
                nMat = biRatGrad[0,:]
                nMat = np.reshape(nMat,(len(nMat),1))

                Fl[globalDOF] += nMat*tvec*jac1*jac2*numquad1d[iquad][1]
            # End iquad loop
        # End ineumann loop
    # End if boundary is not None

    FVec = Fb + Fl
    KMat = TMat + VMat

    print("Memory size of K matrix after assembly: {:.4f} MB".format(KMat.nbytes/(1024**2)))
    print("Memory size of M matrix after assembly: {:.4f} MB".format(MMat.nbytes/(1024**2)))
    print("Memory size of F vector after assembly: {:.4f} MB".format(FVec.nbytes/(1024**2)))
    # print(VMat)
    # print(MMat)
    return KMat,FVec,MMat
# End function