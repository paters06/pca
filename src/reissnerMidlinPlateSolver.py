# Python libraries
import numpy as np
import numpy.linalg

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs

def parametricCoordinate(ua,ub,va,vb,gausspta,gaussptb):
    localpts = np.zeros((1,2))
    localpts[0][0] = 0.5*(ub - ua)*gausspta + 0.5*(ub + ua)
    localpts[0][1] = 0.5*(vb - va)*gaussptb + 0.5*(vb + va)
    return localpts
# End function

def elasticMatrix_RM(E,nu,t):
    db = np.zeros((3,3))
    ds = np.zeros((2,2))

    # Db
    db[0][0] = 1.0/(1.0 - nu)
    db[1][1] = 1.0/(1.0 - nu)
    db[2][2] = 2.0
    db[0][1] = nu/(1.0 - nu)
    db[1][0] = nu/(1.0 - nu)
    db *= (E*t**3)/(12*(1+nu))

    # Ds
    ds[0][0] = 1.0
    ds[1][1] = 1.0
    ds *= (0.5*t*E)/(1+nu)

    return db,ds
# End function

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(surface,surfaceprep,numquad,matprop,boundaryprep):
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()

    K = np.zeros((3*P.shape[0],3*P.shape[0]))
    KBending = np.zeros((3*P.shape[0],3*P.shape[0]))
    KShear = np.zeros((3*P.shape[0],3*P.shape[0]))
    F = np.zeros((3*P.shape[0],1))
    Fb = np.zeros((3*P.shape[0],1))
    Fl = np.zeros((3*P.shape[0],1))
    M = np.zeros((3*P.shape[0],3*P.shape[0]))

    numquad2d = numquad[0]
    numquad1d = numquad[1]

    # Extraction of surface preprocessing
    nonzeroctrlpts = surfaceprep[0]
    surfacespan = surfaceprep[1]
    elementcorners = surfaceprep[2]

    paramGrad = np.zeros((2,2))
    numElems = len(elementcorners)

    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    # Definition of the material matrix
    E = matprop[0]
    nu = matprop[1]
    rho = matprop[2]
    load = matprop[3]
    t = matprop[4]
    dBMat,dSMat = elasticMatrix_RM(E,nu,t)

    bvec = np.zeros((3,1))
    bvec[0][0] = -rho*9.8

    loadvec = np.zeros((3,1))
    loadvec[0][0] = load

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

        # Global degrees of freedom
        globalDOF = np.zeros(3*len(idR),dtype=int)
        dof0 = 3*np.array(idR)
        dof1 = dof0 + 1
        dof2 = dof0 + 2
        globalDOF[0::3] = dof0
        globalDOF[1::3] = dof1
        globalDOF[2::3] = dof2

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
            N2 = biRatGrad[0,:]
            dN2 = biRatGrad[1:3,:]
            dN2dxi = invJac.T@dN2

            bBMat = np.zeros((3,3*dN2dxi.shape[1]))
            bSMat = np.zeros((2,3*dN2dxi.shape[1]))

            # First row - Bb
            bBMat[0,1::3] = dN2dxi[0,:]

            # Second row - Bb
            bBMat[1,2::3] = dN2dxi[1,:]

            # Third row - Bb
            bBMat[2,1::3] = dN2dxi[1,:]
            bBMat[2,2::3] = dN2dxi[0,:]

            # First row - Bs
            bSMat[0,0::3] = dN2dxi[0,:]
            bSMat[0,1::3] = -N2

            # Second row - Bs
            bSMat[1,0::3] = dN2dxi[1,:]
            bSMat[1,2::3] = -N2

            # Global indexing
            KBending[globalDOFx,globalDOFy] += (bBMat.T@dBMat@bBMat)*wJac
            KShear[globalDOFx,globalDOFy] += (bSMat.T@dSMat@bSMat)*wJac

            # Body forces and mass integrals
            if abs(rho) > 1e-5 or abs(load) > 1e-5:
                nMat = np.zeros((3,3*biRatGrad.shape[1]))
                nMat[0,0::3] = biRatGrad[0,:]
                nMat[1,1::3] = biRatGrad[0,:]
                nMat[2,2::3] = biRatGrad[0,:]

                Fb[globalDOF] += (nMat.T@bvec)*wJac
                Fl[globalDOF] += (nMat.T@loadvec)*wJac
                M[globalDOFx,globalDOFy] += rho*(nMat.T@nMat)*wJac
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
        # Load integrals
        print('Computing the load integrals')
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
            globalDOF = np.zeros(3*len(idR),dtype=int)
            dof0 = 3*np.array(idR)
            dof1 = dof0 + 1
            dof2 = dof0 + 2
            globalDOF[0::3] = dof0
            globalDOF[1::3] = dof1
            globalDOF[2::3] = dof2

            for iquad in range(numquad1d.shape[0]):
                coor = parametricCoordinate(aPoint[0],bPoint[0],aPoint[1],bPoint[1],numquad1d[iquad][0],numquad1d[iquad][0])

                biRatGrad = rbs.bivariateRationalGradient(mU,mV,p,q,uspan,vspan,coor[0][0],coor[0][1],U,V,Pw)
                Jac = (biRatGrad[1:3,:]@P[idR,:]).T
                jac1 = np.linalg.norm(Jac[:,paramaxis])
                jac2 = 0.5*np.sum(bPoint-aPoint)

                if jac1 > 1e-6:
                    unitTangetVec = Jac[:,paramaxis]/jac1
                else:
                    unitTangetVec = np.zeros((2,1))

                if neumanntype == "tangent":
                    tvec = neumannval*unitTangetVec
                elif neumanntype == "normal":
                    unitNormalVec = rotMat@unitTangetVec
                    tvec = neumannval*unitNormalVec
                else:
                    print("Wrong load configuration")

                tvec = np.reshape(tvec,(3,1))
                nMat = np.zeros((3,3*biRatGrad.shape[1]))

                nMat[0,0::3] = biRatGrad[0,:]
                nMat[1,1::3] = biRatGrad[0,:]
                nMat[2,2::3] = biRatGrad[0,:]

                Fl[globalDOF] += (nMat.T@tvec)*jac1*jac2*numquad1d[iquad][1]
            # End iquad loop
        # End iload loop
    # End if boundary is not None

    K = KBending + KShear
    F = Fb + Fl
    return K,F,M
# End function