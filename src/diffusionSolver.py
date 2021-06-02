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

def conductivityMatrix(kappa):
    kmat = np.zeros((2,2))
    kmat[0][0] = kappa
    kmat[1][1] = kappa
    return kmat

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(surface,surfaceprep,numquad,matprop,boundaryprep):
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()

    K = np.zeros((P.shape[0],P.shape[0]))
    F = np.zeros((P.shape[0],1))
    Fb = np.zeros((P.shape[0],1))
    Fl = np.zeros((P.shape[0],1))
    M = np.zeros((P.shape[0],P.shape[0]))

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
    kappa = matprop[0]
    rho = matprop[1]
    source = matprop[2]
    kMat = conductivityMatrix(kappa)

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
            K[globalDOFx,globalDOFy] += (dN2dxi.T@kMat@dN2dxi)*wJac

            # Body forces integral
            if abs(source) > 1e-5:
                nMat = biRatGrad[0,None,:]
                Fb[globalDOF] += nMat.T*source*wJac

            # Mass integral
            if abs(rho) > 1e-5:
                nMat = biRatGrad[0,:]
                nMat = np.reshape(nMat,(1,len(nMat)))
                M[globalDOFx,globalDOFy] += rho*(nMat.T@nMat)*wJac
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

    F = Fb + Fl
    return K,F,M

def assemblyMultipatchWeakForm(multisurface,surfaceprep,numquad,matprop,boundaryprep):
    multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints = multisurface.retrieveSurfaceInformation()

    Ktotal = np.zeros((2*fullP.shape[0],2*fullP.shape[0]))
    Ftotal = np.zeros((2*fullP.shape[0],1))
    Fbtotal = np.zeros((2*fullP.shape[0],1))
    Fltotal = np.zeros((2*fullP.shape[0],1))

    numquad2d = numquad[0]
    numquad1d = numquad[1]

    # Definition of the material matrix
    E = matprop[0]
    nu = matprop[1]
    rho = matprop[2]
    dMat = elasticMatrix(E,nu)

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    paramGrad = np.zeros((2,2))

    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8

    numpatches = len(multiU)

    # Patch loop
    print('Computing the strain-energy and body forces integrals')
    for ipatch in range(0,numpatches):
        Ui = multiU[ipatch]
        Vi = multiV[ipatch]

        pi = multip[ipatch]
        qi = multiq[ipatch]

        Pi = fullP[idcontrolpoints[ipatch],:]
        wi = fullw[idcontrolpoints[ipatch],:]

        Kpatch = np.zeros((2*Pi.shape[0],2*Pi.shape[0]))
        Fbpatch = np.zeros((2*Pi.shape[0],1))

        # Extraction of surface preprocessing
        nonzeroctrlpts = surfaceprep[ipatch][0]
        surfacespan = surfaceprep[ipatch][1]
        elementcorners = surfaceprep[ipatch][2]
        numElems = len(elementcorners)

        Pwl = rbs.weightedControlPoints(Pi,wi)
        Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

        # Precomputing info for the nurbs derivatives
        mU = len(Ui) - 1
        mV = len(Vi) - 1
        nU = mU - pi - 1
        nV = mV - qi - 1

        # Global degrees of freedom
        globalDOF = np.zeros(2*len(idcontrolpoints[ipatch]),dtype=int)
        dof0 = 2*np.array(idcontrolpoints[ipatch])
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1

        globalDOFx,globalDOFy = np.meshgrid(globalDOF,globalDOF,indexing='xy')

        # Strain-energy and body force integrals
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

            # Patch degrees of freedom
            patchDOF = np.zeros(2*len(idR),dtype=int)
            dof0 = 2*np.array(idR)
            dof1 = dof0 + 1
            patchDOF[0::2] = dof0
            patchDOF[1::2] = dof1

            patchDOFx,patchDOFy = np.meshgrid(patchDOF,patchDOF,indexing='xy')

            # K stiffness matrix
            for iquad in range(numquad2d.shape[0]):
                coor = parametricCoordinate(aPoint[0],cPoint[0],aPoint[1],cPoint[1],numquad2d[iquad][0],numquad2d[iquad][1])

                # NURBS gradient
                biRatGrad = rbs.bivariateRationalGradient(mU,mV,pi,qi,uspan,vspan,coor[0][0],coor[0][1],Ui,Vi,Pwi)

                # Jacobian
                jac = (biRatGrad[1:3,:]@Pi[idR,:]).T
                wJac = abs(np.linalg.det(jac))*detJac2*numquad2d[iquad][2]

                # Strain displacement matrix
                invJac = np.linalg.inv(jac)
                dN2 = biRatGrad[1:3,:]
                dN2dxi = invJac.T@dN2

                bMat = np.zeros((3,2*dN2dxi.shape[1]))
                #dNx
                bMat[0,0::2] = dN2dxi[0,:]
                bMat[2,0::2] = dN2dxi[1,:]
                #dNy
                bMat[1,1::2] = dN2dxi[1,:]
                bMat[2,1::2] = dN2dxi[0,:]

                # Patch indexing
                Kpatch[patchDOFx,patchDOFy] += (bMat.T@dMat@bMat)*wJac

                # Body forces integral
                if abs(rho) > 1e-5:
                    nMat = np.zeros((2,2*biRatGrad.shape[1]))
                    nMat[0,0::2] = biRatGrad[0,:]
                    nMat[1,1::2] = biRatGrad[0,:]

                    Fbpatch[patchDOF] += (nMat.T@bvec)*wJac
            # End of quadrature loop
        # End of element loop
        Ktotal[globalDOFx,globalDOFy] += Kpatch
        Fbtotal[globalDOF] += Fbpatch
    # End of patch loop

    # Load conditions loop
    # Extraction of boundary preprocessing
    loadedpatches = boundaryprep[0]
    nonzeroctrlptsload = boundaryprep[1]
    boundaryspan = boundaryprep[2]
    boundarycorners = boundaryprep[3]
    axisselector = boundaryprep[4]
    valuesload = boundaryprep[5]
    neumanntype = boundaryprep[6]

    numLoadedElems = len(boundarycorners)

    # Load integrals
    print('Computing the load integrals')
    for ineumann in range(0,numLoadedElems):
        # Extracting the indices of the patch
        ipatch = loadedpatches[ineumann]
        # Extracting the indices of the non-zero control points of the loaded elements
        idR = nonzeroctrlptsload[ineumann]
        # Extracting the indices for the location of the parametric segment
        uspan = boundaryspan[ineumann][0]
        vspan = boundaryspan[ineumann][1]
        # Extracting the corners of the parametric segment
        aPoint = boundarycorners[ineumann][0]
        bPoint = boundarycorners[ineumann][1]
        # Extracting the non-zero column index of the boundary jacobian
        paramaxis = axisselector[ineumann]
        # Extracting the value of the load in the face
        load = valuesload[ineumann]

        Ui = multiU[ipatch]
        Vi = multiV[ipatch]

        pi = multip[ipatch]
        qi = multiq[ipatch]

        Pi = fullP[idcontrolpoints[ipatch],:]
        wi = fullw[idcontrolpoints[ipatch],:]

        Flpatch = np.zeros((2*Pi.shape[0],1))

        mu = len(Ui) - 1
        mv = len(Vi) - 1

        Pwl = rbs.weightedControlPoints(Pi,wi)
        Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

        # Global degrees of freedom
        globalDOF = np.zeros(2*len(idcontrolpoints[ipatch]),dtype=int)
        dof0 = 2*np.array(idcontrolpoints[ipatch])
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1

        globalDOFx,globalDOFy = np.meshgrid(globalDOF,globalDOF,indexing='xy')

        # Patch degrees of freedom
        patchDOF = np.zeros(2*len(idR),dtype=int)
        dof0 = 2*np.array(idR)
        dof1 = dof0 + 1
        patchDOF[0::2] = dof0
        patchDOF[1::2] = dof1

        patchDOFx,patchDOFy = np.meshgrid(patchDOF,patchDOF,indexing='xy')

        for iquad in range(numquad1d.shape[0]):
            coor = parametricCoordinate(aPoint[0],bPoint[0],aPoint[1],bPoint[1],numquad1d[iquad][0],numquad1d[iquad][0])

            biRatGrad = rbs.bivariateRationalGradient(mU,mV,pi,qi,uspan,vspan,coor[0][0],coor[0][1],Ui,Vi,Pwi)
            Jac = (biRatGrad[1:3,:]@Pi[idR,:]).T
            jac1 = np.linalg.norm(Jac[:,paramaxis])
            jac2 = 0.5*np.sum(bPoint-aPoint)

            if jac1 > 1e-6:
                unitTangentVec = Jac[:,paramaxis]/jac1
            else:
                unitTangentVec = np.zeros((2,1))

            if neumanntype[ineumann] == "tangent":
                tvec = load*unitTangentVec
            elif neumanntype[ineumann] == "normal":
                unitNormalVec = rotMat@unitTangentVec
                tvec = load*unitNormalVec
            else:
                print("Wrong load configuration")

            tvec = np.reshape(tvec,(2,1))
            nMat = np.zeros((2,2*biRatGrad.shape[1]))

            nMat[0,0::2] = biRatGrad[0,:]
            nMat[1,1::2] = biRatGrad[0,:]

            Flpatch[patchDOF] += (nMat.T@tvec)*jac1*jac2*numquad1d[iquad][1]
        # End quadrature loop
        Fltotal[globalDOF] += Flpatch
    # End load loop

    Ftotal = Fbtotal + Fltotal
    return Ktotal,Ftotal
