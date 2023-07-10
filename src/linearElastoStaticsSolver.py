# Python libraries
import numpy as np
# import numpy.linalg

# Local project
# import src.basisFunctions as bfunc
import src.nurbs as rbs

def parametricCoordinate(ua,ub,va,vb,gausspta,gaussptb):
    localpts = np.zeros((1,2))
    localpts[0][0] = 0.5*(ub - ua)*gausspta + 0.5*(ub + ua)
    localpts[0][1] = 0.5*(vb - va)*gaussptb + 0.5*(vb + va)
    return localpts

def elasticMatrix(E,nu):
    dmat = np.zeros((3,3))
    dmat[0][0] = 1 - nu
    dmat[1][1] = 1 - nu
    dmat[2][2] = (1 - 2*nu)/2
    dmat[0][1] = nu
    dmat[1][0] = nu
    dmat *= E/((1+nu)*(1-2*nu))
    return dmat

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(surface,surfaceprep,numquad,matprop,boundaryprep):
    U,V,p,q,P,w = surface.retrieveSurfaceInformation()

    K = np.zeros((2*P.shape[0],2*P.shape[0]))
    F = np.zeros((2*P.shape[0],1))
    Fb = np.zeros((2*P.shape[0],1))
    Fl = np.zeros((2*P.shape[0],1))
    M = np.zeros((2*P.shape[0],2*P.shape[0]))

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
    dMat = elasticMatrix(E,nu)

    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8

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
        globalDOF = np.zeros(2*len(idR),dtype=int)
        dof0 = 2*np.array(idR)
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1

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

            bMat = np.zeros((3,2*dN2dxi.shape[1]))
            #dNx
            bMat[0,0::2] = dN2dxi[0,:]
            bMat[2,0::2] = dN2dxi[1,:]
            #dNy
            bMat[1,1::2] = dN2dxi[1,:]
            bMat[2,1::2] = dN2dxi[0,:]

            # Global indexing
            K[globalDOFx,globalDOFy] += (bMat.T@dMat@bMat)*wJac

            # Body forces and mass integrals
            if abs(rho) > 1e-5:
                nMat = np.zeros((2,2*biRatGrad.shape[1]))
                nMat[0,0::2] = biRatGrad[0,:]
                nMat[1,1::2] = biRatGrad[0,:]

                Fb[globalDOF] += (nMat.T@bvec)*wJac
                M[globalDOFx,globalDOFy] += rho*(nMat.T@nMat)*wJac

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
            globalDOF = np.zeros(2*len(idR),dtype=int)
            dof0 = 2*np.array(idR)
            dof1 = dof0 + 1
            globalDOF[0::2] = dof0
            globalDOF[1::2] = dof1

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

                tvec = np.reshape(tvec,(2,1))
                nMat = np.zeros((2,2*biRatGrad.shape[1]))

                nMat[0,0::2] = biRatGrad[0,:]
                nMat[1,1::2] = biRatGrad[0,:]

                Fl[globalDOF] += (nMat.T@tvec)*jac1*jac2*numquad1d[iquad][1]

    F = Fb + Fl
    return K,F,M

def assemblyMultipatchWeakForm(multisurface,surfaceprep,numquad,matprop,boundaryprep):
    multiU,multiV,multip,multiq,multiP,multiw,globalPatchIndices = multisurface.retrieveSurfaceInformation()

    multisurface.createFullControlPolygon()
    fullNumberControlPoints = multisurface.fullP.shape[0]
    Ktotal = np.zeros((2*fullNumberControlPoints,2*fullNumberControlPoints))
    Mtotal = np.zeros((2*fullNumberControlPoints,2*fullNumberControlPoints))
    Ftotal = np.zeros((2*fullNumberControlPoints,1))
    Fbtotal = np.zeros((2*fullNumberControlPoints,1))
    Fltotal = np.zeros((2*fullNumberControlPoints,1))

    # numquad2d = numquad[0]
    # numquad1d = numquad[1]

    # Definition of the material matrix
    # E = matprop[0]
    # nu = matprop[1]
    # rho = matprop[2]
    # dMat = elasticMatrix(E,nu)

    # Rotation matrix for -pi/2
    # rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    # paramGrad = np.zeros((2,2))

    # bvec = np.zeros((2,1))
    # bvec[1][0] = -rho*9.8

    numpatches = len(multiU)

    # Patch loop
    print('Computing the strain-energy and body forces integrals')
    for ipatch in range(numpatches):
        print('Patch #',ipatch)

        Uinit = multiU[ipatch]
        Vinit = multiV[ipatch]
        pinit = multip[ipatch]
        qinit = multiq[ipatch]
        Pinit = multiP[ipatch]
        winit = multiw[ipatch]

        surface_i = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)
        surfaceprep_i = surfaceprep[ipatch]

        print(Pinit)

        if boundaryprep[ipatch] is not None:
            boundaryprep_i = boundaryprep[ipatch][1:]
            # print(boundaryprep_i)
        else:
            boundaryprep_i = boundaryprep[ipatch]

        Kpatch,Fpatch,Mpatch = assemblyWeakForm(surface_i,surfaceprep_i,numquad,matprop,boundaryprep_i)

        # Patch degrees of freedom
        # print(globalPatchIndices[ipatch])
        patchDOF = np.zeros(2*len(globalPatchIndices[ipatch]),dtype=int)
        dof0 = 2*np.array(globalPatchIndices[ipatch])
        dof1 = dof0 + 1.0
        patchDOF[0::2] = dof0
        patchDOF[1::2] = dof1
        # print(patchDOF)

        patchDOFx,patchDOFy = np.meshgrid(patchDOF,patchDOF,indexing='xy')

        Ktotal[patchDOFx,patchDOFy] += Kpatch
        Mtotal[patchDOFx,patchDOFy] += Mpatch
        Ftotal[patchDOF] += Fpatch

    return Ktotal,Ftotal,Mtotal

def assemblyMultipatchWeakFormv2(multisurface,surfaceprep,numquad,matprop,boundaryprep):
    multiU,multiV,multip,multiq,multiP,multiw,globalPatchIndices = multisurface.retrieveSurfaceInformation()

    #fullP is not known
    multisurface.createFullControlPolygon()
    fullNumberControlPoints = multisurface.fullP.shape[0]
    Ktotal = np.zeros((2*fullNumberControlPoints,2*fullNumberControlPoints))
    Mtotal = np.zeros((2*fullNumberControlPoints,2*fullNumberControlPoints))
    Ftotal = np.zeros((2*fullNumberControlPoints,1))
    Fbtotal = np.zeros((2*fullNumberControlPoints,1))
    Fltotal = np.zeros((2*fullNumberControlPoints,1))

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

    print(len(surfaceprep))

    # Patch loop
    print('Computing the strain-energy and body forces integrals')
    for ipatch in range(0,numpatches):
        Ui = multiU[ipatch]
        Vi = multiV[ipatch]

        pi = multip[ipatch]
        qi = multiq[ipatch]

        Pi = multiP[ipatch]
        wi = multiw[ipatch]

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
        globalDOF = np.zeros(2*len(globalPatchIndices[ipatch]),dtype=int)
        dof0 = 2*np.array(globalPatchIndices[ipatch])
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
        Ktotal[globalDOFx,globalDOFy] += Kpatch
        Fbtotal[globalDOF] += Fbpatch

        if boundaryprep[ipatch] is not None:
            # Load conditions loop
            # Extraction of boundary preprocessing
            loadedpatches = boundaryprep[ipatch][0]
            nonzeroctrlptsload = boundaryprep[ipatch][1]
            boundaryspan = boundaryprep[ipatch][2]
            boundarycorners = boundaryprep[ipatch][3]
            axisselector = boundaryprep[ipatch][4]
            valuesload = boundaryprep[ipatch][6]
            loadtype = boundaryprep[ipatch][5]

            numLoadedElems = len(boundarycorners)

            # Load integrals
            print('Computing the load integrals')
            for iload in range(0,numLoadedElems):
                # Extracting the indices of the patch
                # ipatch = loadedpatches[iload]
                ipatch = loadedpatches
                # Extracting the indices of the non-zero control points of the loaded elements
                idR = nonzeroctrlptsload[iload]
                # Extracting the indices for the location of the parametric segment
                uspan = boundaryspan[iload][0]
                vspan = boundaryspan[iload][1]
                # Extracting the corners of the parametric segment
                aPoint = boundarycorners[iload][0]
                bPoint = boundarycorners[iload][1]
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

                Flpatch = np.zeros((2*Pi.shape[0],1))

                mu = len(Ui) - 1
                mv = len(Vi) - 1

                Pwl = rbs.weightedControlPoints(Pi,wi)
                Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

                # Global degrees of freedom
                globalDOF = np.zeros(2*len(globalPatchIndices[ipatch]),dtype=int)
                dof0 = 2*np.array(globalPatchIndices[ipatch])
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
                        unitTangetVec = Jac[:,paramaxis]/jac1
                    else:
                        unitTangetVec = np.zeros((2,1))

                    if loadtype[iload] == "tangent":
                        # tvec = (loadvalue/abs(loadvalue))*unitTangetVec
                        tvec = load*unitTangetVec
                    elif loadtype[iload] == "normal":
                        unitNormalVec = rotMat@unitTangetVec
                        # tvec = (loadvalue/abs(loadvalue))*unitNormalVec
                        tvec = load*unitNormalVec
                    else:
                        print("Wrong load configuration")

                    tvec = np.reshape(tvec,(2,1))
                    nMat = np.zeros((2,2*biRatGrad.shape[1]))

                    nMat[0,0::2] = biRatGrad[0,:]
                    nMat[1,1::2] = biRatGrad[0,:]

                    Flpatch[patchDOF] += (nMat.T@tvec)*jac1*jac2*numquad1d[iquad][1]
                Fltotal[globalDOF] += Flpatch

    Ftotal = Fbtotal + Fltotal
    return Ktotal,Ftotal, Mtotal
