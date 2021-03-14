# Python libraries
import math
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs
from src.linearElastoStaticsSolver import elasticMatrix
import src.plottingScripts as plts

################ POSTPROCESSING ####################

def displacementField(numpoints,U,V,p,q,P,D,w,surfaceprep):
    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1
    
    # Extraction of surface preprocessing
    nonzeroctrlpts = surfaceprep[0]
    surfacespan = surfaceprep[1]
    elementcorners = surfaceprep[2]
    
    numElems = len(elementcorners)
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    
    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)

    # Geometric coordinates
    cpts = np.zeros((2,numelemsu*numpoints,numelemsv*numpoints))

    # Displacements
    upts = np.zeros((2,numelemsu*numpoints,numelemsv*numpoints))

    for ielem in range(0,numElems):
        # Extracting the indices of the non-zero control points
        idR = nonzeroctrlpts[ielem]
        # Extracting the indices for the location of the parametric element
        uspan = surfacespan[ielem][0]
        vspan = surfacespan[ielem][1]
        # Extracting the corners of the parametric element
        apt = elementcorners[ielem][0]
        cpt = elementcorners[ielem][1]

        urank = np.linspace(apt[0],cpt[0],numpoints)
        vrank = np.linspace(apt[1],cpt[1],numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):
                R = rbs.bivariateRationalFunction(mu,mv,p,q,uspan,vspan,urank[i],vrank[j],U,V,Pw)
                
                upts[:,iu*numpoints + i,jv*numpoints + j] = R@D[idR,:]
                cpts[:,iu*numpoints + i,jv*numpoints + j] = R@P[idR,:]

    return upts,cpts

# Improve this function
def stressField(numpoints,U,V,p,q,P,w,dtot,matprop,surfaceprep):
    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1
    
    # Extraction of surface preprocessing
    nonzeroctrlpts = surfaceprep[0]
    surfacespan = surfaceprep[1]
    elementcorners = surfaceprep[2]
    
    numElems = len(elementcorners)
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    
    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)
    
    # Definition of the material matrix
    E = matprop[0]
    nu = matprop[1]
    rho = matprop[2]
    dMat = elasticMatrix(E,nu)
    
    paramgrad = np.zeros((2,2))
    sigma = np.zeros((3,numelemsu*numpoints,numelemsv*numpoints))
    nanArray = []

    for ielem in range(0,numElems):
        # Extracting the indices of the non-zero control points
        idR = nonzeroctrlpts[ielem]
        # Extracting the indices for the location of the parametric element
        uspan = surfacespan[ielem][0]
        vspan = surfacespan[ielem][1]
        # Extracting the corners of the parametric element
        apt = elementcorners[ielem][0]
        cpt = elementcorners[ielem][1]

        urank = np.linspace(apt[0],cpt[0],numpoints)
        vrank = np.linspace(apt[1],cpt[1],numpoints)
        
        # Global degrees of freedom
        globalDOF = np.zeros(2*len(idR),dtype=int)
        dof0 = 2*np.array(idR)
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):
                xpcoor = urank[i]
                ypcoor = vrank[j]
                
                biRatGrad = rbs.bivariateRationalGradient(mu,mv,p,q,uspan,vspan,xpcoor,ypcoor,U,V,Pw)

                jac = (biRatGrad[1:3,:]@P[idR,:]).T
                detJac = np.linalg.det(jac)

                if abs(detJac) > 1e-5:
                    invJac = np.linalg.inv(jac)
                    dN2 = biRatGrad[1:3,:]
                    dN2dxi = invJac.T@dN2
                    
                    numpts = dN2dxi.shape[1]
                    bmat = np.zeros((3,2*numpts))
                    #dNx
                    bmat[0,0::2] = dN2dxi[0,:]
                    bmat[2,0::2] = dN2dxi[1,:]
                    #dNy
                    bmat[1,1::2] = dN2dxi[1,:]
                    bmat[2,1::2] = dN2dxi[0,:]
                    
                    svec = dMat@(bmat@dtot[globalDOF,:])
                else:
                    print("Singularity")
                    print([iu*numpoints + i,jv*numpoints + j])
                    svec = np.empty((3,1))
                    svec[:] = np.NaN
                    print(svec.shape)
                    nanArray.append([iu*numpoints + i,jv*numpoints + j])

                sigma[:,iu*numpoints + i,jv*numpoints + j] = svec.T

    xsize,ysize = numelemsu*numpoints,numelemsv*numpoints
#    print(nanArray)
    # Recomputing NaN values
    for nanA in nanArray:
        print("Recomputing NaN values")
        iu = nanA[0]
        jv = nanA[1]
        
        sumSigma = np.zeros((1,3))
        numneigh = 0
        
        neighbors = [[iu-1,jv-1],
                     [iu-1,jv],
                     [iu-1,jv+1],
                     [iu,jv-1],
                     [iu,jv+1],
                     [iu+1,jv-1],
                     [iu+1,jv],
                     [iu+1,jv+1]]
        
        for neigh in neighbors:
            if neigh[0] >= 0 and neigh[0] <= xsize - 1:
                if neigh[1] >= 0 and neigh[1] <= ysize - 1:
                    if not np.any(np.isnan(sigma[:,neigh[0],neigh[1]])):
                        sumSigma += sigma[:,neigh[0],neigh[1]]
                        numneigh += 1
        
        sigma[:,iu,jv] = sumSigma/numneigh

    return sigma

def plotDisplacementFields(cx,cy,ux,uy):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    xlength = np.amax(np.absolute(cx))
    ylength = np.amax(np.absolute(cy))

    aspectRatio = xlength/ylength

    if aspectRatio > 1.5:
        fig, (ax1,ax2) = plt.subplots(2,1,sharex='col',sharey='row')
    else:
        fig, (ax1,ax2) = plt.subplots(1,2,sharex='col',sharey='row')

    fig.suptitle('Displacement field components')
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    field1 = ax1.pcolormesh(cx,cy,ux,vmin=ux.min(),vmax=ux.max())
    ax1.set_title('Ux')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_aspect('equal')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right",size="5%",pad=0.1)
    cb1 = fig.colorbar(field1,cax=cax,label='[m]')

    field2 = ax2.pcolormesh(cx,cy,uy,vmin=uy.min(),vmax=uy.max())
    ax2.set_title('Uy')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_aspect('equal')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right",size="5%",pad=0.1)
    cb2 = fig.colorbar(field2,cax=cax,label='[m]')

    plt.tight_layout()
    plt.show()

def plotStressFields(cx,cy,sx,sy,sxy,svm):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    xlength = np.amax(np.absolute(cx))
    ylength = np.amax(np.absolute(cy))

    aspectRatio = xlength/ylength

    if aspectRatio > 1.5:
        fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex='col',sharey='row')
    else:
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')

    # Uncomment for cantileverBeam.py
    # fig, axs = plt.subplots(4,1,sharex='col',sharey='row')

    # Uncomment for pressureCylinder.py and plateWithHole.py
    # fig, axs = plt.subplots(2,2,sharex='col',sharey='row')

    fig.suptitle('Stress field components')
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    field1 = ax1.pcolormesh(cx,cy,sx,vmin=sx.min(),vmax=sx.max())
    ax1.set_title('Sx')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_aspect('equal')
    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right",size="5%",pad=0.1)
    cb1 = fig.colorbar(field1,cax=cax1,label='[Pa]')

    field2 = ax2.pcolormesh(cx,cy,sy,vmin=sy.min(),vmax=sy.max())
    ax2.set_title('Sy')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_aspect('equal')
    divider = make_axes_locatable(ax2)
    cax2 = divider.append_axes("right",size="5%",pad=0.1)
    cb2 = fig.colorbar(field2,cax=cax2,label='[Pa]')

    field3 = ax3.pcolormesh(cx,cy,sxy,vmin=sxy.min(),vmax=sxy.max())
    ax3.set_title('Sxy')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_aspect('equal')
    divider = make_axes_locatable(ax3)
    cax3 = divider.append_axes("right",size="5%",pad=0.1)
    cb3 = fig.colorbar(field3,cax=cax3,label='[Pa]')

    field4 = ax4.pcolormesh(cx,cy,svm,vmin=svm.min(),vmax=svm.max())
    ax4.set_title('Von Mises stress')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_aspect('equal')
    divider = make_axes_locatable(ax4)
    cax4 = divider.append_axes("right",size="5%",pad=0.1)
    cb4 = fig.colorbar(field4,cax=cax4,label='[Pa]')

    plt.tight_layout()
    plt.show()

def postProcessing(U,V,p,q,P,D,w,dtot,surfaceprep,matprop):
    elementcorners = surfaceprep[2]
    numelems = len(elementcorners)
#    numelems = nodeselem.shape[0]

    if numelems < 5:
        numpoints = 11
    elif numelems >= 5 and numelems < 10:
        numpoints = 9
    elif numelems >= 10 and numelems < 20:
        numpoints = 7
    else:
        numpoints = 5

    upts,cpts = displacementField(numpoints,U,V,p,q,P,D,w,surfaceprep)
    sigmapts = stressField(numpoints,U,V,p,q,P,w,dtot,matprop,surfaceprep)
    svm = np.sqrt( sigmapts[0,:,:]**2 - 2*sigmapts[0,:,:]*sigmapts[1,:,:] + sigmapts[1,:,:]**2 + 3*sigmapts[2,:,:]**2 )
#    plts.plotting2DField(cpts[0,:,:],cpts[1,:,:],upts[0,:,:],P,["Ux Displacement Field","[m]"])
    plts.plotting2DField(cpts[0,:,:],cpts[1,:,:],sigmapts[0,:,:],P,["Sx Stress Field","[Pa]"])
#    plotDisplacementFields(cpts[0,:,:],cpts[1,:,:],upts[0,:,:],upts[1,:,:])
#    plotStressFields(cpts[0,:,:],cpts[1,:,:],sigmapts[0,:,:],sigmapts[1,:,:],sigmapts[2,:,:],svm)
