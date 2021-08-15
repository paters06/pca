# Python libraries
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

# Local project
import src.nurbs as rbs
from src.reissnerMidlinPlateSolver import elasticMatrix_RM
import src.plottingScripts as plts

################ POSTPROCESSING ####################

def fieldSmoother(mat,nanarray,xsize,ysize):
    smoothmat = mat.copy()

    # Recomputing NaN values
    for nanA in nanarray:
        print("Recomputing NaN values")
        iu = nanA[0]
        jv = nanA[1]

        summat = 0
        numneigh = 0

        neighbors = [[iu-1,jv-1],[iu-1,jv],
                    [iu-1,jv+1],[iu,jv-1],
                    [iu,jv+1],[iu+1,jv-1],
                    [iu+1,jv],[iu+1,jv+1]]

        for neigh in neighbors:
            if neigh[0] >= 0 and neigh[0] <= xsize - 1:
                if neigh[1] >= 0 and neigh[1] <= ysize - 1:
                    # if not np.any(np.isnan(self.sigma[:,neigh[0],neigh[1]])):
                    if not np.any(np.isnan(smoothmat[neigh[0],neigh[1]])):
                        # sumSigma += self.sigma[:,neigh[0],neigh[1]]
                        summat += mat[neigh[0],neigh[1]]
                        numneigh += 1
                    # End if
                # End if
            # End if
        # End for loop

        smoothmat[iu,jv] = summat/numneigh
    # End for loop

    return smoothmat
# End function

#######################################################
################# CLASS DEFINITION ####################
#######################################################

class SolutionField:
    """
    A class for the solution field
    """
    def __init__(self,surface,D,dtot):
        """
        Constructor of the class
        """
        U,V,p,q,P,w = surface.retrieveSurfaceInformation()
        self.Usol = U
        self.Vsol = V
        self.psol = p
        self.qsol = q
        self.Psol = P
        self.wsol = w
        self.Dsol = D
        self.dtot = dtot
    # End constructor method

    def updateControlPoints(self,D_new,dtotal_new):
        self.Dsol = D_new
        self.dtot = dtotal_new
    # Enf function
# End parent class

class ElasticitySolution(SolutionField):
    def __init__(self,surface,t,D,dtot):
        SolutionField.__init__(self,surface,D,dtot)
        self.thickness = t
        # self.Dsol = D
        # self.dtot = dtot
    # End constructor method

    def displacementField(self,numpoints,surfaceprep):
        mu = len(self.Usol) - 1
        mv = len(self.Vsol) - 1
        nu = mu - self.psol - 1
        nv = mv - self.qsol - 1

        # Extraction of surface preprocessing
        nonzeroctrlpts = surfaceprep[0]
        surfacespan = surfaceprep[1]
        elementcorners = surfaceprep[2]

        numElems = len(elementcorners)
        numelemsu = len(np.unique(self.Usol)) - 1
        numelemsv = len(np.unique(self.Vsol)) - 1

        Pwl = rbs.weightedControlPoints(self.Psol,self.wsol)
        Pw = rbs.listToGridControlPoints(Pwl,self.Usol,self.Vsol,self.psol,self.qsol)

        # Geometric coordinates
        self.cpts = np.zeros((2,numelemsu*numpoints,numelemsv*numpoints))

        # Displacements
        cartVec = np.zeros(3)
        posVec = np.zeros((3,1))
        z = 0.5*self.thickness
        self.upts = np.zeros((3,numelemsu*numpoints,numelemsv*numpoints))

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
                    R = rbs.bivariateRationalFunction(mu,mv,self.psol,self.qsol,uspan,vspan,urank[i],vrank[j],self.Usol,self.Vsol,Pw)

                    Usol = R@self.Dsol[idR,:]
                    posVec = R@self.Psol[idR,:]
                    # z = (posVec[2] - 0.5*self.thickness)
                    
                    # Usol[1] is x-axis
                    # Usol[2] is y-axis
                    # Usol[0] is z-axis
                    cartVec[0] = -z*Usol[0,1]
                    cartVec[1] = -z*Usol[0,2]
                    cartVec[2] = Usol[0,0]

                    self.upts[:,iu*numpoints + i,jv*numpoints + j] = cartVec
                    self.cpts[:,iu*numpoints + i,jv*numpoints + j] = posVec
                # End for loop
            # End for loop
        # End for loop

        return self.upts,self.cpts
    # End function

    # Improve this function
    def stressField(self,numpoints,matprop,surfaceprep):
        mu = len(self.Usol) - 1
        mv = len(self.Vsol) - 1
        nu = mu - self.psol - 1
        nv = mv - self.qsol - 1

        # Extraction of surface preprocessing
        nonzeroctrlpts = surfaceprep[0]
        surfacespan = surfaceprep[1]
        elementcorners = surfaceprep[2]

        numElems = len(elementcorners)
        numelemsu = len(np.unique(self.Usol)) - 1
        numelemsv = len(np.unique(self.Vsol)) - 1

        Pwl = rbs.weightedControlPoints(self.Psol,self.wsol)
        Pw = rbs.listToGridControlPoints(Pwl,self.Usol,self.Vsol,self.psol,self.qsol)

        # Definition of the material matrix
        E = matprop[0]
        nu = matprop[1]
        t = matprop[4]
        dBMat,dSMat = elasticMatrix_RM(E,nu,t)

        self.sigma = np.zeros((5,numelemsu*numpoints,numelemsv*numpoints))
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
            globalDOF = np.zeros(3*len(idR),dtype=int)
            dof0 = 3*np.array(idR)
            dof1 = dof0 + 1
            dof2 = dof0 + 2
            globalDOF[0::3] = dof0
            globalDOF[1::3] = dof1
            globalDOF[2::3] = dof2

            jv = ielem//numelemsu
            iu = ielem%numelemsu

            for j in range(0,numpoints):
                for i in range(0,numpoints):
                    xpcoor = urank[i]
                    ypcoor = vrank[j]

                    biRatGrad = rbs.bivariateRationalGradient(mu,mv,self.psol,self.qsol,uspan,vspan,xpcoor,ypcoor,self.Usol,self.Vsol,Pw)

                    jac = (biRatGrad[1:3,:]@self.Psol[idR,:]).T
                    detJac = abs(np.linalg.det(jac))

                    if abs(detJac) > 1e-5:
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

                        bendmvec = -dBMat@(bBMat@self.dtot[globalDOF,:])
                        shearqvec = dSMat@(bSMat@self.dtot[globalDOF,:])
                    else:
                        print("Singularity")
                        print([iu*numpoints + i,jv*numpoints + j])
                        bendmvec = np.empty((3,1))
                        shearqvec = np.empty((2,1))
                        bendmvec[:] = np.NaN
                        shearqvec[:] = np.NaN
                        nanArray.append([iu*numpoints + i,jv*numpoints + j])
                    # End if

                    self.sigma[0:3,iu*numpoints + i,jv*numpoints + j] = bendmvec.T
                    self.sigma[3:5,iu*numpoints + i,jv*numpoints + j] = shearqvec.T
                # End for loop
            # End for loop
        # End for loop
        
        if len(nanArray) > 0:
            xsize,ysize = numelemsu*numpoints,numelemsv*numpoints
            for ii in range(0,5):
                self.sigma[ii,:,:] = fieldSmoother(self.sigma[ii,:,:],xsize,ysize)
            # End for loop
        # End if

        return self.sigma
    # End function

    def plotDisplacementFields(self):
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cx = self.cpts[0,:,:]
        cy = self.cpts[1,:,:]
        ux = self.upts[0,:,:]
        uy = self.upts[1,:,:]

        xlength = np.amax(np.absolute(cx))
        ylength = np.amax(np.absolute(cy))

        aspectRatio = xlength/ylength

        if aspectRatio > 1.5:
            fig, (ax1,ax2) = plt.subplots(2,1,sharex='col',sharey='row')
        else:
            fig, (ax1,ax2) = plt.subplots(1,2,sharex='col',sharey='row')
        # End if

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
    # End function

    def plotStressFields(self):
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cx = self.cpts[0,:,:]
        cy = self.cpts[1,:,:]
        sx = self.sigma[0,:,:]
        sy = self.sigma[1,:,:]
        sxy = self.sigma[2,:,:]

        svm = np.sqrt(sx**2 - 2*sx*sy + sy**2 + 3*sxy**2)

        xlength = np.amax(np.absolute(cx))
        ylength = np.amax(np.absolute(cy))

        aspectRatio = xlength/ylength

        if aspectRatio > 1.5:
            fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex='col',sharey='row')
        else:
            fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row')
        # End if

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
    # End function

    def showExtremaValues(self):
        print("Displacements")
        print("UX ==> Max: {:.5f} m. Min: {:.5f} m".format(np.max(self.upts[0,:,:]),np.min(self.upts[0,:,:])))
        print("UY ==> Max: {:.5f} m. Min: {:.5f} m".format(np.max(self.upts[1,:,:]),np.min(self.upts[1,:,:])))
        print("Stresses")
        print("SXX ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.sigma[0,:,:]),np.min(self.sigma[0,:,:])))
        print("SYY ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.sigma[1,:,:]),np.min(self.sigma[1,:,:])))
        print("SXY ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.sigma[2,:,:]),np.min(self.sigma[2,:,:])))
    # End function
# End child class

######################################################
################# OTHER FUNCTIONS ####################
######################################################

def postProcessing(phenomenon,surface,surfaceprep,dsol,matprop=None):
    elementcorners = surfaceprep[2]
    numelems = len(elementcorners)

    if numelems < 5:
        numpoints = 11
    elif numelems >= 5 and numelems < 10:
        numpoints = 9
    elif numelems >= 10 and numelems < 20:
        numpoints = 7
    else:
        numpoints = 5
    # End if

    # solfield = SolutionField(phenomenon,surface,D,dtot)
    if phenomenon == "Elasticity":
        d1 = dsol[0::3]
        d2 = dsol[1::3]
        d3 = dsol[2::3]
        D = np.hstack((d1,d2))
        D = np.hstack((D,d3))

        # print(D)
        thickness = matprop[4]
        elastfield = ElasticitySolution(surface,thickness,D,dsol)
        cartdisppts,geompts = elastfield.displacementField(numpoints,surfaceprep)
        sigmapts = elastfield.stressField(numpoints,matprop,surfaceprep)
        plts.plotting2DField(geompts[0,:,:],geompts[1,:,:],cartdisppts[2,:,:],["Uz Displacement Field","[mm]"])
        plts.plotting2DField(geompts[0,:,:],geompts[1,:,:],sigmapts[0,:,:],["Mxx bending moment","[N-mm]"])

        # elastfield.plotDisplacementFields()
        # elastfield.plotStressFields()
        # elastfield.showExtremaValues()
    # elif phenomenon == "Heat":
    #     heatfield = HeatSolution(surface,dsol,dsol)
    #     tpts,cpts = heatfield.temperatureField(numpoints,surfaceprep)
    #     heatfield.plotTemperatureField()
    #     heatfield.showExtremaValues()
    # elif phenomenon == "Schrodinger":
    #     wavefield = SchrodingerSolution(surface,dsol,dsol)
    #     wfpts,cpts = wavefield.waveFunctionField(numpoints,surfaceprep)
    #     wavefield.plotWaveFunctionField()
    # End if
# End Function
