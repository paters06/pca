# Python libraries
import math
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
from matplotlib import animation

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs
from src.linearElastoStaticsSolver import elasticMatrix
import src.plottingScripts as plts

################ POSTPROCESSING ####################
class SolutionField:
    """
    A class for the solution field
    """
    def __init__(self,phenomenon,surface,D,dtot):
        """
        Constructor of the class
        """
        U,V,p,q,P,w = surface.retrieveSurfaceInformation()
        self.phenomenon = phenomenon
        self.Usol = U
        self.Vsol = V
        self.psol = p
        self.qsol = q
        self.Psol = P
        self.wsol = w
        self.Dsol = D
        self.dtot = dtot

    def updateControlPoints(self,D_new,dtotal_new):
        self.Dsol = D_new
        self.dtot = dtotal_new

    def temperatureField(self,numpoints,surfaceprep):
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

        # Temperature
        self.tpts = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

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

                    self.tpts[iu*numpoints + i,jv*numpoints + j] = R@self.Dsol[idR,:]
                    self.cpts[:,iu*numpoints + i,jv*numpoints + j] = R@self.Psol[idR,:]

        return self.tpts,self.cpts

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
        self.upts = np.zeros((2,numelemsu*numpoints,numelemsv*numpoints))

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

                    self.upts[:,iu*numpoints + i,jv*numpoints + j] = R@self.Dsol[idR,:]
                    self.cpts[:,iu*numpoints + i,jv*numpoints + j] = R@self.Psol[idR,:]

        return self.upts,self.cpts

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
        rho = matprop[2]
        dMat = elasticMatrix(E,nu)

        paramgrad = np.zeros((2,2))
        self.sigma = np.zeros((3,numelemsu*numpoints,numelemsv*numpoints))
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

                    biRatGrad = rbs.bivariateRationalGradient(mu,mv,self.psol,self.qsol,uspan,vspan,xpcoor,ypcoor,self.Usol,self.Vsol,Pw)

                    jac = (biRatGrad[1:3,:]@self.Psol[idR,:]).T
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

                        svec = dMat@(bmat@self.dtot[globalDOF,:])
                    else:
                        print("Singularity")
                        print([iu*numpoints + i,jv*numpoints + j])
                        svec = np.empty((3,1))
                        svec[:] = np.NaN
                        print(svec.shape)
                        nanArray.append([iu*numpoints + i,jv*numpoints + j])

                    self.sigma[:,iu*numpoints + i,jv*numpoints + j] = svec.T

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
                        if not np.any(np.isnan(self.sigma[:,neigh[0],neigh[1]])):
                            sumSigma += self.sigma[:,neigh[0],neigh[1]]
                            numneigh += 1

            self.sigma[:,iu,jv] = sumSigma/numneigh

        return self.sigma

    def showExtremaValues(self):
        if self.phenomenon == "Elasticity":
            print("Displacements")
            print("UX ==> Max: {:.5f} m. Min: {:.5f} m".format(np.max(self.upts[0,:,:]),np.min(self.upts[0,:,:])))
            print("UY ==> Max: {:.5f} m. Min: {:.5f} m".format(np.max(self.upts[1,:,:]),np.min(self.upts[1,:,:])))
            print("Stresses")
            print("SXX ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.sigma[0,:,:]),np.min(self.sigma[0,:,:])))
            print("SYY ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.sigma[1,:,:]),np.min(self.sigma[1,:,:])))
            print("SXY ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.sigma[2,:,:]),np.min(self.sigma[2,:,:])))
        elif self.phenomenon == "Heat":
            print("Temperature")
            print("T ==> Max: {:.5f} °C. Min: {:.5f} °C".format(np.max(self.tpts),np.min(self.tpts)))
        else:
            print("Check physics selection")

    def plotTemperatureField(self):
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cx = self.cpts[0,:,:]
        cy = self.cpts[1,:,:]
        T = self.tpts

        xlength = np.amax(np.absolute(cx))
        ylength = np.amax(np.absolute(cy))

        aspectRatio = xlength/ylength

        fig,ax1 = plt.subplots()

        # if aspectRatio > 1.5:
        #     fig, (ax1,ax2) = plt.subplots(2,1,sharex='col',sharey='row')
        # else:
        #     fig, (ax1,ax2) = plt.subplots(1,2,sharex='col',sharey='row')

        # fig.suptitle('Temperature field')
        # fig.subplots_adjust(hspace=0.4, wspace=0.4)

        field1 = ax1.pcolormesh(cx,cy,T,vmin=T.min(),vmax=T.max())
        ax1.set_title('Temperature Field')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_aspect('equal')
        # divider = make_axes_locatable(ax1)
        # cax = divider.append_axes("right",size="5%",pad=0.1)
        cb1 = fig.colorbar(field1,label='[°C]')

        # plts.plotting2DField(cpts[0,:,:],cpts[1,:,:],sigmapts[0,:,:],["Sx Stress Field","[Pa]"])

        plt.tight_layout()
        plt.show()

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

def plotTransientField(phenomenon,surface,surfaceprep,matprop,Un,T,dt,savevideo):
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

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig = plt.figure()
    ax = fig.add_subplot(111)

    div = make_axes_locatable(ax)
    cax = div.append_axes('right', '5%', '5%')

    frames = []
    numframes = Un.shape[1]
    print("CALCULATING FIELDS FOR ANIMATION")
    for i in range(numframes):
        solfield = SolutionField(phenomenon,surface,Un[:,i,None],Un[:,i,None])
        tpts,cpts = solfield.temperatureField(numpoints,surfaceprep)
        frames.append(tpts)
        if i == 0:
            xcoor = cpts[0,:,:]
            ycoor = cpts[1,:,:]
        # End if
    # End for loop
    numsteps = int(T/dt + 1)
    tsteps = np.linspace(0,T,numsteps)

    extent = np.min(xcoor), np.max(xcoor), np.min(ycoor), np.max(ycoor)
    print("CALCULATION FINISHED")

    field = ax.imshow(frames[0],origin='lower',extent=extent)
    time_text = ax.text(0.02,1.02,"Time: 0.0 s",transform=ax.transAxes)
    # frame_text = ax.text(0.75,1.02,"Frame # 0",transform=ax.transAxes)
    cb = fig.colorbar(field,cax=cax,label='[°C]')
    vmax = np.max(frames[-1])
    vmin = np.min(frames[0])
    field.set_clim(vmin,vmax)
    tx = ax.set_title('Temperature Field')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    # ax.set_aspect('equal')
    # plt.tight_layout()

    def animate(i):
        arr = frames[i]
        # vmax = np.max(arr)
        # vmin = np.min(arr)
        field.set_data(arr)
        time_text.set_text('Time: %.2f s' % tsteps[i])
        # frame_text.set_text('Frame # %d' % i)
        # field.set_clim(vmin,vmax)
        # return field,

    anim = animation.FuncAnimation(fig,animate,frames=Un.shape[1], interval=10,blit=False)
    if savevideo:
        print("SAVING VIDEO WITH FFMPEG LIBRARY")
        # anim.save('first_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
        writer = animation.FFMpegWriter(fps=30, codec='h264', bitrate=-1)
        anim.save("third_test_impl.mp4", writer=writer)
    else:
        print("SHOWING ANIMATION IN RUNTIME")
        plt.show()

def postProcessing(phenomenon,surface,D,dtot,surfaceprep,matprop):
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

    solfield = SolutionField(phenomenon,surface,D,dtot)
    if phenomenon == "Elasticity":
        upts,cpts = solfield.displacementField(numpoints,surfaceprep)
        sigmapts = solfield.stressField(numpoints,matprop,surfaceprep)
        # plts.plotting2DField(cpts[0,:,:],cpts[1,:,:],upts[0,:,:],["Ux Displacement Field","[m]"])
        # plts.plotting2DField(cpts[0,:,:],cpts[1,:,:],sigmapts[0,:,:],["Sx Stress Field","[Pa]"])

        # solfield.plotDisplacementFields()
        solfield.plotStressFields()
    elif phenomenon == "Heat":
        upts,tpts = solfield.temperatureField(numpoints,surfaceprep)
        solfield.plotTemperatureField()

    solfield.showExtremaValues()
