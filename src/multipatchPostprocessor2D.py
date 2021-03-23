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

def defineNumberOfEvaluationPoints(numelems):
    if numelems < 5:
        numpoints = 11
    elif numelems >= 5 and numelems < 10:
        numpoints = 9
    elif numelems >= 10 and numelems < 20:
        numpoints = 7
    else:
        numpoints = 5

    return numpoints

class SolutionField:
    """
    A class for the multipatch solution field
    """
    def __init__(self,multisurface,D,dtot):
        """
        Constructor of the class
        """
        multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints = multisurface.retrieveSurfaceInformation()
        self.multiU = multiU
        self.multiV = multiV
        self.multip = multip
        self.multiq = multiq
        self.fullP = fullP
        self.fullw = fullw
        self.idcontrolpoints = idcontrolpoints
        self.fullD = D
        self.dtot = dtot

    def displacementFieldGrid(self,surfaceprep):
        """
        Computing the displacement field in grid form
        """
        numpatches = len(self.multiU)

        self.fullcpts = []
        self.fullupts = []

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.fullP[self.idcontrolpoints[ipatch],:]
            wi = self.fullw[self.idcontrolpoints[ipatch],:]
            Di = self.fullD[self.idcontrolpoints[ipatch],:]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu = mu - pi - 1
            nv = mv - qi - 1

            # Extraction of surface preprocessing
            nonzeroctrlpts = surfaceprep[ipatch][0]
            surfacespan = surfaceprep[ipatch][1]
            elementcorners = surfaceprep[ipatch][2]

            numElems = len(elementcorners)
            numelemsu = len(np.unique(Ui)) - 1
            numelemsv = len(np.unique(Vi)) - 1

            numpoints = defineNumberOfEvaluationPoints(numElems)

            Pwl = rbs.weightedControlPoints(Pi,wi)
            Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            # Geometric coordinates
            cpts = np.zeros((Pi.shape[1],numelemsu*numpoints,numelemsv*numpoints))

            # Displacements
            upts = np.zeros((Pi.shape[1],numelemsu*numpoints,numelemsv*numpoints))

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
                        R = rbs.bivariateRationalFunction(mu,mv,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pwi)

                        upts[:,iu*numpoints + i,jv*numpoints + j] = R@Di[idR,:]
                        cpts[:,iu*numpoints + i,jv*numpoints + j] = R@Pi[idR,:]
                    # End i loop
                # End j loop
            # End element loop
            self.fullcpts.append(cpts)
            self.fullupts.append(upts)
        #End patch loop

        return self.fullupts,self.fullcpts

    def displacementFieldList(self,surfaceprep):
        """
        Computing the displacement field in list form
        """
        numpatches = len(self.multiU)

        cptslist = []
        uptslist = []

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.fullP[self.idcontrolpoints[ipatch],:]
            wi = self.fullw[self.idcontrolpoints[ipatch],:]
            Di = self.fullD[self.idcontrolpoints[ipatch],:]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu = mu - pi - 1
            nv = mv - qi - 1

            # Extraction of surface preprocessing
            nonzeroctrlpts = surfaceprep[ipatch][0]
            surfacespan = surfaceprep[ipatch][1]
            elementcorners = surfaceprep[ipatch][2]

            numElems = len(elementcorners)
            numelemsu = len(np.unique(Ui)) - 1
            numelemsv = len(np.unique(Vi)) - 1

            numpoints = defineNumberOfEvaluationPoints(numElems)

            Pwl = rbs.weightedControlPoints(Pi,wi)
            Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            # Geometric coordinates
            cpts = np.zeros((numElems*numpoints*numpoints,Pi.shape[1]))

            # Displacements
            upts = np.zeros((numElems*numpoints*numpoints,Pi.shape[1]))

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

                for j in range(0,numpoints):
                    for i in range(0,numpoints):
                        R = rbs.bivariateRationalFunction(mu,mv,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pwi)

                        upts[ielem*numpoints*numpoints + j*numpoints + i,:] = R@Di[idR,:]
                        cpts[ielem*numpoints*numpoints + j*numpoints + i,:] = R@Pi[idR,:]
                    # End i loop
                # End j loop
            # End element loop
            cptslist.append(cpts)
            uptslist.append(upts)
        #End patch loop

        # Rearranging the full matrix of points
        for ipatch in range(0,numpatches):
            if ipatch == 0:
                self.fullcpts = cptslist[ipatch]
                self.fullupts = uptslist[ipatch]
            else:
                self.fullcpts = np.vstack((self.fullcpts,cptslist[ipatch]))
                self.fullupts = np.vstack((self.fullupts,uptslist[ipatch]))

        return self.fullupts,self.fullcpts

    def stressFieldGrid(self,matprop,surfaceprep):
        """
        Computing the stress field in grid form
        """
        numpatches = len(self.multiU)

        self.fullsigma = []

        # Definition of the material matrix
        E = matprop[0]
        nu = matprop[1]
        rho = matprop[2]
        dMat = elasticMatrix(E,nu)

        paramgrad = np.zeros((2,2))

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.fullP[self.idcontrolpoints[ipatch],:]
            wi = self.fullw[self.idcontrolpoints[ipatch],:]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu = mu - pi - 1
            nv = mv - qi - 1

            # Extraction of surface preprocessing
            nonzeroctrlpts = surfaceprep[ipatch][0]
            surfacespan = surfaceprep[ipatch][1]
            elementcorners = surfaceprep[ipatch][2]

            numElems = len(elementcorners)
            numelemsu = len(np.unique(Ui)) - 1
            numelemsv = len(np.unique(Vi)) - 1

            numpoints = defineNumberOfEvaluationPoints(numElems)

            Pwl = rbs.weightedControlPoints(Pi,wi)
            Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            sigma = np.zeros((3,numelemsu*numpoints,numelemsv*numpoints))
            nanArray = []

            # Global degrees of freedom
            globalDOF = np.zeros(2*len(self.idcontrolpoints[ipatch]),dtype=int)
            dof0 = 2*np.array(self.idcontrolpoints[ipatch])
            dof1 = dof0 + 1
            globalDOF[0::2] = dof0
            globalDOF[1::2] = dof1

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

                # Patch degrees of freedom
                patchDOF = np.zeros(2*len(idR),dtype=int)
                dof0 = 2*np.array(idR)
                dof1 = dof0 + 1
                patchDOF[0::2] = dof0
                patchDOF[1::2] = dof1

                jv = ielem//numelemsu
                iu = ielem%numelemsu

                for j in range(0,numpoints):
                    for i in range(0,numpoints):
                        xpcoor = urank[i]
                        ypcoor = vrank[j]

                        biRatGrad = rbs.bivariateRationalGradient(mu,mv,pi,qi,uspan,vspan,xpcoor,ypcoor,Ui,Vi,Pwi)

                        jac = (biRatGrad[1:3,:]@Pi[idR,:]).T
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

                            svec = dMat@(bmat@self.dtot[globalDOF[patchDOF],:])
                        else:
                            print("Singularity")
                            print([iu*numpoints + i,jv*numpoints + j])
                            svec = np.empty((3,1))
                            svec[:] = np.NaN
                            print(svec.shape)
                            nanArray.append([iu*numpoints + i,jv*numpoints + j])
                        # End if
                        sigma[:,iu*numpoints + i,jv*numpoints + j] = svec.T
                    # End i loop
                # End j loop
            # End element loop
            xsize,ysize = numelemsu*numpoints,numelemsv*numpoints
            # print(nanArray)
            # Recomputing NaN values
            for nanA in nanArray:
                print("Recomputing NaN values")
                iu = nanA[0]
                jv = nanA[1]

                sumSigma = np.zeros((1,3))
                numneigh = 0

                neighbors = [[iu-1,jv-1],[iu-1,jv],
                             [iu-1,jv+1],[iu,jv-1],
                             [iu,jv+1],[iu+1,jv-1],
                             [iu+1,jv],[iu+1,jv+1]]

                for neigh in neighbors:
                    if neigh[0] >= 0 and neigh[0] <= xsize - 1:
                        if neigh[1] >= 0 and neigh[1] <= ysize - 1:
                            if not np.any(np.isnan(sigma[:,neigh[0],neigh[1]])):
                                sumSigma += sigma[:,neigh[0],neigh[1]]
                                numneigh += 1
                            # End if
                        # End if
                    # End if
                # End neighbor loop

                sigma[:,iu,jv] = sumSigma/numneigh
            # End Nan loop
            self.fullsigma.append(sigma)
        # End patch loop

        return self.fullsigma

    def stressFieldList(self,matprop,surfaceprep):
        """
        Computing the stress field in list form
        """
        numpatches = len(self.multiU)

        sigmalist = []

        # Definition of the material matrix
        E = matprop[0]
        nu = matprop[1]
        rho = matprop[2]
        dMat = elasticMatrix(E,nu)

        paramgrad = np.zeros((2,2))

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.fullP[self.idcontrolpoints[ipatch],:]
            wi = self.fullw[self.idcontrolpoints[ipatch],:]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu = mu - pi - 1
            nv = mv - qi - 1

            # Extraction of surface preprocessing
            nonzeroctrlpts = surfaceprep[ipatch][0]
            surfacespan = surfaceprep[ipatch][1]
            elementcorners = surfaceprep[ipatch][2]

            numElems = len(elementcorners)
            numelemsu = len(np.unique(Ui)) - 1
            numelemsv = len(np.unique(Vi)) - 1

            numpoints = defineNumberOfEvaluationPoints(numElems)

            Pwl = rbs.weightedControlPoints(Pi,wi)
            Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            sigma = np.zeros((numElems*numpoints*numpoints,3))
            nanArray = []

            # Global degrees of freedom
            globalDOF = np.zeros(2*len(self.idcontrolpoints[ipatch]),dtype=int)
            dof0 = 2*np.array(self.idcontrolpoints[ipatch])
            dof1 = dof0 + 1
            globalDOF[0::2] = dof0
            globalDOF[1::2] = dof1

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

                # Patch degrees of freedom
                patchDOF = np.zeros(2*len(idR),dtype=int)
                dof0 = 2*np.array(idR)
                dof1 = dof0 + 1
                patchDOF[0::2] = dof0
                patchDOF[1::2] = dof1

                for j in range(0,numpoints):
                    for i in range(0,numpoints):
                        xpcoor = urank[i]
                        ypcoor = vrank[j]

                        biRatGrad = rbs.bivariateRationalGradient(mu,mv,pi,qi,uspan,vspan,xpcoor,ypcoor,Ui,Vi,Pwi)

                        jac = (biRatGrad[1:3,:]@Pi[idR,:]).T
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

                            svec = dMat@(bmat@self.dtot[globalDOF[patchDOF],:])
                        else:
                            print("Singularity")
                            print([iu*numpoints + i,jv*numpoints + j])
                            svec = np.empty((3,1))
                            svec[:] = np.NaN
                            print(svec.shape)
                            nanArray.append([ielem*numpoints*numpoints + j*numpoints + i])
                        # End if
                        sigma[ielem*numpoints*numpoints + j*numpoints + i,:] = svec.T
                    # End i loop
                # End j loop
            # End element loop
            xsize,ysize = numelemsu*numpoints,numelemsv*numpoints
            # print(nanArray)
            nanArray = []
            # Recomputing NaN values
            for nanA in nanArray:
                print("Recomputing NaN values")
                iu = nanA[0]
                jv = nanA[1]

                sumSigma = np.zeros((1,3))
                numneigh = 0

                neighbors = [[iu-1,jv-1],[iu-1,jv],
                             [iu-1,jv+1],[iu,jv-1],
                             [iu,jv+1],[iu+1,jv-1],
                             [iu+1,jv],[iu+1,jv+1]]

                for neigh in neighbors:
                    if neigh[0] >= 0 and neigh[0] <= xsize - 1:
                        if neigh[1] >= 0 and neigh[1] <= ysize - 1:
                            if not np.any(np.isnan(sigma[:,neigh[0],neigh[1]])):
                                sumSigma += sigma[:,neigh[0],neigh[1]]
                                numneigh += 1
                            # End if
                        # End if
                    # End if
                # End neighbor loop

                sigma[:,iu,jv] = sumSigma/numneigh
            # End Nan loop
            sigmalist.append(sigma)
        # End patch loop

        # Rearranging the full matrix of points
        for ipatch in range(0,numpatches):
            if ipatch == 0:
                self.fullsigma = sigmalist[ipatch]
            else:
                self.fullsigma = np.vstack((self.fullsigma,sigmalist[ipatch]))

        return self.fullsigma

    def plotDisplacementFields(self):
        """
        Plot the multipatch displacement field
        FIX THE ROUTINE
        """
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cx = self.fullcpts[:,0]
        cy = self.fullcpts[:,1]
        ux = self.fullupts[:,0]
        uy = self.fullupts[:,1]

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
        """
        Plot the multipatch stress field
        FIX THE ROUTINE
        """
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        cx = self.fullcpts[0,:,:]
        cy = self.fullcpts[1,:,:]
        sx = self.fullsigma[0,:,:]
        sy = self.fullsigma[1,:,:]
        sxy = self.fullsigma[2,:,:]

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


################ MAIN POSTPROCESSING FUNCTION ####################

def postProcessing(multisurface,D,dtot,surfaceprep,matprop):
    solfield = SolutionField(multisurface,D,dtot)
    fullupts,fullcpts = solfield.displacementFieldList(surfaceprep)
    fullsigmapts = solfield.stressFieldList(matprop,surfaceprep)

    # plts.plotMultipatchField(fullcpts,fullupts,0,["Ux Displacement Field","[m]"])
    plts.plotMultipatchField(fullcpts,fullsigmapts,0,["Sx Stress Field","[Pa]"])
    # solfield.plotDisplacementFields()
    # plotStressFields(cpts[0,:,:],cpts[1,:,:],sigmapts[0,:,:],sigmapts[1,:,:],sigmapts[2,:,:],svm)
