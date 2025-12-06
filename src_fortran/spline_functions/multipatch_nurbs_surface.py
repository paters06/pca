import numpy as np
import matplotlib.pyplot as plt

import src.spline_functions.basisFunctions as bfunc
from src.spline_functions.nurbs import NURBSObject

class MultiPatchNURBSSurface(NURBSObject):
    """
    A class that represents an object conformed
    by multiple patches of nurbs surfaces
    """
    def __init__(self,multiU,multiV,multip,multiq,multiP,multiw):
        """
        Constructor for the multipatch object
        It will be implemented as an class of arrays
        """
        self.multiU = multiU
        self.multiV = multiV
        self.multip = multip
        self.multiq = multiq
        self.multiP = multiP
        self.multiw = multiw
        self.createFullControlPolygon()

    def retrieveSurfaceInformation(self):
        """
        Getter method for the multipatch surface class
        """
        return self.multiU,self.multiV,self.multip,self.multiq,self.multiP,self.multiw,self.globalPatchIndices

    def updateMultiPatchInformation(self,multiU,multiV,multip,multiq,multiP,multiw):
        """
        Setter method for the multipatch surface class
        """
        self.multiU = multiU
        self.multiV = multiV
        self.multip = multip
        self.multiq = multiq
        self.multiP = multiP
        self.multiw = multiw
        self.createFullControlPolygon()

    def createFullControlPolygon(self):
        for i in range(len(self.multiP)):
            if i == 0:
                joinedP = self.multiP[i]
                joinedw = self.multiw[i]
            else:
                joinedP = np.vstack((joinedP,self.multiP[i]))
                joinedw = np.vstack((joinedw,self.multiw[i]))
        
        self.fullP,indices = np.unique(joinedP,axis=0,return_index=True)
        self.fullw = joinedw[indices]

        self.globalPatchIndices = []
        for i in range(len(self.multiP)):
            boolCoincidentRows = self.multiP[i][:,None] == self.fullP
            patchIndices = []
            for j in range(len(boolCoincidentRows)):
                idx_i = list(np.where(boolCoincidentRows[j].all(axis=1))[0])
                patchIndices += idx_i
            self.globalPatchIndices.append(patchIndices)

    def createMultipatchSurface(self):
        """
        Create a surface of multiple nurbs patches for further plotting
        """
        numpoints = 11

        numpatches = len(self.multiU)

        self.fullcpts = []

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            urank = np.linspace(Ui.min(),Ui.max(),numpoints)
            vrank = np.linspace(Vi.min(),Vi.max(),numpoints)

            Pi = self.multiP[ipatch]
            wi = self.multiw[ipatch]

            Pwl = self.weightedControlPoints(Pi,wi)
            Pw = self.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            mui = len(Ui) - 1
            mvi = len(Vi) - 1
            nui = mui - pi - 1
            nvi = mvi - qi - 1

            cpts = np.zeros((Pi.shape[1],len(urank),len(vrank)))

            for j in range(len(vrank)):
                for i in range(len(urank)):
                    uspan = bfunc.findKnotInterval(nui,pi,urank[i],Ui)
                    vspan = bfunc.findKnotInterval(nvi,qi,vrank[j],Vi)

                    idR = self.nonZeroIndicesElement(uspan,vspan,pi,qi,nui)

                    R = self.bivariateRationalFunction(pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                    S = R@Pi[idR,:]

                    cpts[:,i,j] = S

            self.fullcpts.append(cpts)

        return self.fullcpts

    def createMultipatchTangentSurface(self):
        """
        Create a surface of the derivatives over multiple nurbs patches
        for further plotting
        """
        numpoints = 11

        numpatches = len(self.multiU)

        self.fullcpu = []
        self.fullcpv = []

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            urank = np.linspace(Ui.min(),Ui.max(),numpoints)
            vrank = np.linspace(Vi.min(),Vi.max(),numpoints)

            Pi = self.multiP[ipatch]
            wi = self.multiw[ipatch]

            Pwl = self.weightedControlPoints(Pi,wi)
            Pw = self.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            mui = len(Ui) - 1
            mvi = len(Vi) - 1
            nui = mui - pi - 1
            nvi = mvi - qi - 1

            # cpts = np.zeros((Pi.shape[1],len(urank),len(vrank)))
            cpu = np.zeros((Pi.shape[1],len(urank),len(vrank)))
            cpv = np.zeros((Pi.shape[1],len(urank),len(vrank)))

            for j in range(len(vrank)):
                for i in range(len(urank)):
                    uspan = bfunc.findKnotInterval(nui,pi,urank[i],Ui)
                    vspan = bfunc.findKnotInterval(nvi,qi,vrank[j],Vi)

                    idR = self.nonZeroIndicesElement(uspan,vspan,pi,qi,nui)

                    # R = self.bivariateRationalFunction(pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                    # S = R@Pi[idR,:]

                    Ralph = self.bivariateRationalGradient(mui,mvi,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                    dS = Ralph@Pi[idR,:]

                    # cpts[:,i,j] = S
                    cpu[:,i,j] = dS[1,:]
                    cpv[:,i,j] = dS[2,:]

            self.fullcpu.append(cpu)
            self.fullcpv.append(cpv)

        return self.fullcpu,self.fullcpv
    
    def create_path_over_multipatch(self, param_pts, id_patches: list[int]):
        numpatches = len(id_patches)
        numpoints = param_pts.shape[0]

        eval_pts = np.zeros((numpatches*numpoints,2))

        id_point = 0

        for ipatch in id_patches:
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.multiP[ipatch]
            wi = self.multiw[ipatch]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu_sp = mu - pi - 1
            nv_sp = mv - qi - 1

            Pwl = self.weightedControlPoints(Pi,wi)
            Pwi = self.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            # Geometric coordinates
            cpts = np.zeros((1,Pi.shape[1]))

            for i in range(0,numpoints):
                uspan = bfunc.findKnotInterval(nu_sp,pi,param_pts[i,0],Ui)
                vspan = bfunc.findKnotInterval(nv_sp,qi,param_pts[i,1],Vi)

                idR = self.nonZeroIndicesElement(uspan,vspan,pi,qi,nu_sp)
                R = self.bivariateRationalFunction(pi,qi,uspan,vspan,param_pts[i,0],param_pts[i,1],Ui,Vi,Pwi)

                cpts = R@Pi[idR,:]

                eval_pts[id_point,:] = cpts[0]
                id_point += 1
        
        self.param_pts = param_pts
        self.patches_on_path = id_patches

        return eval_pts

    def plotMultipatchSurface(self):
        """
        Plot the multipatch surface
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # ax = plt.axes(projection = '3d')
        # ax.set_aspect('equal','box')
        for ipatch in range(0,len(self.fullcpts)):
            cpts = self.fullcpts[ipatch]
            cx = cpts[0,:,:]
            cy = cpts[1,:,:]
            cz = cpts[2,:,:]
            # ax.contour3D(cx, cy, cz, 50, cmap = 'viridis')
            ax.plot_surface(cx, cy, cz, cmap = 'viridis')

        ax.scatter(self.fullP[:,0],self.fullP[:,1],self.fullP[:,2], color = 'red')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

    def plotMultipatchTangentSurface(self,component):
        """
        Plot the multipatch tangent surface
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # ax = plt.axes(projection = '3d')
        for ipatch in range(0,len(self.fullcpts)):
            cpts = self.fullcpts[ipatch]
            if component == "u":
                cppts = self.fullcpu[ipatch]
            else:
                cppts = self.fullcpv[ipatch]

            cx = cpts[0,:,:]
            cy = cpts[1,:,:]
            cz = cpts[2,:,:]

            cpx = cppts[0,:,:]
            cpy = cppts[1,:,:]
            cpz = cppts[2,:,:]

            ax.plot_surface(cx, cy, cz, cmap = 'viridis')
            plt.quiver(cx,cy,cz,cpx,cpy,cpz,color=['b'],length = 0.01,normalize = True)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()