# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

# Local project
import src.basisFunctions as bfunc

def binomial(a,b):
    bc = 1.0
    for j in range(1,b+1):
        bc *= ((a+1-j)//j)

    return bc

#######################################################
################# CLASS DEFINITION ####################
#######################################################

class NURBSCurve:
    """
    A class that represent a nurbs curve
    """
    def __init__(self,P,w,p,U=None):
        """
        Initialize the nurbs object with the control points
        and their respective weights,the degree of the spline,
        and the knot vector

        """
        self.P = P
        self.w = w
        self.p = p

        if U is None:
            self.U = bfunc.generateUniformKnotVector(P.shape[0],p)
        else:
            self.U = U

    # End constructor method

    def createCurve(self):
        """
        Create a nurbs curve for further plotting

        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)

        self.cpts = np.zeros((numpoints,2))

        mu = len(self.U) - 1
        nu = mu - self.p - 1
        idx = np.arange(0,self.p+1)

        for i in range(len(urank)):
            uspan = bfunc.findKnotInterval(nu,self.p,urank[i],self.U)
            idxU = uspan + idx - self.p
            nbas = bfunc.basisFunction(uspan,urank[i],mu,self.p,self.U)

            nbas = np.reshape(nbas,(1,len(nbas)))

            ratFunc = (nbas*self.w[idxU,:].T)/(nbas@self.w[idxU,:])

            self.cpts[i,:] = ratFunc@self.P[idxU,:]
        # End for loop
        return self.cpts

    def createTangentCurve(self):
        """
        Create a nurbs tangent curve for further plotting

        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)

        self.cppts = np.zeros((numpoints,2))

        Pw = weightedControlPoints(self.P,self.w)

        mu = len(self.U) - 1
        nu = mu - self.p - 1
        idx = np.arange(0,self.p+1)

        # Derivative order
        d = 1

        for i in range(len(urank)):
            uspan = bfunc.findKnotInterval(nu,self.p,urank[i],self.U)
            idxU = uspan + idx - self.p
            nbas = bfunc.basisFunction(uspan,urank[i],mu,self.p,self.U)
            dnbasU = bfunc.derBasisFunction(uspan,urank[i],mu,self.p,self.U,d)

            # Hughes' way
            Aders = dnbasU*Pw[idxU,-1].T
            wders = dnbasU@Pw[idxU,-1]
            dRatdU = univariateRationalDerivative(Aders,wders,d)
            Ck = dRatdU@self.P[idxU,:]

            self.cppts[i,:] = Ck[d,:]
        # End for loop
        return self.cppts

    def plotCurve(self):
        """
        Plot the curve

        """
        fig,ax = plt.subplots()
        plt.plot(self.cpts[:,0],self.cpts[:,1])
        ax.set_aspect('equal','box')
        plt.plot(self.P[:,0],self.P[:,1],'ro')
        plt.plot(self.P[:,0],self.P[:,1])
        plt.show()

    def plotTangentCurve(self):
        """
        Plot the tangent curve

        """
        fig = plt.figure()
        plt.plot(self.P[:,0],self.P[:,1],'ro')
        plt.plot(self.P[:,0],self.P[:,1])
        plt.plot(self.cpts[:,0],self.cpts[:,1])
        plt.quiver(self.cpts[:,0],self.cpts[:,1],self.cppts[:,0],self.cppts[:,1],color=['k'])
        plt.show()

class NURBSSurface:
    """
    A class that represent a nurbs surface
    """
    def __init__(self,P,w,p,q,U=None,V=None,gridsize=None):
        """
        Initialize the nurbs object with the control points
        and their respective weights,the degree of the spline,
        and the knot vector for both parametric directions
        """
        self.P = P
        self.w = w
        self.p = p
        self.q = q

        if U is None:
            self.U = bfunc.generateUniformKnotVector(gridsize[0],p)
            self.V = bfunc.generateUniformKnotVector(gridsize[1],q)
        else:
            self.U = U
            self.V = V

    # End constructor method

    def retrieveSurfaceInformation(self):
        """
        Getter method for the surface class
        """
        return self.U,self.V,self.p,self.q,self.P,self.w

    def updateSurfaceInformation(self,U,V,p,q,P,w):
        """
        Setter method for the surface class
        """
        self.U = U
        self.V = V
        self.p = p
        self.q = q
        self.P = P
        self.w = w

    def pointInSurface(self,upt,vpt):
        mu = len(self.U) - 1
        mv = len(self.V) - 1
        nu = mu - self.p - 1
        nv = mv - self.q - 1

        Pwl = weightedControlPoints(self.P,self.w)
        Pw = listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        uspan = bfunc.findKnotInterval(nu,self.p,upt,self.U)
        vspan = bfunc.findKnotInterval(nv,self.q,vpt,self.V)

        idR = nonZeroIndicesSurface(uspan,vspan,self.p,self.q,nu)

        R = bivariateRationalFunction(mu,mv,self.p,self.q,uspan,vspan,upt,vpt,self.U,self.V,Pw)
        spt = R@self.P[idR,:]
        return spt

    def createSurface(self,numpoints=41):
        """
        Create a nurbs surface for further plotting
        """
        # numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)
        vrank = np.linspace(self.V.min(),self.V.max(),numpoints)

        Pwl = weightedControlPoints(self.P,self.w)
        Pw = listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        mu = len(self.U) - 1
        mv = len(self.V) - 1
        nu = mu - self.p - 1
        nv = mv - self.q - 1

        idxu = np.arange(0,self.p+1)
        idxv = np.arange(0,self.q+1)

        self.cpts = np.zeros((self.P.shape[1],len(urank),len(vrank)))

        for j in range(len(vrank)):
            for i in range(len(urank)):
                uspan = bfunc.findKnotInterval(nu,self.p,urank[i],self.U)
                vspan = bfunc.findKnotInterval(nv,self.q,vrank[j],self.V)

                idR = nonZeroIndicesSurface(uspan,vspan,self.p,self.q,nu)

                R = bivariateRationalFunction(mu,mv,self.p,self.q,uspan,vspan,urank[i],vrank[j],self.U,self.V,Pw)
                S = R@self.P[idR,:]

                self.cpts[:,i,j] = S
            # End i loop
        # End j loop
        return self.cpts

    def createTangentSurface(self):
        """
        Create a nurbs tangent surface for further plotting
        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)
        vrank = np.linspace(self.V.min(),self.V.max(),numpoints)

        Pwl = weightedControlPoints(self.P,self.w)
        Pw = listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        mu = len(self.U) - 1
        mv = len(self.V) - 1
        nu = mu - self.p - 1
        nv = mv - self.q - 1

        idxu = np.arange(0,self.p+1)
        idxv = np.arange(0,self.q+1)

        self.cpu = np.zeros((self.P.shape[1],len(urank),len(vrank)))
        self.cpv = np.zeros((self.P.shape[1],len(urank),len(vrank)))

        for j in range(len(vrank)):
            for i in range(len(urank)):
                uspan = bfunc.findKnotInterval(nu,self.p,urank[i],self.U)
                vspan = bfunc.findKnotInterval(nv,self.q,vrank[j],self.V)

                idR = nonZeroIndicesSurface(uspan,vspan,self.p,self.q,nu)

                Ralph = bivariateRationalGradient(mu,mv,self.p,self.q,uspan,vspan,urank[i],vrank[j],self.U,self.V,Pw)
                dS = Ralph@self.P[idR,:]

                self.cpu[:,i,j] = dS[1,:]
                self.cpv[:,i,j] = dS[2,:]
            # End i loop
        # End j loop

        return self.cpu,self.cpv

    def createBoundary(self):
        """
        Create a boundary of the nurbs surface for further plotting
        """
        Pwl = weightedControlPoints(self.P,self.w)
        Pw = listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        mu = len(self.U) - 1
        mv = len(self.V) - 1
        nu = mu - self.p - 1
        nv = mv - self.q - 1

        idxu = np.arange(0,self.p+1)
        idxv = np.arange(0,self.q+1)

        boundarycoor = []

        boundlimits = [[np.array([0.0,0.0]),np.array([1.0,0.0])],
                      [np.array([1.0,0.0]),np.array([1.0,1.0])],
                      [np.array([1.0,1.0]),np.array([0.0,1.0])],
                      [np.array([0.0,1.0]),np.array([0.0,0.0])]]

        numpt = 11
        for bndlim in boundlimits:
            parampath = np.linspace(bndlim[0],bndlim[1],numpt,endpoint=True)
            coor = np.zeros((parampath.shape[0],2))
            ipath = 0
            for ppath in parampath:
                uspan = bfunc.findKnotInterval(nu,self.p,ppath[0],self.U)
                vspan = bfunc.findKnotInterval(nv,self.q,ppath[1],self.V)
                idR = nonZeroIndicesSurface(uspan,vspan,self.p,self.q,nu)

                R = bivariateRationalFunction(mu,mv,self.p,self.q,uspan,vspan,ppath[0],ppath[1],self.U,self.V,Pw)
                S = R@self.P[idR,:]
                coor[ipath,:] = S
                ipath += 1

            boundarycoor.append(coor)

        for bc in range(len(boundarycoor)):
            if bc == 0:
                self.boundarycoor1 = boundarycoor[bc]
            else:
                self.boundarycoor1 = np.vstack((self.boundarycoor1,boundarycoor[bc]))

        return self.boundarycoor1

    def plotSurface(self):
        """
        Plot the surface
        """
        cx = self.cpts[0,:,:]
        cy = self.cpts[1,:,:]
        cz = self.cpts[2,:,:]

        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        # ax.contour3D(cx,cy,cz,cmap='viridis')
        ax.plot_surface(cx,cy,cz,cmap='viridis')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

    def plotTangentSurface(self,component):
        """
        Plot the tangent surface
        """
        cx = self.cpts[0,:,:]
        cy = self.cpts[1,:,:]
        cz = self.cpts[2,:,:]

        if component == "u":
            cpx = self.cpu[0,:,:]
            cpy = self.cpu[1,:,:]
            cpz = self.cpu[2,:,:]
        else:
            cpx = self.cpv[0,:,:]
            cpy = self.cpv[1,:,:]
            cpz = self.cpv[2,:,:]

        fig = plt.figure()
        ax = plt.axes(projection = '3d')
        # ax.contour3D(cx,cy,cz,cmap='viridis')
        ax.plot_surface(cx,cy,cz,cmap='viridis')
        plt.quiver(cx,cy,cz,cpx,cpy,cpz,color=['k'],length = 0.01,normalize = True)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

class MultiPatchNURBSSurface():
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
    # End function

    def retrieveSurfaceInformation(self):
        """
        Getter method for the multipatch surface class
        """
        return self.multiU,self.multiV,self.multip,self.multiq,self.multiP,self.multiw,self.globalPatchIndices
    # End function

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
    # End function

    def createFullControlPolygon(self):
        for i in range(len(self.multiP)):
            if i == 0:
                joinedP = self.multiP[i]
                joinedw = self.multiw[i]
            else:
                joinedP = np.vstack((joinedP,self.multiP[i]))
                joinedw = np.vstack((joinedw,self.multiw[i]))
            # End if
        # End for loop
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
        #End for loop
    # End function

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

            Pwl = weightedControlPoints(Pi,wi)
            Pw = listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            mui = len(Ui) - 1
            mvi = len(Vi) - 1
            nui = mui - pi - 1
            nvi = mvi - qi - 1

            cpts = np.zeros((Pi.shape[1],len(urank),len(vrank)))

            for j in range(len(vrank)):
                for i in range(len(urank)):
                    uspan = bfunc.findKnotInterval(nui,pi,urank[i],Ui)
                    vspan = bfunc.findKnotInterval(nvi,qi,vrank[j],Vi)

                    idR = nonZeroIndicesSurface(uspan,vspan,pi,qi,nui)

                    R = bivariateRationalFunction(mui,mvi,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                    S = R@Pi[idR,:]

                    cpts[:,i,j] = S
                # End i loop
            # End j loop

            self.fullcpts.append(cpts)
        # End patch loop

        return self.fullcpts
    # End function

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

            Pwl = weightedControlPoints(Pi,wi)
            Pw = listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

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

                    idR = nonZeroIndicesSurface(uspan,vspan,pi,qi,nui)

                    R = bivariateRationalFunction(mui,mvi,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                    S = R@Pi[idR,:]

                    Ralph = bivariateRationalGradient(mui,mvi,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                    dS = Ralph@Pi[idR,:]

                    # cpts[:,i,j] = S
                    cpu[:,i,j] = dS[1,:]
                    cpv[:,i,j] = dS[2,:]

            self.fullcpu.append(cpu)
            self.fullcpv.append(cpv)

        return self.fullcpu,self.fullcpv
    # End function

    def plotMultipatchSurface(self):
        """
        Plot the multipatch surface
        """
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
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
    # End function

    def plotMultipatchTangentSurface(self,component):
        """
        Plot the multipatch tangent surface
        """
        fig = plt.figure()
        ax = plt.axes(projection = '3d')
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
    # End function
# End class

#####################################################
################# CONTROL POINTS ####################
#####################################################

#Convert from real control points to homogeneous ones
def weightedControlPoints(P,w):
    Pw = np.hstack((P,np.ones((P.shape[0],1))))
    Pw *= w
    return Pw

#Convert from homogeneous to real projection
def geometricControlPoints(Pw):
    w = Pw[:,-1]
    w = np.reshape(w,(len(w),1))
    P = Pw[:,:-1]/w
    return P,w

#Convert list of control points to a spatial grid
def listToGridControlPoints(Pl,U,V,p,q):
    #Number of control points in the U direction
    NU = len(U) - p - 1

    #Number of control points in the V direction
    NV = len(V) - q - 1

    # print("NU:",NU)
    # print("NV:",NV)

    Pg = np.zeros((Pl.shape[1],NU,NV))

    for i in range(0,Pl.shape[1]):
        Pg[i,:,:] = np.reshape(Pl[:,i],(NU,NV),order='F')

    return Pg

def gridToListControlPoints(Pg):
    Pl = np.zeros((Pg.shape[1]*Pg.shape[2],Pg.shape[0]))

    for i in range(0,Pg.shape[0]):
        Pl[:,i] = np.reshape(Pg[i,:,:],(Pg.shape[1]*Pg.shape[2]),order='F')

    return Pl

#########################################################################
###################### RATIONAL BASIS FUNCTIONS #########################
#########################################################################

def nonZeroIndicesSurface(uspan,vspan,p,q,nu):
    idr = []

    for l in range(0,q+1):
        for k in range(0,p+1):
            i = (vspan-q+l)*(nu+1) + (uspan-p+k)
            idr.append(i)

    return idr

def rationalCurveDerivative(aders,wders,d):
    Ck = np.zeros((d+1,2))

    for k in range(0,d+1):
        v = aders[k,:]
        for i in range(1,k+1):
            v -= binomial(k,i)*wders[i]*Ck[k-i,:]

        Ck[k,:] = v/wders[0]

    return Ck

def univariateRationalDerivative(aders,wders,d):
    drat = np.zeros((d+1,aders.shape[1]))

    for k in range(0,d+1):
        v = aders[k,:]
        for i in range(1,k+1):
            v -= binomial(k,i)*wders[i]*drat[k-i,:]

        drat[k,:] = v/wders[0]

    return drat

def bivariateRationalFunction(mu,mv,p,q,uspan,vspan,u,v,U,V,Pw):
    Nu = bfunc.basisFunction(uspan,u,mu,p,U)
    Nv = bfunc.basisFunction(vspan,v,mv,q,V)

    R = np.zeros((1,(p+1)*(q+1)))
    i = 0

    for l in range(0,q+1):
        for k in range(0,p+1):
            R[0][i] = Nu[k]*Nv[l]*Pw[-1,uspan-p+k,vspan-q+l]
            i += 1

    R /= np.sum(R)

    return R

def bivariateRationalGradient(mu,mv,p,q,uspan,vspan,u,v,U,V,Pw):
    # Derivative order
    d = 1

    # The first row for dNu and dNv has the shape functions
    # Nu and Nv respectively

    dNu = bfunc.derBasisFunction(uspan,u,mu,p,U,d)
    dNv = bfunc.derBasisFunction(vspan,v,mv,q,V,d)

    # The first row for dR has the shape function
    # The second row for dR has the derivative w.r.t u
    # The third row for dR has the derivative w.r.t v
    dR = np.zeros((3,(p+1)*(q+1)))

    # The first row for dA and dW has the derivative w.r.t u
    # The second row for dA and dW has the derivative w.r.t v
    dA = np.zeros((2,(p+1)*(q+1)))
#    dW = np.zeros((2,1))
    i = 0

    for l in range(0,q+1):
        for k in range(0,p+1):
            iu = uspan - p + k
            jv = vspan - q + l
            dR[0][i] = dNu[0][k]*dNv[0][l]*Pw[-1,iu,jv]
            dA[0][i] = dNu[d][k]*dNv[0][l]*Pw[-1,iu,jv]
            dA[1][i] = dNu[0][k]*dNv[d][l]*Pw[-1,iu,jv]
            i += 1

    W = np.sum(dR[0,:])
    dW = np.sum(dA,axis = 1)
    biN = dR[0,:]/W

    dR[0,:] = biN
    dR[1,:] = (dA[0,:] - dW[0]*biN)/W
    dR[2,:] = (dA[1,:] - dW[1]*biN)/W

    return dR
