# Python libraries
import numpy as np
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
    def __init__(self,U,p,P,w):
        """
        Initialize the nurbs object with the control points
        and their respective weights,the degree of the spline,
        and the knot vector
        
        """
        self.U = U
        self.p = p
        self.P = P
        self.w = w
    
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
    def __init__(self,U,V,p,q,P,w):
        """
        Initialize the nurbs object with the control points
        and their respective weights,the degree of the spline,
        and the knot vector for both parametric directions
        
        """
        self.U = U
        self.V = V
        self.p = p
        self.q = q
        self.P = P
        self.w = w
    
    def createSurface(self):
        """
        Create a nurbs surface for further plotting
        
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
    
#    print("NU:",NU)
#    print("NV:",NV)

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

################## NURBS CURVES ######################

def nurbsCurve(U,p,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)

    cpts = np.zeros((numpoints,2))

    mu = len(U) - 1
    nu = mu - p - 1
    idx = np.arange(0,p+1)

    for i in range(len(urank)):
        uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
        idxU = uspan + idx - p
        nbas = bfunc.basisFunction(uspan,urank[i],mu,p,U)
        
        nbas = np.reshape(nbas,(1,len(nbas)))
        
        ratFunc = (nbas*w[idxU,:].T)/(nbas@w[idxU,:])
        
        cpts[i,:] = ratFunc@P[idxU,:]
    return cpts

def nurbsCurveTangent(U,p,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)

    cppts = np.zeros((numpoints,2))
    
    Pw = weightedControlPoints(P,w)

    mu = len(U) - 1
    nu = mu - p - 1
    idx = np.arange(0,p+1)
    
    # Derivative order
    d = 1
    
    for i in range(len(urank)):
        uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
        idxU = uspan + idx - p
        nbas = bfunc.basisFunction(uspan,urank[i],mu,p,U)
        dnbasU = bfunc.derBasisFunction(uspan,urank[i],mu,p,U,d)
        
        """Just using the formulas way"""
#        Aders = dnbasU*Pw[idxU,-1].T
#        wders = dnbasU@Pw[idxU,-1]
#        wders = np.reshape(wders,(len(wders),1))
#        dRatdU = (Aders - (wders*nbas))/wders[0]
#        Ck = dRatdU@P[idxU,:]
        
        """The NURBS Book way"""
        
#        dCw = dnbasU@Pw[idxU,:]
        # Selecting colums from 0 to the previous to the last one
#        Aders = dCw[:,0:-1]
        # Selecting the last column
#        wders = dCw[:,-1]
#        Ck = rationalCurveDerivative(Aders,wders,d)

        """Hughes' way"""
        Aders = dnbasU*Pw[idxU,-1].T
        wders = dnbasU@Pw[idxU,-1]
        dRatdU = univariateRationalDerivative(Aders,wders,d)
        Ck = dRatdU@P[idxU,:]
        
        cppts[i,:] = Ck[d,:]
    
    return cppts

################# NURBS SURFACES #####################

def nurbs2DBoundary(U,V,p,q,P,w):
    Pwl = weightedControlPoints(P,w)
    Pw = listToGridControlPoints(Pwl,U,V,p,q)
    
    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1
    
    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)

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
            uspan = bfunc.findKnotInterval(nu,p,ppath[0],U)
            vspan = bfunc.findKnotInterval(nv,q,ppath[1],V)
            idR = nonZeroIndicesSurface(uspan,vspan,p,q,nu)
            
            R = bivariateRationalFunction(mu,mv,p,q,uspan,vspan,ppath[0],ppath[1],U,V,Pw)
            S = R@P[idR,:]
            coor[ipath,:] = S
            ipath += 1

        boundarycoor.append(coor)

    for bc in range(len(boundarycoor)):
        if bc == 0:
            boundarycoor1 = boundarycoor[bc]
        else:
            boundarycoor1 = np.vstack((boundarycoor1,boundarycoor[bc]))

    return boundarycoor1

def nurbsSurface(U,V,p,q,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)
    
    Pwl = weightedControlPoints(P,w)
    Pw = listToGridControlPoints(Pwl,U,V,p,q)
    
    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1
    
    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)
    
    cpts = np.zeros((P.shape[1],len(urank),len(vrank)))
    
    for j in range(len(vrank)):
        for i in range(len(urank)):
            uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
            vspan = bfunc.findKnotInterval(nv,q,vrank[j],V)
            
            idR = nonZeroIndicesSurface(uspan,vspan,p,q,nu)
            
            R = bivariateRationalFunction(mu,mv,p,q,uspan,vspan,urank[i],vrank[j],U,V,Pw)
            S = R@P[idR,:]

            cpts[:,i,j] = S
    return cpts

def nurbsSurfaceTangent(U,V,p,q,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    Pwl = weightedControlPoints(P,w)
    Pw = listToGridControlPoints(Pwl,U,V,p,q)
    
    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1
    
    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)

    cpu = np.zeros((P.shape[1],len(urank),len(vrank)))
    cpv = np.zeros((P.shape[1],len(urank),len(vrank)))
    
    for j in range(len(vrank)):
        for i in range(len(urank)):
            uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
            vspan = bfunc.findKnotInterval(nv,q,vrank[j],V)
            
            idR = nonZeroIndicesSurface(uspan,vspan,p,q,nu)
            
            Ralph = bivariateRationalGradient(mu,mv,p,q,uspan,vspan,urank[i],vrank[j],U,V,Pw)
            dS = Ralph@P[idR,:]

            cpu[:,i,j] = dS[1,:]
            cpv[:,i,j] = dS[2,:]
    
    return cpu,cpv

################# MULTIPATCH NURBS SURFACES #####################

def multipatchNurbsSurface(mulU,mulV,mulp,mulq,fullP,fullw,idctrlpts):
    numpoints = 11
    
    numpatches = len(mulU)
    
    fullcpts = []
    
    for ipatch in range(0,numpatches):
        Ui = mulU[ipatch]
        Vi = mulV[ipatch]
        
        pi = mulp[ipatch]
        qi = mulq[ipatch]
        
        urank = np.linspace(Ui.min(),Ui.max(),numpoints)
        vrank = np.linspace(Vi.min(),Vi.max(),numpoints)
        
        Pi = fullP[idctrlpts[ipatch],:]
        wi = fullw[idctrlpts[ipatch],:]
        
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
        
        fullcpts.append(cpts)
    # End patch loop
    
    return fullcpts

def multipatchNurbsSurfaceTangent(mulU,mulV,mulp,mulq,fullP,fullw,idctrlpts):
    numpoints = 11
    
    numpatches = len(mulU)
    
    fullcpu = []
    fullcpv = []
    
    for ipatch in range(0,numpatches):
        Ui = mulU[ipatch]
        Vi = mulV[ipatch]
        
        pi = mulp[ipatch]
        qi = mulq[ipatch]
        
        urank = np.linspace(Ui.min(),Ui.max(),numpoints)
        vrank = np.linspace(Vi.min(),Vi.max(),numpoints)
        
        Pi = fullP[idctrlpts[ipatch]]
        wi = fullw[idctrlpts[ipatch]]
        
        Pwl = weightedControlPoints(Pi,wi)
        Pw = listToGridControlPoints(Pwl,Ui,Vi,pi,qi)
        
        mui = len(Ui) - 1
        mvi = len(Vi) - 1
        nui = mui - pi - 1
        nvi = mvi - qi - 1
        
#        cpts = np.zeros((Pi.shape[1],len(urank),len(vrank)))
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
                
#                cpts[:,i,j] = S
                cpu[:,i,j] = dS[1,:]
                cpv[:,i,j] = dS[2,:]
        
        fullcpu.append(cpu)
        fullcpv.append(cpv)
    
    return fullcpu,fullcpv
