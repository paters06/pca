import numpy as np
import matplotlib.pyplot as plt

import src.spline_functions.basisFunctions as bfunc
from src.spline_functions.nurbs import NURBSObject

class NURBSCurve(NURBSObject):
    """
    A class that represent a nurbs curve
    """
    def __init__(self, P: np.ndarray, w: np.ndarray, p: int, U: np.ndarray|None) -> None:
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
            nbas = bfunc.basisFunction(uspan,urank[i],self.p,self.U)

            nbas = np.reshape(nbas,(1,len(nbas)))

            ratFunc = (nbas*self.w[idxU,:].T)/(nbas@self.w[idxU,:])

            self.cpts[i,:] = ratFunc@self.P[idxU,:]
        return self.cpts

    def createTangentCurve(self):
        """
        Create a nurbs tangent curve for further plotting

        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)

        self.cppts = np.zeros((numpoints,2))

        Pw = self.weightedControlPoints(self.P,self.w)

        mu = len(self.U) - 1
        nu = mu - self.p - 1
        idx = np.arange(0,self.p+1)

        # Derivative order
        d = 1

        for i in range(len(urank)):
            uspan = bfunc.findKnotInterval(nu,self.p,urank[i],self.U)
            idxU = uspan + idx - self.p
            # nbas = bfunc.basisFunction(uspan,urank[i],self.p,self.U)
            dnbasU = bfunc.derBasisFunction(uspan,urank[i],mu,self.p,self.U,d)

            # Hughes' way
            Aders = dnbasU*Pw[idxU,-1].T
            wders = dnbasU@Pw[idxU,-1]
            dRatdU = self.univariateRationalDerivative(Aders,wders,d)
            Ck = dRatdU@self.P[idxU,:]

            self.cppts[i,:] = Ck[d,:]
        return self.cppts

    def plotBasisFunctions(self) -> None:
        """
        Put legend outside the plot
        https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot/4700762#4700762
        """
        numpoints = 81
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)
        mu = len(self.U) - 1
        nu = mu - self.p - 1
        # idx = np.arange(0,self.p+1)
        N_bas = np.zeros((nu+1,numpoints))

        for i in range(len(urank)):
            # uspan = bfunc.findKnotInterval(nu,self.p,urank[i],self.U)
            for j in range(0, nu+1):
                Npi = bfunc.oneBasisFunction(self.p, self.U, j, urank[i])
                N_bas[j,i] = Npi
        
        fig, ax = plt.subplots()
        for j in range(0, nu+1):
            label_i = "N("+str(j)+","+str(self.p)+")"
            plt.plot(urank, N_bas[j,:], label=label_i)
        fig.legend(loc=5)
        # plt.subplots_adjust(right=0.8)
        fig.tight_layout(rect=(0., 0., 0.85, 1.))
        plt.show()

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
        # fig = plt.figure()
        plt.plot(self.P[:,0],self.P[:,1],'ro')
        plt.plot(self.P[:,0],self.P[:,1])
        plt.plot(self.cpts[:,0],self.cpts[:,1])
        plt.quiver(self.cpts[:,0],self.cpts[:,1],self.cppts[:,0],self.cppts[:,1],color=['k'])
        plt.show()