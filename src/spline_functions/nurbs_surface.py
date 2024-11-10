import numpy as np

# from src.spline_functions.basisFunctions import find_knot_interval
from src.spline_functions.nurbs import NURBSObject

# F2py imported modules
from bspline_basis_functions_pmod import bspline_basis_functions  # noqa: E402

class NURBSSurface(NURBSObject):
    """
    A class that represent a nurbs surface
    """
    def __init__(self, P:np.ndarray, w: np.ndarray, p: int,
                 q: int, U:np.ndarray, V:np.ndarray):
        """
        Initialize the nurbs object with the control points
        and their respective weights,the degree of the spline,
        and the knot vector for both parametric directions
        """
        self.P = P
        self.w = w
        self.p = p
        self.q = q
        self.U = U
        self.V = V

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

    def point_in_surface(self,upt,vpt):
        mu = len(self.U) - 1
        # mv = len(self.V) - 1
        nu = mu - self.p - 1
        # nv = mv - self.q - 1

        Pwl = self.weightedControlPoints(self.P,self.w)
        Pw = self.listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        # uspan = find_knot_interval(nu,self.p,upt,self.U)
        # vspan = find_knot_interval(nv,self.q,vpt,self.V)
        uspan = bspline_basis_functions.find_span(self.p, upt, self.U)
        vspan = bspline_basis_functions.find_span(self.q, vpt, self.V)

        idR = self.nonZeroIndicesElement(uspan,vspan,self.p,self.q,nu)

        R = self.bivariateRationalFunction(self.p,self.q,uspan,vspan,upt,vpt,self.U,self.V,Pw)
        spt = R@self.P[idR,:]
        return spt

    def create_surface(self, numpoints:int=41) -> np.ndarray:
        """
        Create a nurbs surface for further plotting
        """
        # numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)
        vrank = np.linspace(self.V.min(),self.V.max(),numpoints)

        Pwl = self.weightedControlPoints(self.P,self.w)
        Pw = self.listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        mu = len(self.U) - 1
        # mv = len(self.V) - 1
        nu = mu - self.p - 1
        # nv = mv - self.q - 1

        # idxu = np.arange(0,self.p+1)
        # idxv = np.arange(0,self.q+1)

        self.cpts = np.zeros((self.P.shape[1],len(urank),len(vrank)))

        for j in range(len(vrank)):
            for i in range(len(urank)):
                # uspan = find_knot_interval(nu,self.p,urank[i],self.U)
                # vspan = find_knot_interval(nv,self.q,vrank[j],self.V)
                uspan = bspline_basis_functions.find_span(self.p, urank[i], self.U)
                vspan = bspline_basis_functions.find_span(self.q, vrank[i], self.V)

                idR = self.nonZeroIndicesElement(uspan,vspan,self.p,self.q,nu)

                R = self.bivariateRationalFunction(self.p,self.q,uspan,vspan,urank[i],vrank[j],self.U,self.V,Pw)
                S = R@self.P[idR,:]

                self.cpts[:,i,j] = S
        return self.cpts

    def create_tangent_surface(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Create a nurbs tangent surface for further plotting
        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)
        vrank = np.linspace(self.V.min(),self.V.max(),numpoints)

        Pwl = self.weightedControlPoints(self.P,self.w)
        Pw = self.listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        mu = len(self.U) - 1
        mv = len(self.V) - 1
        nu = mu - self.p - 1
        # nv = mv - self.q - 1

        # idxu = np.arange(0,self.p+1)
        # idxv = np.arange(0,self.q+1)

        self.cpu = np.zeros((self.P.shape[1],len(urank),len(vrank)))
        self.cpv = np.zeros((self.P.shape[1],len(urank),len(vrank)))

        for j in range(len(vrank)):
            for i in range(len(urank)):
                # uspan = find_knot_interval(nu,self.p,urank[i],self.U)
                # vspan = find_knot_interval(nv,self.q,vrank[j],self.V)
                uspan = bspline_basis_functions.find_span(self.p, urank[i], self.U)
                vspan = bspline_basis_functions.find_span(self.q, vrank[i], self.V)

                idR = self.nonZeroIndicesElement(uspan,vspan,self.p,self.q,nu)

                Ralph = self.bivariateRationalGradient(mu,mv,self.p,self.q,uspan,vspan,urank[i],vrank[j],self.U,self.V,Pw)
                dS = Ralph@self.P[idR,:]

                self.cpu[:,i,j] = dS[1,:]
                self.cpv[:,i,j] = dS[2,:]
            # End i loop
        # End j loop

        return self.cpu, self.cpv

    def create_boundary(self):
        """
        Create a boundary of the nurbs surface for further plotting
        """
        Pwl = self.weightedControlPoints(self.P,self.w)
        Pw = self.listToGridControlPoints(Pwl,self.U,self.V,self.p,self.q)

        mu = len(self.U) - 1
        # mv = len(self.V) - 1
        nu = mu - self.p - 1
        # nv = mv - self.q - 1

        # idxu = np.arange(0,self.p+1)
        # idxv = np.arange(0,self.q+1)

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
                # uspan = find_knot_interval(nu,self.p,ppath[0],self.U)
                # vspan = find_knot_interval(nv,self.q,ppath[1],self.V)
                uspan = bspline_basis_functions.find_span(self.p, ppath[0], self.U)
                vspan = bspline_basis_functions.find_span(self.q, ppath[1], self.V)
                idR = self.nonZeroIndicesElement(uspan,vspan,self.p,self.q,nu)

                R = self.bivariateRationalFunction(self.p,self.q,uspan,vspan,ppath[0],ppath[1],self.U,self.V,Pw)
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
