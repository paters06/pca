import numpy as np
import matplotlib.pyplot as plt

# import src.spline_functions.basisFunctions as bfunc
# from src.spline_functions.basisFunctions import find_knot_interval
# from src.spline_functions.basisFunctions import basis_function
# from src.spline_functions.basisFunctions import der_basis_function
from src.spline_functions.basisFunctions import one_basis_function
from src.spline_functions.nurbs import NURBSObject

# F2py imported modules
# from bspline_basis_functions_pmod import bspline_basis_functions  # noqa: E402
# from nurbs_curve_pymod import nurbs_curve

class NURBSCurve(NURBSObject):
    """
    A class that represent a nurbs curve
    """
    def __init__(self, P: np.ndarray, w: np.ndarray, p: int, U: np.ndarray) -> None:
        """
        Initialize the nurbs object with the control points
        and their respective weights,the degree of the spline,
        and the knot vector

        """
        self.P = P
        self.w = w
        self.p = p
        self.U = U

    @property
    def control_points(self) -> np.ndarray:
        return self.P

    def create_curve(self) -> np.ndarray:
        """
        Create a nurbs curve for further plotting
        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)

        self.cpts = np.zeros((numpoints,2))

        # mu = len(self.U) - 1
        # nu = mu - self.p - 1
        idx = np.arange(0,self.p+1)

        for i in range(len(urank)):
            # uspan = find_knot_interval(nu,self.p,urank[i],self.U)
            uspan = bspline_basis_functions.find_span(self.p, urank[i], self.U)
            idxU = uspan + idx - self.p
            # nbas = basis_function(uspan,urank[i],self.p,self.U)
            nbas = bspline_basis_functions.basis_function(uspan, urank[i], self.p, self.U)

            nbas = np.reshape(nbas,(1,len(nbas)))

            ratFunc = (nbas*self.w[idxU,:].T)/(nbas@self.w[idxU,:])

            self.cpts[i,:] = ratFunc@self.P[idxU,:]
        return self.cpts
    
    def create_curve_fortran(self) -> np.ndarray:
        """
        Create a nurbs curve for further plotting through fortran code
        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)

        self.cpts = np.zeros((numpoints,2))

        Pw = self.weightedControlPoints(self.P,self.w)

        # print(nurbs_curve.curve_point.__doc__)

        for i in range(len(urank)):
            Ci = np.zeros((1,2))
            Ci = nurbs_curve.curve_point(self.p, self.U, Pw, urank[i])
            self.cpts[i,:] = Ci

        return self.cpts

    def create_tangent_curve(self) -> np.ndarray:
        """
        Create a nurbs tangent curve for further plotting

        """
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)

        self.cppts = np.zeros((numpoints,2))

        Pw = self.weightedControlPoints(self.P,self.w)

        # mu = len(self.U) - 1
        # nu = mu - self.p - 1
        idx = np.arange(0,self.p+1)

        # Derivative order
        d = 1

        for i in range(len(urank)):
            # uspan = find_knot_interval(nu,self.p,urank[i],self.U)
            uspan = bspline_basis_functions.find_span(self.p, urank[i], self.U)
            idxU = uspan + idx - self.p
            # nbas = basis_function(uspan,urank[i],self.p,self.U)
            # dnbasU = der_basis_function(uspan,urank[i],mu,self.p,self.U,d)
            dnbasU = bspline_basis_functions.der_basis_functions(uspan, urank[i], self.p, d, self.U)

            # Hughes' way
            Aders = dnbasU*Pw[idxU,-1].T
            wders = dnbasU@Pw[idxU,-1]
            dRatdU = self.univariateRationalDerivative(Aders,wders,d)
            Ck = dRatdU@self.P[idxU,:]

            self.cppts[i,:] = Ck[d,:]
        return self.cppts
    
    def create_tangent_curve_fortran(self) -> np.ndarray:
        numpoints = 41
        urank = np.linspace(self.U.min(),self.U.max(),numpoints)
        
        for i in range(len(urank)):
            Cpi = np.zeros((1,2))
            self.cppts[i,:] = Cpi
        
        self.cppts = np.zeros((numpoints,2))

        return self.cppts

    def plot_basis_functions(self) -> None:
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
            # uspan = find_knot_interval(nu,self.p,urank[i],self.U)
            for j in range(0, nu+1):
                Npi = one_basis_function(self.p, self.U, j, urank[i])
                N_bas[j,i] = Npi
        
        fig, ax = plt.subplots()
        for j in range(0, nu+1):
            label_i = "N("+str(j)+","+str(self.p)+")"
            plt.plot(urank, N_bas[j,:], label=label_i)
        fig.legend(loc=5)
        # plt.subplots_adjust(right=0.8)
        fig.tight_layout(rect=(0., 0., 0.85, 1.))
        plt.show()
