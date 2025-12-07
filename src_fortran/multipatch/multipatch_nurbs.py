# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

# Local project
import src.spline_functions.basisFunctions as bfunc
from src.nurbs_curve import NURBSCurve

class MultiPatchNURBSCurve():
    def __init__(self, curve1: NURBSCurve, curve2: NURBSCurve) -> None:
        self.curve1 = curve1
        self.curve2 = curve2
        self.curve_list = [curve1, curve2]
        self.num_patches = len(self.curve_list)

        self._join_knot_vectors()
        self._join_control_points()
        self.cpts = self._create_multipatch_curve()
    
    def _join_knot_vectors(self):
        U1 = self.curve1.U
        U2 = self.curve2.U

        m1 = len(U1)
        p1 = self.curve1.p

        m2 = len(U2)
        p2 = self.curve2.p

        print('U1: {}'.format(U1))
        print('U2: {}'.format(U2))

        w1 = self.curve1.w
        w2 = self.curve2.w

        U1_max = U1.max()
        U2_pre = U2 + U1_max

        self.knot_to_repeat = U1[-1]
        start = U1[:m1-p1-1]
        end = U2_pre[p2+1:m2]
        # joint = self.knot_to_repeat*np.ones(self.curve1.p)
        joint = self.knot_to_repeat

        print(start)
        print(joint)
        print(end)

        start_w = w1[:w1.shape[0]-1,:]
        # print(self.curve2.w)
        end_w = w2[1:w2.shape[0],:]
        joint_w = np.ones((self.curve1.p,1))

        # print(start_w)
        # print(joint_w)
        # print(end_w)

        self.full_U = np.hstack((start, joint, end))
        self.full_p = self.curve1.p
        self.full_w = np.vstack((start_w, joint_w, end_w))
        print(self.full_U)
        # print(self.full_w)

    def _join_control_points(self):
        start = self.curve1.P
        end = self.curve2.P[1:self.curve2.P.shape[0],:]

        self.full_P = np.vstack((start, end))
        print(self.full_P)
    
    def _create_multipatch_curve(self):
        """
        Create a nurbs curve for further plotting
        """
        numpoints = 81
        urank_full = np.linspace(self.full_U.min(),self.full_U.max(),numpoints)

        self.cpts = np.zeros((numpoints,2))

        mu = len(self.full_U) - 1
        nu = mu - self.full_p - 1
        idx = np.arange(0,self.full_p+1)

        for i in range(len(urank_full)):
            uspan = bfunc.findKnotInterval(nu,self.full_p,urank_full[i],self.full_U)
            idxU = uspan + idx - self.full_p
            nbas = bfunc.basisFunction(uspan,urank_full[i],mu,self.full_p,self.full_U)

            nbas = np.reshape(nbas,(1,len(nbas)))

            ratFunc = (nbas*self.full_w[idxU,:].T)/(nbas@self.full_w[idxU,:])

            self.cpts[i,:] = ratFunc@self.full_P[idxU,:]
        # End for loop
        return self.cpts

    def plot_multipatch_basis_functions(self):
        """
        Put legend outside the plot
        https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot/4700762#4700762
        """
        numpoints = 81

        patch_colors = ['r', 'b']

        urank_full = np.linspace(self.full_U.min(),self.full_U.max(),numpoints)
        id_joint = np.max(np.where(urank_full<self.knot_to_repeat))
        mu_full = len(self.full_U) - 1
        nu_full = mu_full - self.full_p - 1

        N1_bas = np.zeros((nu_full+1,numpoints))

        for i in range(len(urank_full)):
            for j in range(0, nu_full+1):
                Npi = bfunc.oneBasisFunction(self.full_p, self.full_U, j, urank_full[i])
                N1_bas[j,i] = Npi

        fig, ax = plt.subplots()
        for j in range(0, nu_full+1):
            ax.plot(urank_full, N1_bas[j,:], c=patch_colors[0])
        plt.show()
    
    def plot_multipatch_curve(self):
        """
        Plot the curve
        """
        # print(self.cpts)
        fig,ax = plt.subplots()
        plt.plot(self.cpts[:,0],self.cpts[:,1])
        ax.set_aspect('equal','box')
        plt.plot(self.full_P[:,0],self.full_P[:,1],'ro')
        plt.plot(self.full_P[:,0],self.full_P[:,1])
        plt.show()

# if __name__ == '__main__':
    # main()
