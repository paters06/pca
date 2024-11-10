import matplotlib.pyplot as plt

from ..spline_functions.nurbs_curve import NURBSCurve
from ..spline_functions.nurbs_surface import NURBSSurface

class NURBSPlotter:
    def __init__(self) -> None:
        pass

    @property
    def nurbs_curve(self) -> NURBSCurve:
        return self._nurbs_curve
    
    @nurbs_curve.setter
    def nurbs_curve(self, nurbs_object: NURBSCurve) -> None:
        self._nurbs_curve = nurbs_object

    @property
    def nurbs_surface(self) -> NURBSSurface:
        return self._nurbs_surface
    
    @nurbs_surface.setter
    def nurbs_surface(self, nurbs_object: NURBSSurface) -> None:
        self._nurbs_surface = nurbs_object

    def plot_curve(self):
        """
        Plot the curve
        """
        P = self.nurbs_curve.control_points
        cpts = self.nurbs_curve.create_curve()
        fig,ax = plt.subplots()
        plt.plot(cpts[:,0],cpts[:,1])
        ax.set_aspect('equal','box')
        plt.plot(P[:,0],P[:,1],'ro')
        plt.plot(P[:,0],P[:,1])
        plt.show()

    def plot_tangent_curve(self):
        """
        Plot the tangent curve
        """
        # fig = plt.figure()
        P = self.nurbs_curve.control_points
        cpts = self.nurbs_curve.create_curve()
        cppts = self.nurbs_curve.create_tangent_curve()
        plt.plot(P[:,0],P[:,1],'ro')
        plt.plot(P[:,0],P[:,1])
        plt.plot(cpts[:,0],cpts[:,1])
        plt.quiver(cpts[:,0],cpts[:,1],cppts[:,0],cppts[:,1],color=['k'])
        plt.show()
    
    def plot_surface(self):
        """
        Plot the surface
        """
        cpts = self.nurbs_surface.create_surface()
        cx = cpts[0,:,:]
        cy = cpts[1,:,:]
        cz = cpts[2,:,:]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # ax = plt.axes(projection = '3d')
        # ax.contour3D(cx,cy,cz,cmap='viridis')
        ax.plot_surface(cx,cy,cz,cmap='viridis')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()

    def plot_tangent_surface(self, component: str):
        """
        Plot the tangent surface
        """
        cpts = self.nurbs_surface.create_surface()
        cx = cpts[0,:,:]
        cy = cpts[1,:,:]
        cz = cpts[2,:,:]

        cpu, cpv = self.nurbs_surface.create_tangent_surface()

        if component == "u":
            cpx = cpu[0,:,:]
            cpy = cpu[1,:,:]
            cpz = cpu[2,:,:]
        else:
            cpx = cpv[0,:,:]
            cpy = cpv[1,:,:]
            cpz = cpv[2,:,:]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # ax = plt.axes(projection = '3d')
        # ax.contour3D(cx,cy,cz,cmap='viridis')
        ax.plot_surface(cx,cy,cz,cmap='viridis')
        plt.quiver(cx,cy,cz,cpx,cpy,cpz,color=['k'],length = 0.01,normalize = True)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()