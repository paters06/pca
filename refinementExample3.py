import numpy as np
import pca_01 as pca
import nurbs as rbs
import plottingScripts as plts
import matplotlib.pyplot as plt
import surfaceRefinements as srfn

Pinp = np.array([[0.5,0],[0.5,0.5],[0,0.5],[1,0],[1,1],[0,1]])
winp = np.array([[1],[0.5*np.sqrt(2)],[1],[1],[0.5*np.sqrt(2)],[1]])
pinp = 2
qinp = 1
Uinp = np.array([0,0,0,1,1,1])
Vinp = np.array([0,0,1,1])

Pwinp = rbs.weightedControlPoints(Pinp,winp)
unew = 0.5
vnew = 0.5

Pwgrid = rbs.listToGridControlPoints(Pwinp,Uinp,Vinp,pinp,qinp)

############# SINGLE KNOT INSERTION EXAMPLE #################

# Uh,Vh,Qwgrid = srfn.knotInsertion(unew,vnew,"VDIR",Uinp,Vinp,pinp,qinp,Pwgrid)
# Qw = rbs.gridToListControlPoints(Qwgrid)
# Qh,wh = rbs.geometricControlPoints(Qw)

# cx,cy = rbs.nurbs2DField(Uinp,Vinp,pinp,qinp,Pinp,winp)
# cxh,cyh = rbs.nurbs2DField(Uh,Vh,pinp,qinp,Qh,wh)

# plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])),Pinp)
# plts.plotting2DField(cxh,cyh,np.zeros((cxh.shape[0],cxh.shape[1])),Qh)
# plts.plotSurfaceKnotInsertion(cx,cy,cxh,cyh,Pinp,Qh)

################ KNOT REFINEMENT EXAMPLE ####################

Ured = Uinp[pinp:-pinp]
# X = 0.5*(Ured[0:-1] + Ured[1:])
# X = np.array([0.25,0.5,0.75])
X = np.array([1/3,2/3])

Ubar,Vbar,Qwgrid = srfn.knotRefinement(X,"UDIR",Uinp,Vinp,pinp,qinp,Pwgrid)
# print(Qwgrid)
Qw = rbs.gridToListControlPoints(Qwgrid)
Qref,wref = rbs.geometricControlPoints(Qw)

cx,cy = rbs.nurbs2DField(Uinp,Vinp,pinp,qinp,Pinp,winp)
cxref,cyref = rbs.nurbs2DField(Ubar,Vbar,pinp,qinp,Qref,wref)

# plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])),Pinp)
# plts.plotting2DField(cxref,cyref,np.zeros((cxref.shape[0],cxref.shape[1])),Qref)
plts.plotSurfaceKnotRefinement(cx,cy,cxref,cyref,Pinp,Qref)
