# Python libraries
import numpy as np

#######################################################################
# DO NOT REMOVE THIS SEGMENT
import os
import sys

# Get current working directory
# Command found in https://note.nkmk.me/en/python-script-file-path/
dir1 = os.path.dirname(os.path.abspath(__file__))
# Insert .. command to go to the upper directory
dir2 = dir1 + '/..'
# Setting the package directory path for the modules execution
sys.path.append(dir2)
#######################################################################

# Local project
import src.nurbs as rbs
import src.plottingScripts as plts
import src.surfaceRefinements as srfn

############# SINGLE KNOT INSERTION EXAMPLE #################
def singleKnotInsertionScript(U,V,p,q,P,w):
    unew = 0.5
    vnew = 0.5

    Pw = rbs.weightedControlPoints(P,w)
    Pwgrid = rbs.listToGridControlPoints(Pw,U,V,p,q)

    Uh,Vh,Qwgrid = srfn.knotInsertion(unew,vnew,"VDIR",U,V,p,q,Pwgrid)
    Qw = rbs.gridToListControlPoints(Qwgrid)
    Qh,wh = rbs.geometricControlPoints(Qw)

    cx,cy = rbs.nurbs2DField(U,V,p,q,P,w)
    cxh,cyh = rbs.nurbs2DField(Uh,Vh,p,q,Qh,wh)

    # plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])),P)
    # plts.plotting2DField(cxh,cyh,np.zeros((cxh.shape[0],cxh.shape[1])),Qh)
    plts.plotSurfaceKnotInsertion(cx,cy,cxh,cyh,P,Qh)

################ KNOT REFINEMENT EXAMPLE ####################
def knotRefinementScript(U,V,p,q,P,w):
    Ured = U[p:-p]
    X = 0.5*(Ured[0:-1] + Ured[1:])
    # X = np.array([0.25,0.5,0.75])
    # X = np.array([1/3,2/3])

    Pw = rbs.weightedControlPoints(P,w)
    Pwgrid = rbs.listToGridControlPoints(Pw,U,V,p,q)

    Ubar,Vbar,Qwgrid = srfn.knotRefinement(X,"VDIR",U,V,p,q,Pwgrid)
    Qw = rbs.gridToListControlPoints(Qwgrid)
    Qref,wref = rbs.geometricControlPoints(Qw)

    cx,cy = rbs.nurbs2DField(U,V,p,q,P,w)
    cxref,cyref = rbs.nurbs2DField(Ubar,Vbar,p,q,Qref,wref)

    # plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])),P)
    # plts.plotting2DField(cxref,cyref,np.zeros((cxref.shape[0],cxref.shape[1])),Qref)
    plts.plotSurfaceKnotRefinement(cx,cy,cxref,cyref,P,Qref)

################ Pre-SPLINE DECOMPOSITION EXAMPLE ####################
def preSplineDecompositionScript(U,V,p,q,P,w):

    Pw = rbs.weightedControlPoints(P,w)
    Pwgrid = rbs.listToGridControlPoints(Pw,U,V,p,q)

    Udecmp,Vdecmp,Qwdecmp = srfn.preSplineDecomposition("UDIR",U,V,p,q,Pwgrid)
    Qw = rbs.gridToListControlPoints(Qwdecmp)
    Qdecmp,wdecmp = rbs.geometricControlPoints(Qw)

    # cxdecmp,cydecmp = rbs.nurbs2DField(Udecmp,Vdecmp,p,q,Qdecmp,wdecmp)
    # plts.plotting2DField(cxdecmp,cydecmp,np.zeros((cxdecmp.shape[0],cxdecmp.shape[1])),Qdecmp)

################ SPLINE DECOMPOSITION EXAMPLE ####################

################ DEGREE ELEVATION FOR SURFACE ####################
def degreeElevationScript(U,V,p,q,P,w):

    Pw = rbs.weightedControlPoints(P,w)
    Pwgrid = rbs.listToGridControlPoints(Pw,U,V,p,q)

    tp = 1
    tq = 1
    Ue,Ve,pe,qe,Qwegrid = srfn.surfaceDegreeElevation("UDIR",U,V,p,q,Pwgrid,tp,tq)
    Qwe = rbs.gridToListControlPoints(Qwegrid)
    Qe,we = rbs.geometricControlPoints(Qwe)

    # print(Qwegrid)

    cx,cy = rbs.nurbs2DField(U,V,p,q,P,w)
    cxe,cye = rbs.nurbs2DField(Ue,Ve,pe,qe,Qe,we)

    # plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])),P)
    # plts.plotting2DField(cxe,cye,np.zeros((cxe.shape[0],cxe.shape[1])),Qe)
    plts.plotSurfaceDegreeElevation(cx,cy,cxe,cye,P,Qe)

def numberToScript(U,V,p,q,P,w,argument):
    switcher = {
    1:singleKnotInsertionScript,
    2:knotRefinementScript,
    3:preSplineDecompositionScript,
    4:degreeElevationScript
    }
    func = switcher.get(argument,lambda:"Invalid option")
    func(U,V,p,q,P,w)

####################################################
################## MAIN PROGRAM ####################
####################################################

Pinp = np.array([[0.5,0],[0.5,0.5],[0,0.5],[1,0],[1,1],[0,1]])
winp = np.array([[1],[0.5*np.sqrt(2)],[1],[1],[0.5*np.sqrt(2)],[1]])
pinp = 2
qinp = 1
Uinp = np.array([0,0,0,1,1,1])
Vinp = np.array([0,0,1,1])

surface1 = rbs.NURBSSurface(Uinp,Vinp,pinp,qinp,Pinp,winp)
cpts = surface1.createSurface()

# option = 4
# numberToScript(Uinp,Vinp,pinp,qinp,Pinp,winp,option)

# Uref,Vref,pref,qref,Pref,wref = srfn.hRefinement("V",Uinp,Vinp,pinp,qinp,Pinp,winp)
# Uref,Vref,pref,qref,Pref,wref = srfn.pRefinement("V",Uinp,Vinp,pinp,qinp,Pinp,winp)
# Uref,Vref,pref,qref,Pref,wref = srfn.kRefinement("U",Uinp,Vinp,pinp,qinp,Pinp,winp)

# reflist = ['h']
# dirlist = ['U']
# reflist = ['h','h']
# dirlist = ['U','V']
reflist = ['k','k','h','h']
dirlist = ['U','V','U','V']
paramlist = [[0.1,0.2,5],[0.3,0.5,3]]
srfn.defaultSurfaceRefinement(surface1,reflist,dirlist)
# srfn.customizedSurfaceRefinement(surface1,reflist,dirlist,paramlist)

cptsref = surface1.createSurface()

# plts.plottingSurface(cpts[0,:,:],cpts[1,:,:],np.zeros((cpts.shape[1],cpts.shape[2])))
# plts.plottingSurface(cptsref[0,:,:],cptsref[1,:,:],np.zeros((cptsref.shape[1],cptsref.shape[2])))
plts.plotSurfaceRefinement(cpts[0,:,:],cpts[1,:,:],cptsref[0,:,:],cptsref[1,:,:],Pinp,surface1.P)
