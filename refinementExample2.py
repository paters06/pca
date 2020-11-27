import numpy as np
import pca_01 as pca
import nurbs as rbs
import plottingScripts as plts
import matplotlib.pyplot as plt
import curveRefinements as crfn

############# SINGLE KNOT INSERTION EXAMPLE #################

def singleKnotInsertionScript(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)
    Uh,Pwh = crfn.knotInsertion(U,p,Pw,0.5)
    Ph,wh = rbs.geometricControlPoints(Pwh)

    cx,cy = rbs.nurbsCurve(U,p,P,w)
    cxh,cyh = rbs.nurbsCurve(Uh,p,Ph,wh)

    # plts.plotCurve2d(cx,cy,P)
    # plts.plotCurve2d(cxh,cyh,Ph)
    plts.plotKnotInsertion(cx,cy,cxh,cyh,P,Ph)

################ KNOT REFINEMENT EXAMPLE ####################

def knotRefinementScript(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)

    Ured = U[p:-p]
    X = 0.5*(Ured[0:-1] + Ured[1:])

    Qw,Ubar = crfn.knotRefinement(U,X,p,Pw)
    Pref,wref = rbs.geometricControlPoints(Qw)

    cx,cy = rbs.nurbsCurve(U,p,P,w)
    cxref,cyref = rbs.nurbsCurve(Ubar,p,Pref,wref)

    # plts.plotCurve2d(cx,cy,P)
    # plts.plotCurve2d(cxref,cyref,Pref)
    plts.plotKnotRefinement(cx,cy,cxref,cyref,P,Pref)

################ Pre-SPLINE DECOMPOSITION EXAMPLE ####################

def preSplineDecompositionScript(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)
    cx,cy = rbs.nurbsCurve(U,p,P,w)
    # plts.plotCurve2d(cx,cy,P)

    Qdecmp,Udecmp = crfn.preSplineDecomposition(U,p,Pw)
    Pdecmp,wdecmp = rbs.geometricControlPoints(Qdecmp)

    cxdecmp,cydecmp = rbs.nurbsCurve(Udecmp,p,Pdecmp,wdecmp)
    # plts.plotCurve2d(cxdecmp,cydecmp,Pdecmp)

################ SPLINE DECOMPOSITION EXAMPLE ####################

def splineDecompositionScript(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)
    Qdecmp,Usplit = crfn.preSplineDecomposition(U,p,Pw)
    Qsplit = crfn.splineSplitting(U,p,Pw)
    Psplit,wsplit = rbs.geometricControlPoints(Qsplit)
    crfn.splineSplittingv2(U,p,Pw)

    cx,cy = rbs.nurbsCurve(U,p,P,w)
    cxs,cys = rbs.nurbsCurve(Usplit,p,P,w)

################ DEGREE ELEVATION EXAMPLE ####################

def degreeElevationScript(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)
    t = 1
    Qe,Ue,pe = crfn.degreeElevation(U,p,Pw,t)
    Pe,we = rbs.geometricControlPoints(Qe)

    cx,cy = rbs.nurbsCurve(U,p,P,w)
    cxe,cye = rbs.nurbsCurve(Ue,pe,Pe,we)

    # plts.plotCurve2d(cx,cy,P)
    # plts.plotCurve2d(cxe,cye,Pe)
    plts.plotDegreeElevation(cx,cy,cxe,cye,P,Pe)

def numberToScript(U,p,P,w,argument):
    switcher = {
    1:singleKnotInsertionScript,
    2:knotRefinementScript,
    3:preSplineDecompositionScript,
    4:splineDecompositionScript,
    5:degreeElevationScript
    }
    func = switcher.get(argument,lambda:"Invalid option")
    func(U,p,P,w)

################ GENERAL REFINEMENT EXAMPLE ####################

#Knot refinement
def hRefinement(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)

    Ured = U[p:-p]
    X = 0.5*(Ured[0:-1] + Ured[1:])

    Qwh,Uh = crfn.knotRefinement(U,X,p,Pw)
    Ph,wh = rbs.geometricControlPoints(Qwh)

    ph = p

    return Uh,ph,Ph,wh

#Degree Elevation
def pRefinement(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)
    t = 1
    Qp,Up,pp = crfn.degreeElevation(U,p,Pw,t)
    Pp,wp = rbs.geometricControlPoints(Qp)

    return Up,pp,Pp,wp

def kRefinement(U,p,P,w):
    Up,pp,Pp,wp = pRefinement(U,p,P,w)
    Uk,pk,Pk,wk = hRefinement(Up,pp,Pp,wp)

    return Uk,pk,Pk,wk

def curveRefinement(refinementlist,U,p,P,w):
    href = 0
    pref = 0
    kref = 0

    Uin = U
    pin = p
    Pin = P
    win = w

    for rfnlist in refinementlist:
        if rfnlist == 'h':
            Uout,pout,Pout,wout = hRefinement(Uin,pin,Pin,win)
            href += 1
        elif rfnlist == 'p':
            Uout,pout,Pout,wout = pRefinement(Uin,pin,Pin,win)
            pref += 1
        elif rfnlist == 'k':
            Uout,pout,Pout,wout = kRefinement(Uin,pin,Pin,win)
            kref += 1
        else:
            print("Invalid option")

        Uin = Uout
        pin = pout
        Pin = Pout
        win = wout

    print('Number of h-refinements')
    print(href)
    print('Number of p-refinements')
    print(pref)
    print('Number of k-refinements')
    print(kref)

    return Uout,pout,Pout,wout

####################################################
###################MAIN PROGRAM#####################
####################################################

Pinp = np.array([[1,0],[1,1],[0,1]])
winp = np.array([[1],[1/np.sqrt(2)],[1]])
pinp = 2
Uinp = np.array([0,0,0,1,1,1])
# Pwinp = rbs.weightedControlPoints(P,w)

# option = 1
# numberToScript(Uinp,pinp,Pinp,winp,option)

# reflist = ['k']
# reflist = ['k','h']
reflist = ['h','h','h']
Uref,pref,Pref,wref = curveRefinement(reflist,Uinp,pinp,Pinp,winp)
# print(Uref)
# print(pref)

cx,cy = rbs.nurbsCurve(Uinp,pinp,Pinp,winp)
cxref,cyref = rbs.nurbsCurve(Uref,pref,Pref,wref)

# plts.plotCurve2d(cx,cy,Pinp)
# plts.plotCurve2d(cxref,cyref,Pref)
plts.plotCurveRefinement(cx,cy,cxref,cyref,Pinp,Pref)
