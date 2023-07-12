import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

####################################################
##################### PLOTS ########################
####################################################

##################### IGA 1D ########################

def plotCurve1d(cx,uy,*argv):
    fig = plt.figure()
    plt.plot(cx,uy)
    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotComparisonCurves(cx,cy,cey,field):
    fig = plt.figure()
    plt.plot(cx,cy)
    plt.plot(cx,cey)

    plt.xlabel('Distance [m]')
    if field == 'u':
        plt.ylabel('Displacement [m]')
        plt.title('Numerical displacement in x-direction')
    elif field == 's':
        plt.ylabel('Stress [Pa]')
        plt.title('Numerical axial stress in x-direction')
    else:
        print('Wrong field selected')

    plt.legend(['IGA','Exact'])
    plt.show()

def plot1DField(cx,uy,field,*argv):
    fig = plt.figure()
    plt.plot(cx,uy)

    plt.xlabel('Distance [m]')
    if field == 'u':
        plt.ylabel('Displacement [m]')
        plt.title('Numerical displacement in x-direction')
    elif field == 's':
        plt.ylabel('Stress [Pa]')
        plt.title('Numerical axial stress in x-direction')
    else:
        print('Wrong field selected')

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

##################### INTERPOLATION ########################

def plotInterpolatedCurve(cx,cy,P,Q):
    fig = plt.figure()
    plt.plot(cx,cy)
    plt.plot(P[0,:],P[1,:],'ro')
    plt.plot(P[0,:],P[1,:])
    plt.plot(Q[0,:],Q[1,:],'ko')
    plt.show()

##################### 3D-CURVES ########################

def plotting3d(cx,cy,cz,*argv):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(cx, cy, cz, 'gray')
    if argv != ():
        if argv[0]=='on':
            if len(argv)==4:
                ax.plot3D(argv[1],argv[2],argv[3], 'red')
            elif len(argv)> 4:
                sys.exit("Too much arguments, please delete one or more")
            else:
                sys.exit("Missing arguments to plot control points")

##################### 2D-FIELDS ########################

def plotting2DField(cx,cy,fz,*argv):
    fig = plt.figure()
    ax = plt.axes()
    titlestring = ""
    colorbarstring = "value"
    # field = ax.pcolormesh(cx,cy,fz)
    field = ax.pcolormesh(cx,cy,fz,vmin=fz.min(),vmax=fz.max())
    # field = ax.pcolormesh(cx,cy,fz,shading='gouraud')
    # field = ax.pcolormesh(cx,cy,fz,shading='gouraud',vmin=fz.min(),vmax=fz.max())
    if argv != ():
        stringlegends = argv[0]
        titlestring = stringlegends[0]
        colorbarstring = stringlegends[1]

    cb = fig.colorbar(field)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(titlestring)
    # cb.set_label(colorbarstring)
    plt.show()

def plotMultipatchField(fullc,fullf,comp,*argv):
    fig = plt.figure()
    ax = plt.axes()
    titlestring = ""
    colorbarstring = "value"

    cmap = 'viridis'
    # cmap = 'RdBu_r'
    num_levels = 20
    field = ax.tricontour(fullc[:,0], fullc[:,1], fullf[:,comp], levels=num_levels, linewidths=0.5, colors='k')
    cntr2 = ax.tricontourf(fullc[:,0], fullc[:,1], fullf[:,comp], levels=num_levels, cmap=cmap)
    # ax.plot(fullc[:,0], fullc[:,1], 'ko', ms=3)
    
    if argv != ():
        stringlegends = argv[0]
        titlestring = stringlegends[0]
        colorbarstring = stringlegends[1]

    cb = fig.colorbar(cntr2)
    # cb = fig.colorbar(field,extend='both')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(titlestring)
    # cb.set_label(colorbarstring)
    plt.show()

##################### SURFACE REFINEMENT ########################

def plotKnotInsertion(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Knot insertion')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].plot(cx,cy)
    ax[0].plot(P[:,0],P[:,1],'ro')
    ax[0].plot(P[:,0],P[:,1],'k')

    ax[1].set_title('After insertion')
    ax[1].set_aspect('equal','box')
    ax[1].plot(cxnew,cynew)
    ax[1].plot(Pnew[:,0],Pnew[:,1],'ro')
    ax[1].plot(Pnew[:,0],Pnew[:,1],'k')

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotKnotRefinement(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Knot refinement')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].plot(cx,cy)
    ax[0].plot(P[:,0],P[:,1],'ro')
    ax[0].plot(P[:,0],P[:,1],'k')

    ax[1].set_title('After refinement')
    ax[1].set_aspect('equal','box')
    ax[1].plot(cxnew,cynew)
    ax[1].plot(Pnew[:,0],Pnew[:,1],'ro')
    ax[1].plot(Pnew[:,0],Pnew[:,1],'k')

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotDegreeElevation(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Degree Elevation')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].plot(cx,cy)
    ax[0].plot(P[:,0],P[:,1],'ro')
    ax[0].plot(P[:,0],P[:,1],'k')

    ax[1].set_title('After elevation')
    ax[1].set_aspect('equal','box')
    ax[1].plot(cxnew,cynew)
    ax[1].plot(Pnew[:,0],Pnew[:,1],'ro')
    ax[1].plot(Pnew[:,0],Pnew[:,1],'k')

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotCurveRefinement(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Curve refinement')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].plot(cx,cy)
    ax[0].plot(P[:,0],P[:,1],'ro')
    ax[0].plot(P[:,0],P[:,1],'k')

    ax[1].set_title('After refinement')
    ax[1].set_aspect('equal','box')
    ax[1].plot(cxnew,cynew)
    ax[1].plot(Pnew[:,0],Pnew[:,1],'ro')
    ax[1].plot(Pnew[:,0],Pnew[:,1],'k')

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotSurfaceKnotInsertion(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Knot insertion')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].pcolormesh(cx,cy,np.zeros((cx.shape[0],cx.shape[1])))
    ax[0].scatter(P[:,0],P[:,1])

    ax[1].set_title('After insertion')
    ax[1].set_aspect('equal','box')
    ax[1].pcolormesh(cxnew,cynew,np.zeros((cxnew.shape[0],cxnew.shape[1])))
    ax[1].scatter(Pnew[:,0],Pnew[:,1])

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotSurfaceKnotRefinement(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Knot refinement')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].pcolormesh(cx,cy,np.zeros((cx.shape[0],cx.shape[1])))
    ax[0].scatter(P[:,0],P[:,1])

    ax[1].set_title('After refinement')
    ax[1].set_aspect('equal','box')
    ax[1].pcolormesh(cxnew,cynew,np.zeros((cxnew.shape[0],cxnew.shape[1])))
    ax[1].scatter(Pnew[:,0],Pnew[:,1])

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotSurfaceDegreeElevation(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Degree Elevation')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].pcolormesh(cx,cy,np.zeros((cx.shape[0],cx.shape[1])))
    ax[0].scatter(P[:,0],P[:,1])

    ax[1].set_title('After elevation')
    ax[1].set_aspect('equal','box')
    ax[1].pcolormesh(cxnew,cynew,np.zeros((cxnew.shape[0],cxnew.shape[1])))
    ax[1].scatter(Pnew[:,0],Pnew[:,1])

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotSurfaceRefinement(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Surface refinement')

    ax[0].set_title('Original')
    ax[0].set_aspect('equal','box')
    ax[0].pcolormesh(cx,cy,np.zeros((cx.shape[0],cx.shape[1])))
    ax[0].scatter(P[:,0],P[:,1])

    ax[1].set_title('After refinement')
    ax[1].set_aspect('equal','box')
    ax[1].pcolormesh(cxnew,cynew,np.zeros((cxnew.shape[0],cxnew.shape[1])))
    ax[1].scatter(Pnew[:,0],Pnew[:,1])

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()
