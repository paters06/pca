# import numpy as np
import matplotlib.pyplot as plt

####################################################
######################PLOTS#########################
####################################################

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

def plotCurve2d(cx,cy,P,*argv):
    fig = plt.figure()
    plt.plot(cx,cy)
    plt.plot(P[0],P[1],'ro')
    plt.plot(P[0],P[1])
    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

def plotInterpolatedCurve(cx,cy,P,Q):
    fig = plt.figure()
    plt.plot(cx,cy)
    plt.plot(P[0,:],P[1,:],'ro')
    plt.plot(P[0,:],P[1,:])
    plt.plot(Q[0,:],Q[1,:],'ko')
    plt.show()

def plotTangentCurve2d(cx,cy,cpx,cpy,P,*argv):
    fig = plt.figure()
    plt.plot(P[0],P[1],'ro')
    plt.plot(P[0],P[1])
    plt.quiver(cx,cy,cpx,cpy,color=['k'])
    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()

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

def plottingSurface(cx,cy,cz,*argv):
    fig = plt.figure()
    ax = plt.axes(projection = '3d')
    ax.contour3D(cx, cy, cz, 50, cmap = 'binary')
    if len(argv)==3:
        ax.plot_wireframe(argv[0], argv[1], argv[2], color = 'red')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z');
    plt.show()

def plotTangentSurface(cx,cy,cz,cpx,cpy,cpz,px,py,pz,*argv):
    fig = plt.figure()
    ax = plt.axes(projection = '3d')
    #ax.contour3D(cx, cy, cz, 50, cmap = 'binary')
    ax.plot_wireframe(px,py,pz, color = 'red')
    plt.quiver(cx,cy,cz,cpx,cpy,cpz,color=['k'],length = 1.0,normalize = True)

    if argv != ():
        ax.set_title(argv[0])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z');
    plt.show()

def plotNodeInsertion(cx,cy,cxnew,cynew,P,Pnew,*argv):
    fig,ax = plt.subplots(1,2,sharey=True,figsize=(8,4.8))
    fig.suptitle('Knot insertion')
    ax[0].set_title('Original')
    ax[0].plot(cx,cy)
    ax[0].plot(P[0],P[1],'ro')
    ax[0].plot(P[0],P[1],'k')
    ax[1].set_title('After insertion')
    ax[1].plot(cxnew,cynew)
    ax[1].plot(Pnew[0],Pnew[1],'ro')
    ax[1].plot(Pnew[0],Pnew[1],'k')

    if argv != ():
        if argv[0] == 'yes':
            plt.savefig(argv[1]+'.png')
        else:
            plt.show()
    else:
        plt.show()
