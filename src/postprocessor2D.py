# Python libraries
import numpy as np
import matplotlib.pyplot as plt

# Local project
import src.nurbs as rbs
import src.plottingScripts as plts
import src.linearElastoStaticsSolver as linElastStat

################ POSTPROCESSING ####################

def displacementField(numpoints,U,V,p,q,P,D,w,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    dx = np.reshape(D[:,0],(D.shape[0],1))
    dy = np.reshape(D[:,1],(D.shape[0],1))

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

    # Geometric coordinates
    cx = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    cy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    # Displacements
    ux = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    uy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):
                ratFunc = rbs.ratFunction(U,V,w,p,q,urank[i],vrank[j])
                # ratFunc,dn2du,dn2dv = rbs.rationalFunctionAndGradient(U,V,w,p,q,urank[i],vrank[j])

                ux[iu*numpoints + i,jv*numpoints + j] = ratFunc@dx
                uy[iu*numpoints + i,jv*numpoints + j] = ratFunc@dy

                cx[iu*numpoints + i,jv*numpoints + j] = ratFunc@px
                cy[iu*numpoints + i,jv*numpoints + j] = ratFunc@py

    return ux,uy,cx,cy

# Improve this function
def stressField(numpoints,U,V,p,q,P,w,dtot,dmat,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    paramgrad = np.zeros((2,2))

    sx = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    sy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    sxy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    maxSx = -1e20
    maxu = 0
    maxv = 0
    maxelem = 0

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):

                xpcoor = urank[i]
                ypcoor = vrank[j]

                jac = linElastStat.jacobian(U,V,w,p,q,xpcoor,ypcoor,P[:,0],P[:,1])
                detJac = np.linalg.det(jac)

                if abs(detJac) > 1e-5:
                    bmat = linElastStat.strainDisplacementMatrix(U,V,w,p,q,xpcoor,ypcoor,jac)
                    svec = dmat@(bmat@dtot)
                else:
                    print("Singularity")
                    print(detJac)

                    xleft = urank[i] - 0.05
                    xright = urank[i] + 0.05
                    ydown = vrank[j] - 0.05
                    yup = vrank[j] + 0.05

                    if xleft < 1e-5:
                        xleft = 0.0

                    if xright > (1.0 + 1e-5):
                        xright = 1.0

                    if ydown < 1e-5:
                        ydown = 0.0

                    if yup > (1.0 + 1e-5):
                        yup = 1.0

                    apt = np.array([xleft,yup])
                    bpt = np.array([xright,yup])
                    cpt = np.array([xleft,ydown])
                    dpt = np.array([xright,ydown])

                    jaca = linElastStat.jacobian(U,V,w,p,q,apt[0],apt[1],P[:,0],P[:,1])
                    bmata = linElastStat.strainDisplacementMatrix(U,V,w,p,q,apt[0],apt[1],jaca)
                    sveca = dmat@(bmata@dtot)

                    jacb = linElastStat.jacobian(U,V,w,p,q,bpt[0],bpt[1],P[:,0],P[:,1])
                    bmatb = linElastStat.strainDisplacementMatrix(U,V,w,p,q,bpt[0],bpt[1],jacb)
                    svecb = dmat@(bmatb@dtot)

                    jacc = linElastStat.jacobian(U,V,w,p,q,cpt[0],cpt[1],P[:,0],P[:,1])
                    bmatc = linElastStat.strainDisplacementMatrix(U,V,w,p,q,cpt[0],cpt[1],jacc)
                    svecc = dmat@(bmatc@dtot)

                    jacd = linElastStat.jacobian(U,V,w,p,q,dpt[0],dpt[1],P[:,0],P[:,1])
                    bmatd = linElastStat.strainDisplacementMatrix(U,V,w,p,q,dpt[0],dpt[1],jacd)
                    svecd = dmat@(bmatd@dtot)

                    svec = 0.25*(sveca + svecb + svecc + svecd)

                sx[iu*numpoints + i,jv*numpoints + j] = svec[0]
                sy[iu*numpoints + i,jv*numpoints + j] = svec[1]
                sxy[iu*numpoints + i,jv*numpoints + j] = svec[2]

    # Maximum stress value
    # maxSx = np.amax(sx)
    # print('Maximum stress value',maxSx)
    # # Index of the minimum stress value
    # imaxSx = np.where(sx == np.amax(sx))
    # listOfCordinates = list(zip(imaxSx[0], imaxSx[1]))
    # for cord in listOfCordinates:
    #     print(cord)

    return sx,sy,sxy

def plotDisplacementFields(cx,cy,ux,uy):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Uncomment for cantileverBeam.py
    # fig, axs = plt.subplots(2,1,sharex='col',sharey='row')

    # Uncomment for pressureCylinder.py and plateWithHole.py
    fig, axs = plt.subplots(1,2,sharex='col',sharey='row')

    fig.suptitle('Displacement field components')
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    field1 = axs[0].pcolormesh(cx,cy,ux,vmin=ux.min(),vmax=ux.max())
    axs[0].set_title('Ux')
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('y')
    axs[0].set_aspect('equal')
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right",size="5%",pad=0.1)
    cb1 = fig.colorbar(field1,cax=cax,label='[m]')

    field2 = axs[1].pcolormesh(cx,cy,uy,vmin=uy.min(),vmax=uy.max())
    axs[1].set_title('Uy')
    axs[1].set_xlabel('x')
    axs[1].set_ylabel('y')
    axs[1].set_aspect('equal')
    divider = make_axes_locatable(axs[1])
    cax = divider.append_axes("right",size="5%",pad=0.1)
    cb2 = fig.colorbar(field2,cax=cax,label='[m]')

    plt.show()

def plotStressFields(cx,cy,sx,sy,sxy,svm):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Uncomment for cantileverBeam.py
    # fig, axs = plt.subplots(4,1,sharex='col',sharey='row')

    # Uncomment for pressureCylinder.py and plateWithHole.py
    fig, axs = plt.subplots(2,2,sharex='col',sharey='row')

    fig.suptitle('Stress field components')
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    field1 = axs[0,0].pcolormesh(cx,cy,sx,vmin=sx.min(),vmax=sx.max())
    axs[0,0].set_title('Sx')
    axs[0,0].set_xlabel('x')
    axs[0,0].set_ylabel('y')
    axs[0,0].set_aspect('equal')
    divider = make_axes_locatable(axs[0,0])
    cax1 = divider.append_axes("right",size="5%",pad=0.1)
    cb1 = fig.colorbar(field1,cax=cax1,label='[Pa]')

    field2 = axs[0,1].pcolormesh(cx,cy,sy,vmin=sy.min(),vmax=sy.max())
    axs[0,1].set_title('Sy')
    axs[0,1].set_xlabel('x')
    axs[0,1].set_ylabel('y')
    axs[0,1].set_aspect('equal')
    divider = make_axes_locatable(axs[0,1])
    cax2 = divider.append_axes("right",size="5%",pad=0.1)
    cb2 = fig.colorbar(field2,cax=cax2,label='[Pa]')

    field3 = axs[1,0].pcolormesh(cx,cy,sxy,vmin=sxy.min(),vmax=sxy.max())
    axs[1,0].set_title('Sxy')
    axs[1,0].set_xlabel('x')
    axs[1,0].set_ylabel('y')
    axs[1,0].set_aspect('equal')
    divider = make_axes_locatable(axs[1,0])
    cax3 = divider.append_axes("right",size="5%",pad=0.1)
    cb3 = fig.colorbar(field3,cax=cax3,label='[Pa]')

    field4 = axs[1,1].pcolormesh(cx,cy,svm,vmin=svm.min(),vmax=svm.max())
    axs[1,1].set_title('Von Mises stress')
    axs[1,1].set_xlabel('x')
    axs[1,1].set_ylabel('y')
    axs[1,1].set_aspect('equal')
    divider = make_axes_locatable(axs[1,1])
    cax4 = divider.append_axes("right",size="5%",pad=0.1)
    cb4 = fig.colorbar(field4,cax=cax4,label='[Pa]')

    plt.show()

def postProcessing(U,V,p,q,P,D,w,paramnodes,nodeselem,dtot,dmat):
    numelems = nodeselem.shape[0]

    if numelems < 5:
        numpoints = 11
    elif numelems >= 5 and numelems < 10:
        numpoints = 9
    elif numelems >= 10 and numelems < 20:
        numpoints = 7
    else:
        numpoints = 5

    ux,uy,cx,cy = displacementField(numpoints,U,V,p,q,P,D,w,paramnodes,nodeselem)
    sx,sy,sxy = stressField(numpoints,U,V,p,q,P,w,dtot,dmat,paramnodes,nodeselem)
    svm = np.sqrt( sx**2 - 2*sx*sy + sy**2 + 3*sxy**2 )
    # plts.plotting2DField(cx,cy,ux,P,["Ux Displacement Field","[m]"])
    # plts.plotting2DField(cx,cy,sx,P,["Sx Stress Field","[Pa]"])
    # plotDisplacementFields(cx,cy,ux,uy)
    plotStressFields(cx,cy,sx,sy,sxy,svm)
