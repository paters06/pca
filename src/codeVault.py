import numpy as np

################## NURBS CURVES ######################

def nurbsCurve(U,p,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)

    cpts = np.zeros((numpoints,2))

    mu = len(U) - 1
    nu = mu - p - 1
    idx = np.arange(0,p+1)

    for i in range(len(urank)):
        uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
        idxU = uspan + idx - p
        nbas = bfunc.basisFunction(uspan,urank[i],mu,p,U)

        nbas = np.reshape(nbas,(1,len(nbas)))

        ratFunc = (nbas*w[idxU,:].T)/(nbas@w[idxU,:])

        cpts[i,:] = ratFunc@P[idxU,:]
    return cpts

def nurbsCurveTangent(U,p,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)

    cppts = np.zeros((numpoints,2))

    Pw = weightedControlPoints(P,w)

    mu = len(U) - 1
    nu = mu - p - 1
    idx = np.arange(0,p+1)

    # Derivative order
    d = 1

    for i in range(len(urank)):
        uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
        idxU = uspan + idx - p
        nbas = bfunc.basisFunction(uspan,urank[i],mu,p,U)
        dnbasU = bfunc.derBasisFunction(uspan,urank[i],mu,p,U,d)

        """Just using the formulas way"""
       # Aders = dnbasU*Pw[idxU,-1].T
       # wders = dnbasU@Pw[idxU,-1]
       # wders = np.reshape(wders,(len(wders),1))
       # dRatdU = (Aders - (wders*nbas))/wders[0]
       # Ck = dRatdU@P[idxU,:]

        """The NURBS Book way"""

       # dCw = dnbasU@Pw[idxU,:]
       #  Selecting colums from 0 to the previous to the last one
       # Aders = dCw[:,0:-1]
       #  Selecting the last column
       # wders = dCw[:,-1]
       # Ck = rationalCurveDerivative(Aders,wders,d)

        """Hughes' way"""
        Aders = dnbasU*Pw[idxU,-1].T
        wders = dnbasU@Pw[idxU,-1]
        dRatdU = univariateRationalDerivative(Aders,wders,d)
        Ck = dRatdU@P[idxU,:]

        cppts[i,:] = Ck[d,:]

    return cppts

################# NURBS SURFACES #####################

def nurbs2DBoundary(U,V,p,q,P,w):
    Pwl = weightedControlPoints(P,w)
    Pw = listToGridControlPoints(Pwl,U,V,p,q)

    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1

    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)

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
            uspan = bfunc.findKnotInterval(nu,p,ppath[0],U)
            vspan = bfunc.findKnotInterval(nv,q,ppath[1],V)
            idR = nonZeroIndicesSurface(uspan,vspan,p,q,nu)

            R = bivariateRationalFunction(mu,mv,p,q,uspan,vspan,ppath[0],ppath[1],U,V,Pw)
            S = R@P[idR,:]
            coor[ipath,:] = S
            ipath += 1

        boundarycoor.append(coor)

    for bc in range(len(boundarycoor)):
        if bc == 0:
            boundarycoor1 = boundarycoor[bc]
        else:
            boundarycoor1 = np.vstack((boundarycoor1,boundarycoor[bc]))

    return boundarycoor1

def nurbsSurface(U,V,p,q,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    Pwl = weightedControlPoints(P,w)
    Pw = listToGridControlPoints(Pwl,U,V,p,q)

    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1

    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)

    cpts = np.zeros((P.shape[1],len(urank),len(vrank)))

    for j in range(len(vrank)):
        for i in range(len(urank)):
            uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
            vspan = bfunc.findKnotInterval(nv,q,vrank[j],V)

            idR = nonZeroIndicesSurface(uspan,vspan,p,q,nu)

            R = bivariateRationalFunction(mu,mv,p,q,uspan,vspan,urank[i],vrank[j],U,V,Pw)
            S = R@P[idR,:]

            cpts[:,i,j] = S
    return cpts

def nurbsSurfaceTangent(U,V,p,q,P,w):
    numpoints = 41
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    Pwl = weightedControlPoints(P,w)
    Pw = listToGridControlPoints(Pwl,U,V,p,q)

    mu = len(U) - 1
    mv = len(V) - 1
    nu = mu - p - 1
    nv = mv - q - 1

    idxu = np.arange(0,p+1)
    idxv = np.arange(0,q+1)

    cpu = np.zeros((P.shape[1],len(urank),len(vrank)))
    cpv = np.zeros((P.shape[1],len(urank),len(vrank)))

    for j in range(len(vrank)):
        for i in range(len(urank)):
            uspan = bfunc.findKnotInterval(nu,p,urank[i],U)
            vspan = bfunc.findKnotInterval(nv,q,vrank[j],V)

            idR = nonZeroIndicesSurface(uspan,vspan,p,q,nu)

            Ralph = bivariateRationalGradient(mu,mv,p,q,uspan,vspan,urank[i],vrank[j],U,V,Pw)
            dS = Ralph@P[idR,:]

            cpu[:,i,j] = dS[1,:]
            cpv[:,i,j] = dS[2,:]

    return cpu,cpv

################# MULTIPATCH NURBS SURFACES #####################

def multipatchNurbsSurface(mulU,mulV,mulp,mulq,fullP,fullw,idctrlpts):
    numpoints = 11

    numpatches = len(mulU)

    fullcpts = []

    for ipatch in range(0,numpatches):
        Ui = mulU[ipatch]
        Vi = mulV[ipatch]

        pi = mulp[ipatch]
        qi = mulq[ipatch]

        urank = np.linspace(Ui.min(),Ui.max(),numpoints)
        vrank = np.linspace(Vi.min(),Vi.max(),numpoints)

        Pi = fullP[idctrlpts[ipatch],:]
        wi = fullw[idctrlpts[ipatch],:]

        Pwl = weightedControlPoints(Pi,wi)
        Pw = listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

        mui = len(Ui) - 1
        mvi = len(Vi) - 1
        nui = mui - pi - 1
        nvi = mvi - qi - 1

        cpts = np.zeros((Pi.shape[1],len(urank),len(vrank)))

        for j in range(len(vrank)):
            for i in range(len(urank)):
                uspan = bfunc.findKnotInterval(nui,pi,urank[i],Ui)
                vspan = bfunc.findKnotInterval(nvi,qi,vrank[j],Vi)

                idR = nonZeroIndicesSurface(uspan,vspan,pi,qi,nui)

                R = bivariateRationalFunction(mui,mvi,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                S = R@Pi[idR,:]

                cpts[:,i,j] = S
            # End i loop
        # End j loop

        fullcpts.append(cpts)
    # End patch loop

    return fullcpts

def multipatchNurbsSurfaceTangent(mulU,mulV,mulp,mulq,fullP,fullw,idctrlpts):
    numpoints = 11

    numpatches = len(mulU)

    fullcpu = []
    fullcpv = []

    for ipatch in range(0,numpatches):
        Ui = mulU[ipatch]
        Vi = mulV[ipatch]

        pi = mulp[ipatch]
        qi = mulq[ipatch]

        urank = np.linspace(Ui.min(),Ui.max(),numpoints)
        vrank = np.linspace(Vi.min(),Vi.max(),numpoints)

        Pi = fullP[idctrlpts[ipatch]]
        wi = fullw[idctrlpts[ipatch]]

        Pwl = weightedControlPoints(Pi,wi)
        Pw = listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

        mui = len(Ui) - 1
        mvi = len(Vi) - 1
        nui = mui - pi - 1
        nvi = mvi - qi - 1

        # cpts = np.zeros((Pi.shape[1],len(urank),len(vrank)))
        cpu = np.zeros((Pi.shape[1],len(urank),len(vrank)))
        cpv = np.zeros((Pi.shape[1],len(urank),len(vrank)))

        for j in range(len(vrank)):
            for i in range(len(urank)):
                uspan = bfunc.findKnotInterval(nui,pi,urank[i],Ui)
                vspan = bfunc.findKnotInterval(nvi,qi,vrank[j],Vi)

                idR = nonZeroIndicesSurface(uspan,vspan,pi,qi,nui)

                R = bivariateRationalFunction(mui,mvi,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                S = R@Pi[idR,:]

                Ralph = bivariateRationalGradient(mui,mvi,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pw)
                dS = Ralph@Pi[idR,:]

                # cpts[:,i,j] = S
                cpu[:,i,j] = dS[1,:]
                cpv[:,i,j] = dS[2,:]

        fullcpu.append(cpu)
        fullcpv.append(cpv)

    return fullcpu,fullcpv

##################### PLOTS ########################

def plotCurve2d(cpts,P,*argv):
   fig,ax = plt.subplots()
   plt.plot(cpts[:,0],cpts[:,1])
   ax.set_aspect('equal','box')
   plt.plot(P[:,0],P[:,1],'ro')
   plt.plot(P[:,0],P[:,1])
   if argv != ():
       if argv[0] == 'yes':
           plt.savefig(argv[1]+'.png')
       else:
           plt.show()
   else:
       plt.show()

def plotTangentCurve2d(cpts,cppts,P,*argv):
   fig = plt.figure()
   plt.plot(P[:,0],P[:,1],'ro')
   plt.plot(P[:,0],P[:,1])
   plt.plot(cpts[:,0],cpts[:,1])
   plt.quiver(cpts[:,0],cpts[:,1],cppts[:,0],cppts[:,1],color=['k'])
   if argv != ():
       if argv[0] == 'yes':
           plt.savefig(argv[1]+'.png')
       else:
           plt.show()
   else:
       plt.show()

def plottingSurface(cx,cy,cz,*argv):
   fig = plt.figure()
   ax = plt.axes(projection = '3d')
   # ax.contour3D(cx, cy, cz, 50, cmap = 'viridis')
   ax.plot_surface(cx, cy, cz, cmap = 'viridis')
   if len(argv)==3:
       px = np.reshape(argv[0],(len(argv[0]),1))
       py = np.reshape(argv[1],(len(argv[1]),1))
       pz = np.reshape(argv[2],(len(argv[2]),1))
       ax.plot_wireframe(px,py,pz, color = 'red')
   ax.set_xlabel('x')
   ax.set_ylabel('y')
   ax.set_zlabel('z')
   plt.show()

def plotTangentSurface(cx,cy,cz,cpx,cpy,cpz,*argv):
   fig = plt.figure()
   ax = plt.axes(projection = '3d')
   #ax.contour3D(cx, cy, cz, 50, cmap = 'binary')
   ax.plot_surface(cx, cy, cz, cmap = 'viridis')
   plt.quiver(cx,cy,cz,cpx,cpy,cpz,color=['k'],length = 0.5,normalize = True)
   # plt.quiver(cx,cy,cz,cpx,cpy,cpz,color=['k'],normalize = True)
   if len(argv)==3:
       px = np.reshape(argv[0],(len(argv[0]),1))
       py = np.reshape(argv[1],(len(argv[1]),1))
       pz = np.reshape(argv[2],(len(argv[2]),1))
       ax.plot_wireframe(px,py,pz, color = 'red')

   if argv != ():
       ax.set_title(argv[0])
   ax.set_xlabel('x')
   ax.set_ylabel('y')
   ax.set_zlabel('z');
   plt.show()
