# Python libraries
import numpy as np

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs

tol = 1e-5

################ KNOT REFINEMENT FOR SURFACE ####################

def knotRefinement(dir,X,U,V,p,q,Pwg):
    if dir == "UDIR":
        a = bfunc.findKnotInterval(len(U)-p-1,p,X[0],U)
        b = bfunc.findKnotInterval(len(U)-p-1,p,X[-1],U)
        b += 1

        n = Pwg.shape[1] - 1
        r = len(X) - 1
        m = n + p + 1
        mcol = Pwg.shape[2] - 1
        Qwg = np.zeros((Pwg.shape[0],n + 1 + r + 1,Pwg.shape[2]))
        Ubar = np.zeros(len(U) + r + 1)

        #Initializing Ubar
        for j in range(0,a+1):
            Ubar[j] = U[j]

        for j in range(b+p,m+1):
            Ubar[j+r+1] = U[j]

        #Copying V into Vbar
        Vbar = V

        #Save unaltered control points
        for row in range(0,mcol+1):

            for k in range(0,a-p+1):
                Qwg[:,k,row] = Pwg[:,k,row]

            for k in range(b-1,n+1):
                Qwg[:,k+r+1,row] = Pwg[:,k,row]

        i = b + p - 1
        k = b + p + r
        for j in range(r,-1,-1):

            while X[j] <= U[i] and i > a:
                Ubar[k] = U[i]
                for row in range(0,mcol+1):
                    Qwg[:,k-p-1,row] = Pwg[:,i-p-1,row]

                k -= 1
                i -= 1

            for row in range(0,mcol+1):
                Qwg[:,k-p-1,row] = Qwg[:,k-p,row]

            for l in range(1,p+1):
                ind = k - p + l
                alpha = Ubar[k+l] - X[j]
                if abs(alpha) < 1e-5:
                    for row in range(0,mcol+1):
                        Qwg[:,ind-1,row] = Qwg[:,ind,row]
                else:
                    alpha /= (Ubar[k+l] - U[i-p+l])
                    for row in range(0,mcol+1):
                        Qwg[:,ind-1,row] = alpha*Qwg[:,ind-1,row] + (1.0 - alpha)*Qwg[:,ind,row]

            Ubar[k] = X[j]
            k -= 1

    if dir == "VDIR":
        a = bfunc.findKnotInterval(len(V)-q-1,q,X[0],V)
        b = bfunc.findKnotInterval(len(V)-q-1,q,X[-1],V)
        b += 1

        n = Pwg.shape[2] - 1
        r = len(X) - 1
        m = n + q + 1
        nrow = Pwg.shape[1] - 1
        Qwg = np.zeros((Pwg.shape[0],Pwg.shape[1],n + 1 + r + 1))
        Vbar = np.zeros(len(V) + r + 1)

        #Initializing Ubar
        for j in range(0,a+1):
            Vbar[j] = V[j]

        for j in range(b+q,m+1):
            Vbar[j+r+1] = V[j]

        #Copying U into Ubar
        Ubar = U

        #Save unaltered control points
        for col in range(0,nrow+1):

            for k in range(0,a-q+1):
                Qwg[:,col,k] = Pwg[:,col,k]

            for k in range(b-1,n+1):
                Qwg[:,col,k+r+1] = Pwg[:,col,k]

        i = b + q - 1
        k = b + q + r
        for j in range(r,-1,-1):

            while X[j] <= V[i] and i > a:
                Vbar[k] = V[i]
                for col in range(0,nrow+1):
                    Qwg[:,col,k-q-1] = Pwg[:,col,i-q-1]

                k -= 1
                i -= 1

            for col in range(0,nrow+1):
                Qwg[:,col,k-q-1] = Qwg[:,col,k-q]

            for l in range(1,q+1):
                ind = k - q + l
                alpha = Vbar[k+l] - X[j]
                if abs(alpha) < 1e-5:
                    for col in range(0,nrow+1):
                        Qwg[:,col,ind-1] = Qwg[:,col,ind]
                else:
                    alpha /= (Vbar[k+l] - V[i-q+l])
                    for col in range(0,nrow+1):
                        Qwg[:,col,ind-1] = alpha*Qwg[:,col,ind-1] + (1.0 - alpha)*Qwg[:,col,ind]

            Vbar[k] = X[j]
            k -= 1

    return Ubar,Vbar,Qwg

################ DEGREE ELEVATION FOR SURFACE ####################

def binomialCoefficient(a,b):
    bc = 1.0
    for j in range(1,b+1):
        bc *= ((a+1-j)/j)

    return bc

def surfaceDegreeElevation(dir,U,V,p,q,Pwg,tp,tq):
    if dir == "UDIR":
        n = Pwg.shape[1] - 1 #n is the number of rows of Pwg minus 1
        mN = n + p + 1 #mN is the m value from the U knot vector

        ph = p + tp
        ph2 = ph//2

        s0 = len(np.unique(U))-2

        bezalfs = np.zeros((p+tp+1,p+1))
        bpts = np.zeros((Pwg.shape[0],p+1,Pwg.shape[2]))
        ebpts = np.zeros((Pwg.shape[0],p+tp+1,Pwg.shape[2]))
        Nextbpts = np.zeros((Pwg.shape[0],p-1,Pwg.shape[2]))
        alfs = np.zeros(p-1)
        Uh = np.zeros(len(U) + len(np.unique(U))*tp)
        Qwg = np.zeros((Pwg.shape[0],len(Uh) - (p+tp) - 1,Pwg.shape[2]))

        Vh = V
        qh = q

        #Compute bezier degree elevation coefficients
        bezalfs[0][0] = 1.0
        bezalfs[ph][p] = 1.0

        for i in range(1,ph2+1):
            inv = 1.0/binomialCoefficient(ph,i)
            mpi = min(p,i)
            for j in range(max(0,i-tp),mpi+1):
                bezalfs[i][j] = inv*binomialCoefficient(p,j)*binomialCoefficient(tp,i-j)

        for i in range(ph2+1,ph):
            mpi = min(p,i)
            for j in range(max(0,i-tp),mpi+1):
                bezalfs[i][j] = bezalfs[ph-i][p-j]

        mh = ph
        kind = ph + 1
        r = -1
        a = p
        b = p + 1
        cind = 1
        ua = Uh[0]
        # Qw[0,:] = Pw[0,:]
        Qwg[:,0,:] = Pwg[:,0,:]

        for i in range(0,ph+1):
            Uh[i] = ua

        #Initialize first bezier segment
        for i in range(0,p+1):
            # bpts[i,:] = Pw[i,:]
            bpts[:,i,:] = Pwg[:,i,:]

        #Big loop through knot vector
        while b < mN:
            i = b
            while b < mN and abs(U[b+1] - U[b]) < 1e-5:
                b += 1

            mul = b - i + 1
            mh += mul + tp
            ub = U[b]
            oldr = r
            r = p - mul
            #Insert knot u(b) r times
            if oldr > 0:
                lbz = (oldr + 2)//2
            else:
                lbz = 1

            if r > 0:
                rbz = ph - (r + 1)//2
            else:
                rbz = ph

            #Insert knot to get bezier segment
            if r > 0:
                numer = ub - ua
                for k in range(p,mul,-1):
                    alfs[k-mul-1] = numer/(U[a+k] - ua)

                for j in range(1,r+1):
                    save = r - j
                    s = mul + j
                    for k in range(p,s-1,-1):
                        bpts[:,k,:] = alfs[k-s]*bpts[:,k,:] + (1.0 - alfs[k-s])*bpts[:,k-1,:]

                    Nextbpts[:,save,:] = bpts[:,p,:]

            #End of insert knot
            #Degree elevate bezier
            #Only points lbz,...,ph are used below
            # ebpts = np.zeros((2,p+t+1))
            for i in range(lbz,ph+1):
                ebpts[:,i,:] = np.zeros((Pwg.shape[0],Pwg.shape[2]))
                mpi = min(p,i)
                for j in range(max(0,i-tp),mpi+1):
                    ebpts[:,i,:] += bezalfs[i][j]*bpts[:,j,:]

            # print('Bezier + t control points of current segment')
            # print(ebpts)
            #End of degree elevating bezier
            #Must remove knot u = U[a] oldr times
            if oldr > 1:
                first = kind - 2
                last = kind
                den = ub - ua
                bet = (ub - Uh[kind-1])/den
                #Knot removal loop
                for tr in range(1,oldr):
                    i = first
                    j = last
                    kj = j - kind + 1
                    #Loop and compute the new
                    #Control points for one removal step
                    while (j - i) > tr:
                        if i < cind:
                            alf = (ub - Uh[i])/(ua - Uh[i])
                            Qwg[:,i,:] = alf*Qwg[:,i,:] + (1.0 - alf)*Qwg[:,i-1,:]

                        if j >= lbz:
                            if (j - tr) <= (kind - ph + oldr):
                                gam = (ub - Uh[j - tr])/den
                                ebpts[:,kj,:] = gam*ebpts[:,kj,:] + (1.0 - gam)*ebpts[:,kj+1,:]
                            else:
                                ebpts[:,kj,:] = bet*ebpts[:,kj,:] + (1.0 - bet)*ebpts[:,kj+1,:]

                        i += 1
                        j -= 1
                        kj -= 1

                    first -= 1
                    last += 1

            #End of removing knot u = U[a]
            #Load the knot ua
            if a != p:
                for i in range(0,ph-oldr):
                    Uh[kind] = ua
                    kind += 1

            #Load control points into Qw
            for j in range(lbz,rbz+1):
                Qwg[:,cind,:] = ebpts[:,j,:]
                cind += 1

            if b < mN:
                #Set up for the next loop through
                for j in range(0,r):
                    bpts[:,j,:] = Nextbpts[:,j,:]

                for j in range(r,p+1):
                    bpts[:,j,:] = Pwg[:,b-p+j,:]

                a = b
                b += 1
                ua = ub
            else:
                #End knot
                for i in range(0,ph+1):
                    Uh[kind+i] = ub

        #End of while loop (b < m)
        nh = mh - ph - 1

    if dir == "VDIR":
        m = Pwg.shape[2] - 1 #m is the number of columns of Pwg minus 1
        mM = m + q + 1 #mM is the m value from the V knot vector

        qh = q + tq
        qh2 = qh//2

        s0 = len(np.unique(V))-2

        bezalfs = np.zeros((q+tq+1,q+1))
        bpts = np.zeros((Pwg.shape[0],Pwg.shape[1],q+1))
        ebpts = np.zeros((Pwg.shape[0],Pwg.shape[1],q+tq+1))
        Nextbpts = np.zeros((Pwg.shape[0],Pwg.shape[1],q-1))
        alfs = np.zeros(q-1)
        Vh = np.zeros(len(V) + len(np.unique(V))*tq)
        Qwg = np.zeros((Pwg.shape[0],Pwg.shape[1],len(Vh) - (q+tq) - 1))

        Uh = U
        ph = p

        #Compute bezier degree elevation coefficients
        bezalfs[0][0] = 1.0
        bezalfs[qh][q] = 1.0

        for i in range(1,qh2+1):
            inv = 1.0/binomialCoefficient(qh,i)
            mqi = min(q,i)
            for j in range(max(0,i-tq),mqi+1):
                bezalfs[i][j] = inv*binomialCoefficient(q,j)*binomialCoefficient(tq,i-j)

        for i in range(qh2+1,qh):
            mqi = min(q,i)
            for j in range(max(0,i-tq),mqi+1):
                bezalfs[i][j] = bezalfs[qh-i][q-j]

        mh = qh
        kind = qh + 1
        r = -1
        a = q
        b = q + 1
        cind = 1
        va = Vh[0]
        # Qw[0,:] = Pw[0,:]
        Qwg[:,:,0] = Pwg[:,:,0]

        for i in range(0,qh+1):
            Vh[i] = va

        #Initialize first bezier segment
        for i in range(0,q+1):
            # bpts[i,:] = Pw[i,:]
            bpts[:,:,i] = Pwg[:,:,i]

        #Big loop through knot vector
        while b < mM:
            i = b
            while b < mM and abs(V[b+1] - V[b]) < 1e-5:
                b += 1

            mul = b - i + 1
            mh += mul + tq
            vb = V[b]
            oldr = r
            r = q - mul
            #Insert knot u(b) r times
            if oldr > 0:
                lbz = (oldr + 2)//2
            else:
                lbz = 1

            if r > 0:
                rbz = qh - (r + 1)//2
            else:
                rbz = qh

            #Insert knot to get bezier segment
            if r > 0:
                numer = vb - va
                for k in range(q,mul,-1):
                    alfs[k-mul-1] = numer/(V[a+k] - va)

                for j in range(1,r+1):
                    save = r - j
                    s = mul + j
                    for k in range(q,s-1,-1):
                        bpts[:,:,k] = alfs[k-s]*bpts[:,:,k] + (1.0 - alfs[k-s])*bpts[:,:,k-1]

                    Nextbpts[:,:,save] = bpts[:,:,q]

            #End of insert knot
            #Degree elevate bezier
            #Only points lbz,...,ph are used below
            # ebpts = np.zeros((2,p+t+1))
            for i in range(lbz,qh+1):
                ebpts[:,:,i] = np.zeros((Pwg.shape[0],Pwg.shape[1]))
                mqi = min(q,i)
                for j in range(max(0,i-tq),mqi+1):
                    ebpts[:,:,i] += bezalfs[i][j]*bpts[:,:,j]

            # print('Bezier + t control points of current segment')
            # print(ebpts)
            #End of degree elevating bezier
            #Must remove knot u = U[a] oldr times
            if oldr > 1:
                first = kind - 2
                last = kind
                den = vb - va
                bet = (vb - Vh[kind-1])/den
                #Knot removal loop
                for tr in range(1,oldr):
                    i = first
                    j = last
                    kj = j - kind + 1
                    #Loop and compute the new
                    #Control points for one removal step
                    while (j - i) > tr:
                        if i < cind:
                            alf = (vb - Vh[i])/(va - Vh[i])
                            Qwg[:,:,i] = alf*Qwg[:,:,i] + (1.0 - alf)*Qwg[:,:,i-1]

                        if j >= lbz:
                            if (j - tr) <= (kind - qh + oldr):
                                gam = (vb - Vh[j - tr])/den
                                ebpts[:,:,kj] = gam*ebpts[:,:,kj] + (1.0 - gam)*ebpts[:,:,kj+1]
                            else:
                                ebpts[:,:,kj] = bet*ebpts[:,:,kj] + (1.0 - bet)*ebpts[:,:,kj+1]

                        i += 1
                        j -= 1
                        kj -= 1

                    first -= 1
                    last += 1

            #End of removing knot v = V[a]
            #Load the knot va
            if a != q:
                for i in range(0,qh-oldr):
                    Vh[kind] = va
                    kind += 1

            #Load control points into Qw
            for j in range(lbz,rbz+1):
                Qwg[:,:,cind] = ebpts[:,:,j]
                cind += 1

            if b < mM:
                #Set up for the next loop through
                for j in range(0,r):
                    bpts[:,:,j] = Nextbpts[:,:,j]

                for j in range(r,q+1):
                    bpts[:,:,j] = Pwg[:,:,b-q+j]

                a = b
                b += 1
                va = vb
            else:
                #End knot
                for i in range(0,qh+1):
                    Vh[kind+i] = vb

        #End of while loop (b < m)
        nh = mh - qh - 1

    return Uh,Vh,ph,qh,Qwg

################ GENERAL REFINEMENT ####################

# Knot refinement
def hRefinement(dir,X,U,V,p,q,Pwg):
    Uh,Vh,Qwhg = knotRefinement(dir,X,U,V,p,q,Pwg)
    return Uh,Vh,p,q,Qwhg

# Degree Elevation
def pRefinement(dir,U,V,p,q,Pwg):
    tp = 1
    tq = 1
    Up,Vp,pp,qp,Qwpg = surfaceDegreeElevation(dir,U,V,p,q,Pwg,tp,tq)

    return Up,Vp,pp,qp,Qwpg

# K refinement
def kRefinement(dir,X,U,V,p,q,Pwg):
    Up,Vp,pp,qp,Pwpg = pRefinement(dir,U,V,p,q,Pwg)
    Uk,Vk,pk,qk,Pwkg = hRefinement(dir,X,Up,Vp,pp,qp,Pwpg)

    return Uk,Vk,pk,qk,Pwkg

def surfaceRefinement(surface,rfnnum,rfntype,rfndir,paramlist=None):
    Uin,Vin,pin,qin,Pin,win = surface.retrieveSurfaceInformation()

    for irfn in range(rfnnum):
        if paramlist is not None:
            # [paramstart,paramend,numknots]
            parlist = paramlist[numindex]
            X = np.linspace(parlist[0],parlist[1],parlist[2])

        Pw = rbs.weightedControlPoints(Pin,win)
        Pwgrid = rbs.listToGridControlPoints(Pw,Uin,Vin,pin,qin)

        if rfndir == "U":
            dir = "UDIR"
            if paramlist == None:
                Ured = Uin[pin:-pin]
                X = 0.5*(Ured[0:-1] + Ured[1:])
        elif rfndir == "V":
            dir = "VDIR"
            if paramlist == None:
                Vred = Vin[qin:-qin]
                X = 0.5*(Vred[0:-1] + Vred[1:])
        else:
            print("Wrong parameter direction")

        if rfntype == 'h':
            Uout,Vout,pout,qout,Qwout = hRefinement(dir,X,Uin,Vin,pin,qin,Pwgrid)
        elif rfntype == 'p':
            Uout,Vout,pout,qout,Qwout = pRefinement(dir,Uin,Vin,pin,qin,Pwgrid)
        elif rfntype == 'k':
            Uout,Vout,pout,qout,Qwout = kRefinement(dir,X,Uin,Vin,pin,qin,Pwgrid)
        else:
            print("Invalid option")

        Qw = rbs.gridToListControlPoints(Qwout)
        Pout,wout = rbs.geometricControlPoints(Qw)

        Uin = Uout
        Vin = Vout
        pin = pout
        qin = qout
        Pin = Pout
        win = wout

    print("{0} {1} refinements in {2} direction".format(rfnnum,rfntype,rfndir))

    surface.updateSurfaceInformation(Uout,Vout,pout,qout,Pout,wout)

def localPatchRefinement(multisurface,patchlist,numreflist,reflist,dirlist):
    multiU,multiV,multip,multiq,multiP,multiw,globalPatchIndices = multisurface.retrieveSurfaceInformation()

    for idpatch in range(len(patchlist)):
        print('Patch #',idpatch)

        idPatchToRefine = patchlist[idpatch]
        # Select info of the desired patch to refine
        Uinit = multiU[idPatchToRefine]
        Vinit = multiV[idPatchToRefine]
        pinit = multip[idPatchToRefine]
        qinit = multiq[idPatchToRefine]
        Pinit = multiP[idPatchToRefine]
        winit = multiw[idPatchToRefine]
        globalPatchIndicesinit = globalPatchIndices[idPatchToRefine]

        numrefs = numreflist[idPatchToRefine]
        reftypes = reflist[idPatchToRefine]
        refdir = dirlist[idPatchToRefine]

        surface_i = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)

        for dir in refdir:
            for ref in reftypes:
                surfaceRefinement(surface_i,numrefs,ref,dir)

        Uinp,Vinp,pinp,qinp,Pinp,winp = surface_i.retrieveSurfaceInformation()

        # Insert new information about refined patch
        multiU[idPatchToRefine] = Uinp
        multiV[idPatchToRefine] = Vinp
        multip[idPatchToRefine] = pinp
        multiq[idPatchToRefine] = qinp
        multiP[idPatchToRefine] = Pinp
        multiw[idPatchToRefine] = winp
    # End idpatch loop

    multisurface.updateMultiPatchInformation(multiU,multiV,multip,multiq,multiP,multiw)
# End function
