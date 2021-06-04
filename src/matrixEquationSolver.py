# Python libraries
import numpy as np
import numpy.linalg

import src.debugScripts as dbg_scrpt

################ MATRIX EQUATION SOLUTION ####################

def dirichletBCEnforcement_Modified(M,K,F,enforceddof,enforcedvalues):
    mRows = K.shape[0]
    mCols = K.shape[1]
    numdof = len(enforceddof)

    Mmod = M.copy()
    Kmod = K.copy()
    Fmod = F.copy()

    print("First modification")
    Kmod[enforceddof,:] = np.zeros((numdof,mCols))
    Mmod[enforceddof,:] = np.zeros((numdof,mCols))

    print("Second modification")
    for i in range(numdof):
        # Second modification
        Kcol = Kmod[:,enforceddof[i]]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        Fmod -= Kcol*enforcedvalues[i]

    print("Third modification")
    Kmod[:,enforceddof] = np.zeros((mRows,numdof))
    Kmod[enforceddof,enforceddof] = np.ones((numdof))

    Mmod[:,enforceddof] = np.zeros((mRows,numdof))
    Mmod[enforceddof,enforceddof] = np.ones((numdof))

    Fmod[enforceddof] = np.reshape(enforcedvalues,(numdof,1))

    totaldofs = np.arange(Fmod.shape[0])

    return Mmod,Kmod,Fmod,totaldofs

def dirichletBCEnforcement_Reduced(K,F,enforceddof,enforcedvalues):
    print("First reduction")
    Fred = np.delete(F,enforceddof,0)
    Kred = np.delete(K,enforceddof,0)

    print("Modification of Freduced")
    for i in range(len(enforceddof)):
        Kcol = Kred[:,enforceddof[i]]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        Fred -= Kcol*enforcedvalues[i]

    print("Second reduction")
    Kred = np.delete(Kred,enforceddof,1)

    totaldofs = np.arange(F.shape[0])

    return Kred,Fred,totaldofs

def solveModifiedMatrixEquations(phenomenon,Kmod,Fmod,totaldofs):
    # Checking full rank in matrix
    mRank = np.linalg.matrix_rank(Kmod)
    mRows = Kmod.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        fullRank = True
    else:
        fullRank = False

    # fullRank = True
    if fullRank:
        print("The matrix has full rank. It is invertible")
        dsol = np.linalg.solve(Kmod,Fmod)
    else:
        print("The matrix has not full rank. It is not invertible")
        dsol = np.linalg.lstsq(Kmod,Fmod,rcond=None)[0]

    if phenomenon == "Elasticity":
        dx = dsol[0::2]
        dy = dsol[1::2]
        D = np.hstack((dx,dy))
        return dsol,D
    elif phenomenon == "Heat":
        return dsol,dsol
    else:
        print("Check physics selection")
        return dsol,dsol
    # End if

def solveReducedMatrixEquations(phenomenon,Kred,Fred,totaldofs,remdofs,values):
    # Checking full rank in matrix
    mRank = np.linalg.matrix_rank(Kred)
    mRows = Kred.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        fullRank = True
    else:
        fullRank = False

    # fullRank = True
    if fullRank:
        print("The matrix has full rank. It is invertible")
        dred = np.linalg.solve(Kred,Fred)
    else:
        print("The matrix has not full rank. It is not invertible")
        dred = np.linalg.lstsq(Kred,Fred,rcond=None)[0]

    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    dtotal = np.zeros((totaldofs.shape[0],1))
    dtotal[reduceddofs,:] = dred

    values = np.reshape(values,(len(values),1))
    dtotal[remdofs,:] = values

    if phenomenon == "Elasticity":
        dx = dtotal[0::2]
        dy = dtotal[1::2]
        D = np.hstack((dx,dy))
        return dtotal,D
    elif phenomenon == "Heat":
        return dtotal,dtotal
    else:
        print("Check physics selection")
        return dtotal,dtotal
    # End if
