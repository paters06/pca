# Python libraries
import numpy as np

# Local project
import src.debugScripts as dbg_scrpt

################ MATRIX EQUATION SOLUTION ####################

def dirichletBCEnforcement_Modified(M,K,F,enforceddof,enforcedvalues):
    mRows = K.shape[0]
    mCols = K.shape[1]

    # enforceddof = np.array(enforceddof)
    # enforceddof,iddof = np.unique(enforceddof,return_index=True)
    # enforcedvalues = np.array(enforcedvalues)[iddof]

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
    # End for loop

    print("Third modification")
    Kmod[:,enforceddof] = np.zeros((mRows,numdof))
    Kmod[enforceddof,enforceddof] = np.ones((numdof))

    Mmod[:,enforceddof] = np.zeros((mRows,numdof))
    Mmod[enforceddof,enforceddof] = np.ones((numdof))

    Fmod[enforceddof] = np.reshape(enforcedvalues,(numdof,1))

    totaldofs = np.arange(Fmod.shape[0])

    return Mmod,Kmod,Fmod,totaldofs
# End function

def dirichletBCEnforcement_Reduced(M,K,F,enforceddof,enforcedvalues):
    print("First reduction")

    # Mred = M.copy()
    # Kred = K.copy()
    # Fred = F.copy()

    Fred = np.delete(F,enforceddof,0)
    Kred = np.delete(K,enforceddof,0)
    Mred = np.delete(M,enforceddof,0)

    print("Modification of Freduced")
    for i in range(len(enforceddof)):
        Kcol = Kred[:,enforceddof[i]]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        Fred -= Kcol*enforcedvalues[i]
    # End for loop

    print("Second reduction")
    Kred = np.delete(Kred,enforceddof,1)
    Mred = np.delete(Mred,enforceddof,1)

    totaldofs = np.arange(F.shape[0])

    return Mred,Kred,Fred,totaldofs
# End function

def solveModifiedMatrixEquations(Kmod,Fmod,totaldofs):
    # Checking full rank in matrix
    mRank = np.linalg.matrix_rank(Kmod)
    mRows = Kmod.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        fullRank = True
    else:
        fullRank = False
    # End if

    # fullRank = True
    if fullRank:
        print("The matrix has full rank. It is invertible")
        dsol = np.linalg.solve(Kmod,Fmod)
    else:
        print("The matrix has not full rank. It is not invertible")
        dsol = np.linalg.lstsq(Kmod,Fmod,rcond=None)[0]
    # End if

    return dsol
# End function

def solveReducedMatrixEquations(Kred,Fred,totaldofs,remdofs,values):
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
    # End if

    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    dtotal = np.zeros((totaldofs.shape[0],1))
    dtotal[reduceddofs,:] = dred

    values = np.reshape(values,(len(values),1))
    dtotal[remdofs,:] = values
    return dtotal
# End Function

def solveReducedEigenvalueProblem(Mred,Kred,numeig,totaldofs,remdofs,values):
    # from scipy.linalg import eig
    # from scipy import linalg
    from scipy.sparse.linalg import eigs,eigsh

    # dbg_scrpt.checkingSymmetricMatrix(Mred)
    # dbg_scrpt.plotSparsity(Kred)
    # for m in Mred:
    #     print(np.sum(m))
    # End for loop

    # eigenvalues,eigenvectors = linalg.eig(Kred,b=Mred)
    # print(eigenvalues)

    # eigenvalues,eigenvectors = eigs(Kred,M=Mred,k=10,which='SM')
    # print(eigenvalues)
    # print(np.real(eigenvalues))

    eigenvalues,eigenvectors = eigsh(Kred,M=Mred,k=10,which='SM')
    # print(eigenvalues)
    # print(np.real(eigenvalues))

    redeigensol = eigenvectors[:,numeig].reshape(-1,1)
    # print(eigensol.T)

    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    eigensol = np.zeros((totaldofs.shape[0],1))
    eigensol[reduceddofs,:] = redeigensol

    values = np.reshape(values,(len(values),1))
    eigensol[remdofs,:] = values

    return eigensol
# End Function
