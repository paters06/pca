# Python libraries
import numpy as np

################ MATRIX EQUATION SOLUTION ####################

def dirichletBCEnforcement(M,K,F,enforceddof,enforcedvalues):
    print("First reduction")

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

def solveMatrixEquations(Kred,Fred,totaldofs,remdofs,values):
    # Checking full rank in matrix
    mRank = np.linalg.matrix_rank(Kred)
    mRows = Kred.shape[0]
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
    from scipy.sparse.linalg import eigs,eigsh

    eigenvalues,eigenvectors = eigsh(Kred,M=Mred,k=10,which='SM')
    print(eigenvalues)
    # print(np.real(eigenvalues))

    redeigensol = eigenvectors[:,numeig].reshape(-1,1)
    # print(redeigensol.T)

    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    eigensol = np.zeros((totaldofs.shape[0],1))
    eigensol[reduceddofs,:] = redeigensol

    values = np.reshape(values,(len(values),1))
    eigensol[remdofs,:] = values

    return eigensol
# End Function
