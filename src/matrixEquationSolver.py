# Python libraries
import numpy as np
import numpy.linalg

################ MATRIX EQUATION SOLUTION ####################

def dirichletBCEnforcement(phenomenon,K,F,dirichletconds):
    enforcednodes = []
    restricteddofs = []
    values = []

    for drchcond in dirichletconds:
        if phenomenon == "Elasticity":
            if len(drchcond) == 5:
                # Beware the first index
                ipatch = drchcond[0]
                inode = drchcond[0]
                value = drchcond[1][0]
                localdofs = drchcond[1][1]
            else:
                inode = drchcond[0]
                value = drchcond[1][0]
                localdofs = drchcond[1][1]
            # End if

            if len(localdofs) == 2:
                # On clamped condition, the node and the value
                # are replicated as many spatial dimensions are
                enforcednodes.append(2*inode)
                enforcednodes.append(2*inode)
                values.append(value)
                values.append(value)
            elif len(localdofs) == 1:
                enforcednodes.append(2*inode)
                values.append(value)
            else:
                print("Wrong restriction")
            # End if
            restricteddofs += localdofs
        elif phenomenon == "Heat":
            if len(drchcond) == 4:
                ipatch = drchcond[0]
                inode = drchcond[1]
                localdofs = drchcond[1]
                value = drchcond[2]
            else:
                inode = drchcond[0]
                localdofs = drchcond[0]
                value = drchcond[1]
            # End if

            enforcednodes.append(inode)
            values.append(value)
            restricteddofs.append(localdofs)
        else:
            print("Wrong physics")
            enforcednodes.append(2*inode)
            values.append(value)
        # End if
    #End dirichlet conditions loop

    if phenomenon == "Elasticity":
        # Calculating the global dofs to be removed
        enforcednodes = np.array(enforcednodes)
        restricteddofs = np.array(restricteddofs)
        restricteddofs += enforcednodes
    elif phenomenon == "Heat":
        restricteddofs = np.array(restricteddofs)
    else:
        print("Check physics selection")
        restricteddofs = np.array(restricteddofs)
    # End if

    print("First reduction")
    Fred = np.delete(F,restricteddofs,0)
    Kred = np.delete(K,restricteddofs,0)

    print("Modification of Freduced")
    for i in range(len(restricteddofs)):
        Kcol = Kred[:,restricteddofs[i]]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        # print(values[i])
        Fred -= Kcol*values[i]

    print("Second reduction")
    Kred = np.delete(Kred,restricteddofs,1)

    totaldofs = np.arange(F.shape[0])

    return Kred,Fred,totaldofs,restricteddofs,values

def solveMatrixEquations(phenomenon,Kred,Fred,totaldofs,remdofs,values):
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
    # for idof in range(len(remdofs)):
        # dtotal[remdofs[idof]] = values[idof]

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
