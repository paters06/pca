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
                ipatch = drchcond[0]
                inode = drchcond[1]
                restriction = drchcond[2]
                localdofs = drchcond[3]
                value = drchcond[4]
            else:
                inode = drchcond[0]
                restriction = drchcond[1]
                localdofs = drchcond[2]
                value = drchcond[3]
            # End if

            if restriction == "C":
                # On clamped condition, the node and the value
                # are replicated as many spatial dimensions are
                enforcednodes.append(2*inode)
                enforcednodes.append(2*inode)
                values.append(value)
                values.append(value)
            elif restriction == "S":
                enforcednodes.append(2*inode)
                values.append(value)
            else:
                print("Wrong restriction")
            # End if

            # print(localdofs)
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

        # print(localdofs)
        # restricteddofs.append(localdofs)
    #End dirichlet conditions loop
    # print(restricteddofs)

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

    return Kred,Fred,restricteddofs,totaldofs

def solveMatrixEquations(phenomenon,Kred,Fred,totaldofs,remdofs,dirichletconds):
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

    if phenomenon == "Elasticity":
        for drchcond in dirichletconds:
            dtotal[drchcond[0]] = drchcond[3]

        dx = dtotal[0::2]
        dy = dtotal[1::2]
        D = np.hstack((dx,dy))

        return dtotal,D
    elif phenomenon == "Heat":
        for drchcond in dirichletconds:
            dtotal[drchcond[0]] = drchcond[1]

        return dtotal,dtotal
    else:
        print("Check physics selection")
        return dtotal,dtotal
    # End if
