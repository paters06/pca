import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def checkingSymmetricMatrix(A):
    check = np.allclose(A, A.T, rtol=1e-2, atol=1e-3)
    if check:
        print("The given matrix is symmetric")
    else:
        print("The given matrix is not symmetric")

def checkFullRankMatrix(A):
    mRank = np.linalg.matrix_rank(A)
    mRows = A.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        return True
    else:
        return False

def plotSparsity(A):
    fig = plt.figure()
    plt.spy(A,markersize=5)
    # plt.imshow(A,cmap=cm.viridis)
    # plt.colorbar()
    plt.show()
