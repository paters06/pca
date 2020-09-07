import numpy as np
import pca_01 as pca
import plottingScripts as plts
import matplotlib.pyplot as plt

tol = 1e-5

def knotInsertion(U,p,P,unew):
    Pnew = np.zeros((2,len(P[0]) + 1))
    Unew = np.zeros(len(U) + 1)

    for k in range(len(U)):
        if U[k] < (unew + tol) and (unew + tol) < (U[k+1]):
            kindex = k

        if abs(unew - U.max())<tol:
            if U[k] < (unew - tol) and (unew - tol) < U[k+1]:
                kindex = k

    for i in range(len(Unew)):
        if i <= kindex:
            Unew[i] = U[i]
        elif i == kindex + 1:
            Unew[i] = unew
        elif i > kindex + 1:
            Unew[i] = U[i-1]
        else:
            print('Index error')

    for i in range(len(Pnew[0])):
        if i <= kindex - p:
            Pnew[:,i] = P[:,i]
        elif i >= kindex - p + 1 and i <= kindex:
            alpha_i = (unew - U[i])/(U[i+p] - U[i])
            Pnew[:,i] = alpha_i*P[:,i] + (1.0 - alpha_i)*P[:,i-1]
        elif i >= kindex + 1:
            Pnew[:,i] = P[:,i-1]
        else:
            print('Index error')

    return Unew,Pnew


P = np.array([[0.0,0.5,1.0],[0.0,1.0,0.0]])
p = 2
U = np.array([0,0,0,1,1,1])

Ubar,Pbar = knotInsertion(U,p,P,0.5)

cx,cy = pca.bSplineCurve(U,p,P)
cxnew,cynew = pca.bSplineCurve(Ubar,p,Pbar)

# plts.plotCurve2d(cx,cy,P)
# plts.plotCurve2d(cxnew,cynew,Pbar)
plts.plotNodeInsertion(cx,cy,cxnew,cynew,P,Pbar)
