# Python libraries
import numpy as np

# Local project
from src.spline_functions.basisFunctions import basis_function
from src.spline_functions.basisFunctions import der_basis_function

def binomial(a: int, b: int) -> float:
    bc = 1.0
    for j in range(1, b+1):
        bc *= ((a+1-j)//j)
    return bc

#####################################################
################# CONTROL POINTS ####################
#####################################################

class NURBSObject:
    def __init__(self) -> None:
        pass

    #Convert from real control points to homogeneous ones
    def weightedControlPoints(self, P: np.ndarray, w: np.ndarray) -> np.ndarray:
        Pw = np.hstack((P,np.ones((P.shape[0],1))))
        Pw *= w
        return Pw

    #Convert from homogeneous to real projection
    def geometricControlPoints(self, Pw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        w = Pw[:,-1]
        w = np.reshape(w,(len(w),1))
        P = Pw[:,:-1]/w
        return P, w

    #Convert list of control points to a spatial grid
    def listToGridControlPoints(self, Pl: np.ndarray, U: np.ndarray,
                                V: np.ndarray, p: int, q: int) -> np.ndarray:
        #Number of control points in the U direction
        NU = len(U) - p - 1

        #Number of control points in the V direction
        NV = len(V) - q - 1

        # print("NU:",NU)
        # print("NV:",NV)

        Pg = np.zeros((Pl.shape[1],NU,NV))

        for i in range(0,Pl.shape[1]):
            Pg[i,:,:] = np.reshape(Pl[:,i],(NU,NV),order='F')

        return Pg

    def gridToListControlPoints(self, Pg: np.ndarray) -> np.ndarray:
        Pl = np.zeros((Pg.shape[1]*Pg.shape[2],Pg.shape[0]))

        for i in range(0,Pg.shape[0]):
            Pl[:,i] = np.reshape(Pg[i,:,:],(Pg.shape[1]*Pg.shape[2]),order='F')

        return Pl

    #########################################################################
    ###################### RATIONAL BASIS FUNCTIONS #########################
    #########################################################################

    def nonZeroIndicesElement(self, uspan: int, vspan: int,
                              p: int, q: int, nu: int) -> list[int]:
        idr = []

        for ll in range(0,q+1):
            for k in range(0,p+1):
                i = (vspan-q+ll)*(nu+1) + (uspan-p+k)
                idr.append(i)

        return idr

    def nonZeroIndiceBoundary(self, apt: np.ndarray, bpt: np.ndarray,
                              uspan: int, vspan: int, p: int, q: int,
                              nu: int) -> list[int]:
        idr = []
        cpt = bpt - apt

        if abs(cpt[0]) < 1e-5:
            for ll in range(0,q+1):
                if abs(apt[0]) < 1e-5:
                    i = (vspan - q + ll)*(nu + 1) + (uspan - p)
                elif abs(apt[0] - 1.0) < 1e-5:
                    i = (vspan - q + ll)*(nu + 1) + (uspan)
                idr.append(i)
        elif abs(cpt[1]) < 1e-5:
            for k in range(0,p+1):
                if abs(apt[1]) < 1e-5:
                    i = (vspan - q)*(nu + 1) + (uspan - p + k)
                elif abs(apt[1] - 1.0) < 1e-5:
                    i = (vspan)*(nu + 1) + (uspan - p + k)
                idr.append(i)
        else:
            print("Wrong dimensions")

        return idr

    def rationalCurveDerivative(self, aders: np.ndarray, wders: np.ndarray,
                                d: int) -> np.ndarray:
        Ck = np.zeros((d+1,2))

        for k in range(0,d+1):
            v = aders[k,:]
            for i in range(1,k+1):
                v -= binomial(k,i)*wders[i]*Ck[k-i,:]

            Ck[k,:] = v/wders[0]

        return Ck

    def univariateRationalDerivative(self, aders: np.ndarray,
                                     wders: np.ndarray, d: int) -> np.ndarray:
        drat = np.zeros((d+1,aders.shape[1]))

        for k in range(0,d+1):
            v = aders[k,:]
            for i in range(1,k+1):
                v -= binomial(k,i)*wders[i]*drat[k-i,:]

            drat[k,:] = v/wders[0]

        return drat

    def bivariateRationalFunction(self, p: int, q: int, uspan: int, vspan: int,
                                  u: float, v: float, U: np.ndarray,
                                  V: np.ndarray, Pw: np.ndarray) -> np.ndarray:
        Nu = basis_function(uspan,u,p,U)
        Nv = basis_function(vspan,v,q,V)

        R = np.zeros((1,(p+1)*(q+1)))
        i = 0

        for ll in range(0,q+1):
            for k in range(0,p+1):
                R[0][i] = Nu[k]*Nv[ll]*Pw[-1,uspan-p+k,vspan-q+ll]
                i += 1

        R /= np.sum(R)

        return R

    def bivariateRationalGradient(self, mu: int, mv: int, p: int,
                                  q: int, uspan: int, vspan: int,
                                  u: float, v: float, U: np.ndarray,
                                  V: np.ndarray, Pw: np.ndarray) -> np.ndarray:
        # Derivative order
        d = 1

        # The first row for dNu and dNv has the shape functions
        # Nu and Nv respectively

        dNu = der_basis_function(uspan,u,mu,p,U,d)
        dNv = der_basis_function(vspan,v,mv,q,V,d)

        # The first row for dR has the shape function
        # The second row for dR has the derivative w.r.t u
        # The third row for dR has the derivative w.r.t v
        dR = np.zeros((3,(p+1)*(q+1)))

        # The first row for dA and dW has the derivative w.r.t u
        # The second row for dA and dW has the derivative w.r.t v
        dA = np.zeros((2,(p+1)*(q+1)))
        # dW = np.zeros((2,1))
        i = 0

        for ll in range(0,q+1):
            for k in range(0,p+1):
                iu = uspan - p + k
                jv = vspan - q + ll
                dR[0][i] = dNu[0][k]*dNv[0][ll]*Pw[-1,iu,jv]
                dA[0][i] = dNu[d][k]*dNv[0][ll]*Pw[-1,iu,jv]
                dA[1][i] = dNu[0][k]*dNv[d][ll]*Pw[-1,iu,jv]
                i += 1

        W = np.sum(dR[0,:])
        dW = np.sum(dA,axis = 1)
        biN = dR[0,:]/W

        dR[0,:] = biN
        dR[1,:] = (dA[0,:] - dW[0]*biN)/W
        dR[2,:] = (dA[1,:] - dW[1]*biN)/W

        return dR
