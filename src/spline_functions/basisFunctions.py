#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:11:32 2020
Modified on Sun Aug 02 15:29:17 2020

@author: hernando
@collaborator: paternina
"""
#"NURBS"
import numpy as np
import sys

#tolerance for knot vector
tol = 1e-5

def generateUniformKnotVector(N,p):
    if N < p:
        print ("Fatal error")
        return np.zeros(p+1)
    else:
        n = N - 1
        M = N + p + 1
        m = M - 1
        uknot = np.zeros(M)

        for j in range(p+1,m-p):
            uknot[j] = (j-p)/(n)

        for j in range(n+1,M):
            uknot[j] = 1.0

        return uknot

def knotGeneratorChord(n,p,tv):
    if n<p:
        print ("Fatal error")
        sys.exit("No such vector curve exists, please try n greater than p")
        return np.zeros(p+1)
    else:
        ustart = np.zeros(p+1)
        uend = np.ones(p+1)
        umid = np.zeros(n-p)
        for j in range(1,n-p+1):
            ui = 0
            for i in range(j,j+p):
                ui += tv[i]
            ui *= (1.0/p)
            umid[j-1] = ui
        return np.concatenate([ustart,umid,uend])

def findKnotInterval(n: int, p: int, u: float, U: np.ndarray) -> int:
    """
    U: knot vector
    u: parameter between 0 and 1
    """
    # Special case
    if abs(u - U[n+1]) < 1e-5:
        return n

    # Starting binary search
    low = p
    high = n + 1
    mid = (low + high)//2

    while u < U[mid] or (u > U[mid+1] or abs(u - U[mid+1]) < 1e-5):
        if u < U[mid]:
            high = mid
        else:
            low = mid

        mid = (low + high)//2

    return mid

####################################################
#################BASIS FUNCTIONS####################
####################################################

def basisFunction(i: int, u: float, m: int, p: int, U: np.ndarray) -> np.ndarray:
    """
    Computing the basis functions Nip (modified version)
    
    ----------
    Input: i,u,m,p,U
    Output: nbasis

    Variables
    ----------

    U: knot vector
    u: parameter between 0 and 1
    """
    N = np.zeros(p+1)
    left = np.zeros(p+1)
    right = np.zeros(p+1)

    N[0] = 1.0
    for j in range(1,p+1):
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved = 0.0
        for r in range(0,j):
            temp = N[r]/(right[r+1] + left[j-r])
            N[r] = saved + right[r+1]*temp
            saved = left[j-r]*temp

        N[j] = saved

    return N

def derBasisFunction(i,u,m,p,U,nd):
    """
    Computing the derivatives of basis functions Nip (modified version)
    Input: i,u,m,p,U
    Output: ders
    """

    ndu = np.zeros((p+1,p+1))
    a = np.zeros((2,p+1))
    left = np.zeros(p+1)
    right = np.zeros(p+1)
    ders = np.zeros((nd+1,p+1))

    ndu[0][0] = 1.0
    
    for j in range(1,p+1):
        left[j] = u - U[i+1-j]
        right[j] = U[i+j] - u
        saved = 0.0
        for r in range(0,j):
            # Lower triangle
            ndu[j][r] = right[r+1] + left[j-r]
            temp = ndu[r][j-1]/ndu[j][r]
            # Upper triangle
            ndu[r][j] = saved + right[r+1]*temp
            saved = left[j-r]*temp

        ndu[j][j] = saved

    # Load the basis functions
    for j in range(0,p+1):
        ders[0][j] = ndu[j][p]

    # This section computes the derivatives
    # according to Eq (2.9) from The NURBS Book
    for r in range(0,p+1):
        # Loop over function index
        s1 = 0
        s2 = 1
        # Alternate rows in array a
        a[0][0] = 1.0
        # Look to compute kth derivative
        for k in range(1,nd+1):
            d = 0.0
            rk = int(r - k)
            pk = int(p - k)

            if r >= k:
                a[s2][0] = a[s1][0]/ndu[pk+1][rk]
                d = a[s2][0]*ndu[rk][pk]

            if rk >= -1:
                j1 = 1
            else:
                j1 = -rk

            if (r-1) <= pk:
                j2 = k - 1
            else:
                j2 = p - r

            for j in range(j1,j2+1):
                a[s2][j] = (a[s1][j] - a[s1][j-1])/ndu[pk+1][rk+j]
                d += a[s2][j]*ndu[rk+j][pk]

            if r <= pk:
                a[s2][k] = -a[s1][k-1]/ndu[pk+1][r]
                d += a[s2][k]*ndu[r][pk]

            ders[k][r] = d
            # Switch rows
            j = s1
            s1 = s2
            s2 = j
    # Multiply through by the correct factors [Eq (2.9)]
    r = p
    for k in range(1,nd+1):
        for j in range(0,p+1):
            ders[k][j] *= r

        r *= (p-k)

    # derk = ders[k,:]
    return ders


def oneBasisFunction(p: int, U: np.ndarray, i: int, u: float) -> float:
    """
    Compute the basis function Nip
    Input: p, U, i, u
    Output: Nip
    """
    Nip = 0.0
    m = len(U) - 1
    N = np.zeros(m-p-1+1)

    if ((i == 0 and u == U[0]) or (i == m - p - 1 and u == U[m])):
        Nip = 1.0
        return Nip
    
    if (u < U[i] or u >= U[i+p+1]):
        Nip = 0.0
        return Nip
    
    for j in range(0,p+1):
        if (u >= U[i+j] and u < U[i+j+1]):
            N[j] = 1.0
        else:
            N[j] = 0.0

    print(N)

    for k in range(1,p+1):
        if N[0] == 0.0:
            saved = 0.0
        else:
            saved = ((u - U[i])*N[0]/(U[i+k] - U[i]))
        
        for j in range(0, p-k+1):
            Uleft = U[i+j+1]
            Uright = U[i+j+k+1]
            if N[j+1] == 0.0:
                N[j] = saved
                saved = 0.0
            else:
                temp = N[j+1]/(Uright - Uleft)
                N[j] = saved + (Uright - u)*temp
                saved = (u - Uleft)*temp
    
    Nip = N[0]

    return Nip

def test_basis_functions():
    U = np.array([0,0,0,1,2,3,4,4,5,5,5])
    p = 2
    u = 5./2
    i = 4
    Nip = oneBasisFunction(p,U,i,u)
    print(Nip)

def test_basis_functions_2():
    U = np.array([0,0,0,1,2,3,4,4,5,5,5])
    p = 2
    u = 5./2
    m = len(U)-1
    n = m - p - 1
    mid = findKnotInterval(n,p,u,U)
    print(mid)

def test_basis_functions_3():
    U = np.array([0,0,0,1,2,3,4,4,5,5,5])
    p = 2
    u = 5./2
    m = len(U)-1
    n = m - p - 1
    i = 4
    nd = 1
    ders = derBasisFunction(i,u,m,p,U,nd)
    print(ders)

if __name__ == '__main__':
    test_basis_functions_3()
