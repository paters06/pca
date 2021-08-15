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

#U: knot vector
#u: parameter between 0 and 1
def findKnotInterval(n,p,u,U):
    # Special case
    if abs(u - U[n+1]) < 1e-5:
        return n
    # End if

    # Starting binary search
    low = p
    high = n + 1
    mid = (low + high)//2

    while u < U[mid] or (u > U[mid+1] or abs(u - U[mid+1]) < 1e-5):
        if u < U[mid]:
            high = mid
        else:
            low = mid
        # End if

        mid = (low + high)//2
    # End while loop

    return mid
# End function

####################################################
#################BASIS FUNCTIONS####################
####################################################

def basisFunction(i,u,m,p,U):
    """
    Computing the basis functions Nip (modified version)
    Input: i,u,m,p,U
    Output: nbasis
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

#    derk = ders[k,:]
    return ders
