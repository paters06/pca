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
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import sys

#tolerance for knot vector
tol = 1e-5

def knotGeneratorUniform(n,p):
    if n<p:
        print ("Fatal error")
        sys.exit("No such vector curve exists, please try n greater than p")
        return np.zeros(p+1)
    else:
        ustart = np.zeros(p+1)
        umid = np.zeros(n-p)
        for j in range(1,n-p+1):
            umid[j-1] = j/(n-p+1)
        uend = np.ones(p+1)
        return np.concatenate([ustart,umid,uend])

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
                # print(tv[i])
                ui += tv[i]
            ui *= (1.0/p)
            umid[j-1] = ui
        return np.concatenate([ustart,umid,uend])

#U: knot vector
#u: parameter between 0 and 1
def findKnotInterval(U,u):
    for i in range(len(U)-1):
        if U[i] < (u + tol) and (u + tol) < (U[i+1]):
            nbasis[i] = 1.0

        if abs(u - U.max())<tol:
            if U[i] < (u - tol) and (u - tol) < U[i+1]:
                nbasis[i] = 1.0
    return index

####################################################
#################BASIS FUNCTIONS####################
####################################################

def nFunction(U,p,u):
    m = len(U) - 1
    for pi in range(p+1):
        # print("p-Order: ",pi)
        if pi != 0:
            nbas = np.zeros((1,len(nbasis[0])-1))
            for j in range(m-pi):

                num1 = u - U[j]
                den1 = U[j+pi] - U[j]

                num2 = U[j+pi+1] - u
                den2 = U[j+pi+1] - U[j+1]

                if abs(den1)<tol:
                    A = 0.0
                else:
                    A = num1/den1

                if abs(den2)<tol:
                    B = 0.0
                else:
                    B = num2/den2

                nbas[0][j] = A*nbasis[0][j] + B*nbasis[0][j+1]

            nbasis = nbas
        else:
            nbasis = np.zeros((1,len(U)-1))
            for i in range(len(nbasis[0])):
                if U[i] < (u + tol) and (u + tol) < (U[i+1]):
                    nbasis[0][i] = 1.0

                if abs(u - U.max())<tol:
                    if U[i] < (u - tol) and (u - tol) < U[i+1]:
                        nbasis[0][i] = 1.0
    return nbasis

def derNFunction(U,p,u):
    m = len(U) - 1
    dnbasis = nFunction(U,p-1,u)
    dnbas = np.zeros((1,len(dnbasis[0])-1))
    for j in range(m-p):

        num1 = 1.0
        den1 = U[j+p] - U[j]

        num2 = 1.0
        den2 = U[j+p+1] - U[j+1]

        if abs(den1)<tol:
            A = 0.0
        else:
            A = num1/den1

        if abs(den2)<tol:
            B = 0.0
        else:
            B = num2/den2

        dnbas[0][j] = p*(A*dnbasis[0][j] - B*dnbasis[0][j+1])

    dnbasis = dnbas
    return dnbasis

def secDerNFunction(U,p,u):
    m = len(U) - 1
    d2nbasis = derNFunction(U,p-1,u)
    d2nbas = np.zeros((1,len(d2nbasis[0])-1))
    for j in range(m-p):

        num1 = 1.0
        den1 = U[j+p] - U[j]

        num2 = 1.0
        den2 = U[j+p+1] - U[j+1]

        if abs(den1)<tol:
            A = 0.0
        else:
            A = num1/den1

        if abs(den2)<tol:
            B = 0.0
        else:
            B = num2/den2

        d2nbas[0][j] = p*(A*d2nbasis[0][j] - B*d2nbasis[0][j+1])

    d2nbasis = d2nbas
    return d2nbasis

####################################################
####################B-SPLINES#######################
####################################################

################B-SPLINES CURVES####################

def bSpline1DCurve(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))

    for i in range(len(urank)):
        nVec = nFunction(U,p,urank[i])
        cx[i] = nVec@px.T
    return cx

def bSplineCurve(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        nVec = nFunction(U,p,urank[i])
        cx[i] = nVec@px.T
        cy[i] = nVec@py.T
    return cx,cy

def bSplineCurveDerivative(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cprimex = np.zeros(len(urank))
    cprimey = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        dnVec = derNFunction(U,p,urank[i])
        cprimex[i] = dnVec@px.T
        cprimey[i] = dnVec@py.T
    return cprimex,cprimey

################B-SPLINES SURFACES###################

def bSplineSurface(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = nFunction(U,p,urank[i])
            nVecj = nFunction(V,q,vrank[j])

            cx[i,j] = nVecj.T@px@nVeci
            cy[i,j] = nVecj.T@py@nVeci
            cz[i,j] = nVecj.T@pz@nVeci
    return cx,cy,cz

def bSplineSurfaceDerivativeU(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            dnVeci = derNFunction(U,p,urank[i])
            nVecj = nFunction(V,q,vrank[j])

            cprimex[i,j] = nVecj.T@px@dnVeci
            cprimey[i,j] = nVecj.T@py@dnVeci
            cprimez[i,j] = nVecj.T@pz@dnVeci
    return cprimex,cprimey,cprimez

def bSplineSurfaceDerivativeV(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = nFunction(U,p,urank[i])
            dnVecj = derNFunction(V,q,vrank[j])

            cprimex[i,j] = dnVecj.T@px@nVeci
            cprimey[i,j] = dnVecj.T@py@nVeci
            cprimez[i,j] = dnVecj.T@pz@nVeci
    return cprimex,cprimey,cprimez

def bSplineSurfaceDerivativeUV(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            dnVeci = derNFunction(U,p,urank[j])
            dnVecj = derNFunction(V,q,vrank[j])

            cprimex[i,j] = dnVecj.T@px@dnVeci
            cprimey[i,j] = dnVecj.T@py@dnVeci
            cprimez[i,j] = dnVecj.T@pz@dnVeci
    return cprimex,cprimey,cprimez

####################################################
######################NURBS#########################
####################################################

##################NURBS CURVES######################

def nurbsCurve(U,p,P,w):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        nVec = nFunction(U,p,urank[i])
        cx[i] = (px@(nVec*w).transpose())/(nVec@w.transpose())
        cy[i] = (py@(nVec*w).transpose())/(nVec@w.transpose())
    return cx,cy

#################NURBS SURFACES#####################

def nurbsSurface(U,V,p,q,px,py,pz,w):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = nFunction(U,p,urank[i])
            nVecj = nFunction(V,q,vrank[j])

            cx[i,j] = (nVecj@(px*w)@nVeci.transpose())/(nVecj@(w)@nVeci.transpose())
            cy[i,j] = (nVecj@(py*w)@nVeci.transpose())/(nVecj@(w)@nVeci.transpose())
            cz[i,j] = (nVecj@(pz*w)@nVeci.transpose())/(nVecj@(w)@nVeci.transpose())
    return cx,cy,cz
