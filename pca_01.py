#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:11:32 2020
Modified on Sun Aug 02 15:29:17 2020

@author: hernando
@collaborator: paternina
"""
"NURBS"

import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import sys

def knotVector(n,p):

    ustart = np.zeros(p+1)
    
    if n<p:
        print ("Fatal error")
        sys.exit("No such vector curve exists, please try n greater than p")
        return ustart
    else:
        umid = np.arange(1,n-p)/(n-p+1)
        uend = np.ones(p+1)
        return np.concatenate([ustart,umid,uend])

#U: knot vector
#u: parameter between 0 and 1
def findKnotInterval(U,u):
    for i in range(len(U)-1):
        if U[i]<=u and U[i+1]>u:
            index = i

    if U[-1] == U.max():
        index = -1

    return index

def uGenerator(a,b,step):
    if 0<=a and b<=1:
        u = np.arange(a,b,step,dtype=np.float64)
    else:
        print('Ranks do not belong to the domain of the function')
        sys.exit("Error message")
    return u  

def nFunction(U,p,u):
    m = len(U) - 1
    for pi in range(p+1):
        # print("p-Order: ",pi)
        if pi != 0:
            nbas = np.zeros(len(nbasis)-1)
            for j in range(m-pi):
                
                num1 = u - U[j]
                den1 = U[j+pi] - U[j]

                num2 = U[j+pi+1] - u
                den2 = U[j+pi+1] - U[j+1]

                if den1==0:
                    A = 0
                else:
                    A = num1/den1
                if den2==0:
                    B = 0
                else:
                    B = num2/den2

                nbas[j] = A*nbasis[j] + B*nbasis[j+1]

            nbasis = nbas
        else:
            index = findKnotInterval(U,u)
            nbasis = np.zeros(len(U)-1)
            for i in range(len(nbasis)):
                if U[i]<=u and U[i+1]>u:
                    nbasis[i] = 1.0

        if u == U.max():
            nbasis[-1:]=1

    return nbasis
        
def bSplineCurve(U,p,px,py):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    for i in range(len(urank)):
        nVec = nFunction(U,p,urank[i])
        cx[i] = px@nVec
        cy[i] = py@nVec
    return cx,cy

def nurbsCurve(U,p,px,py,w):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))
    
    for i in range(len(urank)):
        nVec = nFunction(U,p,urank[i])
        cx[i] = (px@(nVec*w))/(nVec@w)
        cy[i] = (py@(nVec*w))/(nVec@w)
    return cx,cy

def bSplineSurface(U,V,p,q,px,py,pz):
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

            cx[i,j] = nVecj.transpose()@px@nVeci
            cy[i,j] = nVecj.transpose()@py@nVeci
            cz[i,j] = nVecj.transpose()@pz@nVeci
    return cx,cy,cz

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

            cx[i,j] = (nVecj.transpose()@(px*w)@nVeci)/(nVecj.transpose()@(w)@nVeci)
            cy[i,j] = (nVecj.transpose()@(py*w)@nVeci)/(nVecj.transpose()@(w)@nVeci)
            cz[i,j] = (nVecj.transpose()@(pz*w)@nVeci)/(nVecj.transpose()@(w)@nVeci)
    return cx,cy,cz

def plotCurve2d(cx,cy,px,py):
    fig = plt.figure()
    plt.plot(cx,cy)
    plt.plot(px,py,'ro')
    plt.plot(px,py)
    plt.show()

def plotting3d(cx,cy,cz,*argv):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(cx, cy, cz, 'gray')
    if argv != ():
        if argv[0]=='on':
            if len(argv)==4:
                ax.plot3D(argv[1],argv[2],argv[3], 'red')
            elif len(argv)> 4:
                sys.exit("Too much arguments, please delete one or more")
            else:
                sys.exit("Missing arguments to plot control points")
                
def plottingSurface(cx,cy,cz,*argv):
    fig = plt.figure()
    ax = plt.axes(projection = '3d')
    ax.contour3D(cx, cy, cz, 50, cmap = 'binary')
    if len(argv)==3:
        ax.plot_wireframe(argv[0], argv[1], argv[2], color = 'red')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z');
    plt.show()
