#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 16:10:13 2020

@author: hernando
@SEGUNDO INTENTO DE HACER B-SPLINES
"""
import numpy as np
import NURBS as bs
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

px=np.array([1,2,3,4,5,6,7])
py=np.array([-1,1,-1,1,-1,1,-1])
pz=np.array([1,2,-1,1,-1,1,1])

p=2

U=bs.knot(len(px)-1,p)

u=bs.ugen(0.0,1.0,0.01)

cx=bs.b_spline(u,px,p)
cy=bs.b_spline(u,py,p)
cz=bs.b_spline(u,pz,p)



#plt.plot(cx,cy)
#plt.plot(px,py)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(cx, cy, cz, 'gray')



