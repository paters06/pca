#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 19:23:10 2020

@author: hernando
Este codigo genera superficies B-splines
"""

import numpy as np
import NURBS as bs
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


#px=np.array([[1,1,1,1],[2,2,2,2]])
#py=np.array([[1,1,1,1],[2,2,2,2]])

px = np.linspace(-6, 6, 5)
py = np.linspace(-6, 6, 5)

px, py = np.meshgrid(px, py)

pz=px*px*px-py*py#np.random.randn(30, 30)

p=1
q=1

#Ui=knot(len(px)-1,p)
#Uj=knot(len(px)-1,q)

u=bs.ugen(0.0,1.0,0.01)
v=bs.ugen(0.0,1.0,0.01)

cx=bs.b_spline2d(u,v,px,p,q)
cy=bs.b_spline2d(u,v,py,p,q)
cz=bs.b_spline2d(u,v,pz,p,q)


fig = plt.figure()
ax = plt.axes(projection='3d')
#ax.plot3D(px, py, pz, 'gray')
#ax.contour3D(px, py, pz, 50, cmap='binary')
ax.plot_wireframe(px, py, pz, color='red')
ax.contour3D(cx, cy, cz, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z');