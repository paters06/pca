#Primera curva B-spline
import numpy as np
import matplotlib.pyplot as plt
#exec('F_BASE.py')
import F_BASE
import knot
Px=np.array([0,1,4,6,8,10])
Py=np.array([-1,0,1,0,-1,0])

print (Px)
print (Py)
#Grado
p=2

rango=len(Px)
nn=rango-1
U=knot.knot(nn,p)

rango=np.arange(0,1,0.01)

CX=rango*0
CY=rango*0
i=0

for u in rango: 
 N0=F_BASE.F_BASE(U,u,p)
 #N0=np.array(N0)
 CX[i]=N0@Px
 CY[i]=N0@Py
 i=i+1
 
plt.plot(CX,CY) 
#plt.plot(Px,Py) 
 