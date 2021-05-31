# Python libraries
import numpy as np
import numpy.linalg

def explicitScheme(M,K,F,u,dt,T):
    numsteps = T/dt + 1
    tsteps = np.linspace(0,T,numsteps)
    u_n = np.zeros((F.shpae[0],numsteps))
    # Initial Conditions
    invM = np.linalg.inv(M)
    A = M - dt*K
    B = dt*F
    # 
    # u_n[:,0] = u

    # Following timesteps
    for i in range(1,numsteps+1):
        Q_n = A@u_n[:,i-1] + B
        u_n[:,i] = invM@Q_n
    return u_n
