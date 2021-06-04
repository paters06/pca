# Python libraries
import numpy as np
import numpy.linalg

# Local project
import src.matrixEquationSolver as matEqnSol

################ EXPLICIT TIME INTEGRATION SCHEME ####################

def explicitScheme(M,K,F,u0,dt,T,enforceddof,enforcedvalues):
    print("EXPLICIT TIME INTEGRATION SCHEME")
    numsteps = int(T/dt + 1)
    tsteps = np.linspace(0,T,numsteps)
    numdof = F.shape[0]
    u_n = np.zeros((numdof,numsteps))

    eigvalues,eigenvectors = np.linalg.eig(K)
    maxEig = np.amax(eigvalues)
    dt_cr = 2/maxEig
    print("Maximum eigenvalue: {0}".format(maxEig))
    print("Maximum time step for explicit scheme: {0} seconds".format(dt_cr))

    if dt < dt_cr:
        Mmod,Kmod,Fmod,totalDofs = matEqnSol.dirichletBCEnforcement_Modified(M,K,F,enforceddof,enforcedvalues)

        # Initial Conditions
        print("STARTING TIME ITERATIONS")
        print("NUMBER OF ITERATIONS: {0}".format(numsteps))
        invM = np.linalg.inv(M)
        A = M - dt*Kmod
        B = dt*Fmod
        # print(invM)
        # print(B)

        u0[enforceddof] = np.reshape(enforcedvalues,(len(enforcedvalues),1))
        u_n[:,0,None] = u0

        # Following timesteps
        for i in range(1,numsteps):
            u_i = u_n[:,i-1,None]
            # print(u_i.T)
            Q_n = A@u_i + B
            # print(Q_n.T)
            u_n[:,i,None] = np.reshape(invM@Q_n,(-1,1))
        # End for loop
        print("TIME INTEGRATION FINISHED")
    else:
        print("Time step too big. The time integration scheme is not stable")
    # End if dt_cr

    return u_n

def implicitScheme(M,K,F,u0,dt,T,enforceddof,enforcedvalues):
    print("IMPLICIT TIME INTEGRATION SCHEME")
    numsteps = int(T/dt + 1)
    tsteps = np.linspace(0,T,numsteps)
    numdof = F.shape[0]
    u_n = np.zeros((numdof,numsteps))

    Mmod,Kmod,Fmod,totalDofs = matEqnSol.dirichletBCEnforcement_Modified(M,K,F,enforceddof,enforcedvalues)

    # print(np.vstack((enforceddof,enforcedvalues)))
    u0[enforceddof] = np.reshape(enforcedvalues,(len(enforcedvalues),1))

    # Initial Conditions
    print("STARTING TIME ITERATIONS")
    print("NUMBER OF ITERATIONS: {0}".format(numsteps))
    E = Mmod + dt*Kmod
    D = Mmod
    invE = np.linalg.inv(E)

    u_n[:,0,None] = u0

    # Following timesteps
    for i in range(1,numsteps):
        u_i = u_n[:,i-1,None]
        # print(u_i.T)
        Q_n = D@u_i + dt*Fmod
        # print(Q_n.T)
        u_n[:,i,None] = np.reshape(invE@Q_n,(-1,1))
    # End for loop
    print("TIME INTEGRATION FINISHED")

    return u_n
