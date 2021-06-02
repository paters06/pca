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
        Kmod,Fmod,totalDofs = matEqnSol.dirichletBCEnforcement_Modified(K,F,enforceddof,enforcedvalues)

        # Initial Conditions
        print("Iteration #0")
        invM = np.linalg.inv(M)
        A = M - dt*Kmod
        B = dt*Fmod
        # print(invM)
        # print(B)

        u_n[:,0,None] = u0

        # Following timesteps
        for i in range(1,numsteps):
            u_i = u_n[:,i-1,None]
            # print(u_i.T)
            Q_n = A@u_i + B
            # print(Q_n.T)
            u_n[:,i,None] = np.reshape(invM@Q_n,(-1,1))
            if i%100 == 0:
                print("Iteration #{0}".format(i))
                # print(u_i.T)
                # print(u_n[:,i,None].T)
        # End for loop
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

    Kmod,Fmod,totalDofs = matEqnSol.dirichletBCEnforcement_Modified(K,F,enforceddof,enforcedvalues)

    # Initial Conditions
    print("Iteration #0")
    E = M + dt*Kmod
    D = M
    invE = np.linalg.inv(E)

    u_n[:,0,None] = u0

    # Following timesteps
    for i in range(1,numsteps):
        u_i = u_n[:,i-1,None]
        # print(u_i.T)
        Q_n = D@u_i + dt*Fmod
        # print(Q_n.T)
        u_n[:,i,None] = np.reshape(invE@Q_n,(-1,1))
        if i%10 == 0:
            print("Iteration #{0}".format(i))
            # print(u_i.T)
            # print(u_n[:,i,None].T)
    # End for loop

    return u_n
