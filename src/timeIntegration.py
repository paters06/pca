# Python libraries
import numpy as np
import math
import numpy.linalg

# Local project
import src.matrixEquationSolver as matEqnSol

###################################################
################ PARABOLIC PDE ####################
###################################################

################ EXPLICIT TIME INTEGRATION SCHEME ####################

def parabolicExplicitScheme(M,K,F,u0,dt,T,enforceddof,enforcedvalues):
    print("EXPLICIT TIME INTEGRATION SCHEME")
    numsteps = int(T/dt + 1)
    numdof = F.shape[0]
    u_n = np.zeros((numdof,numsteps))

    eigvalues,eigenvectors = np.linalg.eig(K)
    maxEig = np.amax(eigvalues)
    dt_cr = 2/maxEig
    print("Maximum eigenvalue: {0}".format(maxEig))
    print("Maximum time step for explicit scheme: {0} seconds".format(dt_cr))

    if dt < dt_cr:
        Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforceddof,enforcedvalues)
        
        reduceddofs = np.setdiff1d(totalDofs,enforceddof)
        # u0[enforceddof] = np.reshape(enforcedvalues,(len(enforcedvalues),1))
        # u_n[:,0,None] = u0
        ured_n = np.zeros((len(reduceddofs),numsteps))
        ured_n[:,0,None] = u0[reduceddofs]

        # Initial Conditions
        print("STARTING TIME ITERATIONS")
        print("NUMBER OF ITERATIONS: {0}".format(numsteps))
        invM = np.linalg.inv(Mred)
        A = Mred - dt*Kred
        B = dt*Fred
        # print(invM)
        # print(B)
        
        enforcedvalues = np.reshape(enforcedvalues,(len(enforcedvalues),1))
        u_n[enforceddof,:] = enforcedvalues

        # Following timesteps
        for i in range(1,numsteps):
            ured_i = ured_n[:,i-1,None]
            # print(u_i.T)
            Q_n = A@ured_i + B
            # print(Q_n.T)
            ured_n[:,i,None] = np.reshape(invM@Q_n,(-1,1))
        # End for loop
        print("TIME INTEGRATION FINISHED")

        u_n[reduceddofs,:] = ured_n
        # print(u_n[:,-1].T)
        # print(ured_n[:,-1].T)
    else:
        print("Time step too big. The time integration scheme is not stable")
    # End if dt_cr

    return u_n
# End function

def parabolicImplicitScheme(M,K,F,u0,dt,T,enforceddof,enforcedvalues):
    print("IMPLICIT TIME INTEGRATION SCHEME")
    numsteps = int(T/dt + 1)
    numdof = F.shape[0]
    u_n = np.zeros((numdof,numsteps))

    Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforceddof,enforcedvalues)

    reduceddofs = np.setdiff1d(totalDofs,enforceddof)
    ured_n = np.zeros((len(reduceddofs),numsteps))
    ured_n[:,0,None] = u0[reduceddofs]

    # print(np.vstack((enforceddof,enforcedvalues)))
    # u0[enforceddof] = np.reshape(enforcedvalues,(len(enforcedvalues),1))

    # Initial Conditions
    print("STARTING TIME ITERATIONS")
    print("NUMBER OF ITERATIONS: {0}".format(numsteps))
    E = Mred + dt*Kred
    D = Mred
    invE = np.linalg.inv(E)

    # u_n[:,0,None] = u0

    # Following timesteps
    for i in range(1,numsteps):
        ured_i = ured_n[:,i-1,None]
        # print(u_i.T)
        Q_n = D@ured_i + dt*Fred
        # print(Q_n.T)
        ured_n[:,i,None] = np.reshape(invE@Q_n,(-1,1))
    # End for loop
    print("TIME INTEGRATION FINISHED")
    u_n[reduceddofs,:] = ured_n
    return u_n
# End function

####################################################
################ HYPERBOLIC PDE ####################
####################################################

def hyperbolicExplicitScheme(M,K,F,d0,v0,dt,T,enforceddof,enforcedvalues):
    # from scipy.sparse.linalg import eigs,eigsh
    from scipy.linalg import eigh
    
    print("EXPLICIT TIME INTEGRATION SCHEME")
    beta = 0.0/2.0
    gamma = 1.0/2.0

    numsteps = int(T/dt + 1)
    numdof = F.shape[0]
    d_n = np.zeros((numdof,numsteps))

    # eigvalues = np.linalg.eigvals(K)
    Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforceddof,enforcedvalues)
    
    # eigvalues = np.linalg.eigvals(Kred)
    # eigvalues,eigvectors = eigsh(Kred,M=Mred,k=10,which='SM')
    eigvalues,eigvectors = eigh(Kred,b=Mred)
    # print(eigvalues)
    print(np.sqrt(eigvalues))
    # print(np.sqrt(np.sort(eigvalues)))
    maxEig = np.amax(eigvalues)
    dt_cr = 1e15
    # if abs(beta -gamma) > 1e-5:
    #     dt_cr = 1.0/math.sqrt(0.5*maxEig*(beta - gamma))
    # else:
    #     dt_cr = 1e15
    # End if
    print("Maximum eigenvalue: {0}".format(maxEig))
    print("Maximum time step for explicit scheme: {0} seconds".format(dt_cr))

    if dt < dt_cr:
        Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforceddof,enforcedvalues)

        reduceddofs = np.setdiff1d(totalDofs,enforceddof)
        dred_n = np.zeros((len(reduceddofs),numsteps))
        dred_n[:,0,None] = d0[reduceddofs]

        Fred_tot = Fred - Kred@dred_n[:,0,None]
        dn = dred_n[:,0,None]
        vn = v0[reduceddofs]
        an = np.linalg.solve(Mred,Fred_tot)

        numsteps = int(T/dt + 1)

        # Initial Conditions
        print("STARTING TIME ITERATIONS")
        print("NUMBER OF ITERATIONS: {0}".format(numsteps))
        for i in range(1,numsteps):
            # Calculation of predictors
            dn_predictor = dn + dt*vn + 0.5*dt**2*(1.0 - 2*beta)*an
            vn_predictor = vn + (1.0 - gamma)*dt*an
            # Calculation of updated acceleration
            Mred_tot = Mred + beta*(dt**2)*Kred
            Fred_tot = Fred - Kred@dn_predictor
            an_predictor = np.linalg.solve(Mred_tot,Fred_tot)

            # Calculation of correctors
            dn_corrector = dn_predictor + beta*(beta**2)*an_predictor
            vn_corrector = vn_predictor + gamma*dt*an_predictor

            # Update of dn and vn
            dn = dn_corrector
            vn = vn_corrector

            dred_n[:,i,None] = dn
        # End for loop
        print("TIME INTEGRATION FINISHED")

        # Assembly of the reduced solution
        d_n[reduceddofs,:] = dred_n
    else:
        print("Time step too big. The time integration scheme is not stable")
    # End if dt_cr

    return d_n
# End function

def hyperbolicImplicitScheme(M,K,F,d0,v0,dt,T,enforceddof,enforcedvalues):
    print("IMPLICIT TIME INTEGRATION SCHEME")
    numsteps = int(T/dt + 1)
    numdof = F.shape[0]
    d_n = np.zeros((numdof,numsteps))

    Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforceddof,enforcedvalues)

    reduceddofs = np.setdiff1d(totalDofs,enforceddof)
    dred_n = np.zeros((len(reduceddofs),numsteps))
    dred_n[:,0,None] = d0[reduceddofs]

    Fred_tot = Fred - Kred@dred_n[:,0,None]
    dn = dred_n[:,0,None]
    vn = v0[reduceddofs]
    an = np.linalg.solve(Mred,Fred_tot)

    beta = 1.0/4
    gamma = 1.0/2

    numsteps = int(T/dt + 1)
    tsteps = np.linspace(0,T,numsteps)

    # Initial Conditions
    print("STARTING TIME ITERATIONS")
    print("NUMBER OF ITERATIONS: {0}".format(numsteps))
    for i in range(1,numsteps):
        # Calculation of predictors
        dn_predictor = dn + dt*vn + 0.5*(dt**2)*(1.0 - 2*beta)*an
        vn_predictor = vn + (1.0 - gamma)*dt*an
        # Calculation of updated acceleration
        Mred_tot = Mred + beta*(dt**2)*Kred
        Fred_tot = Fred - Kred@dn_predictor
        an_predictor = np.linalg.solve(Mred_tot,Fred_tot)

        # Calculation of correctors
        dn_corrector = dn_predictor + beta*(dt**2)*an_predictor
        vn_corrector = vn_predictor + gamma*dt*an_predictor

        # Update of dn and vn
        dn = dn_corrector
        vn = vn_corrector

        dred_n[:,i,None] = dn
    # End for loop
    print("TIME INTEGRATION FINISHED")

    # Assembly of the reduced solution
    d_n[reduceddofs,:] = dred_n
    return d_n
# End function