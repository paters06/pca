# Python libraries
import numpy as np

#######################################################################
# DO NOT REMOVE THIS SEGMENT
import os
import sys

# Get current working directory
# Command found in https://note.nkmk.me/en/python-script-file-path/
dir1 = os.path.dirname(os.path.abspath(__file__))
# Insert .. command to go to the upper directory
dir2 = dir1 + '/..'
# Setting the package directory path for the modules execution
sys.path.append(dir2)
#######################################################################

# Local project
import src.nurbs as rbs
import src.preprocessor2D as pre2D
import src.waveEquationSolver as waveEq
import src.timeIntegration as timeIntg
import src.matrixEquationSolver as matEqnSol
import src.postprocessor2D as post2D
import src.surfaceRefinements as srfn
import src.debugScripts as dbg_scrpt

def defineInitialConditions(P):
    icu = np.zeros((P.shape[0],1))
    icv = np.zeros((P.shape[0],1))
    for i in range(P.shape[0]):
        xt = P[i][0]# + 2.0
        yt = P[i][1]# + 1.0
        icu[i][0] = 0.1*(4*xt - xt**2)*(2*yt - yt**2)
    # End for loop
    return icu,icv
# End function

####################################################
################## MAIN PROBLEM ####################
####################################################

def mainProgram():
    #Data
    phenomenon = "Heat"
    a = 4.0
    b = 2.0
    kappa = 12.5 #Pa
    rho = 2.5
    source = 0.0 #kg/m3
    materialProperties = [kappa,rho,source]

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0.0,0.0],[a,0.0],[0.0,b],[a,b]])

    winit = np.ones((Pinit.shape[0],1))

    #Isogeometric routines
    Uinit = np.array([0,0,1,1])
    Vinit = np.array([0,0,1,1])

    pinit = 1
    qinit = 1

    geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)

    doRefinement = 'Y'

    if doRefinement == 'Y':
        # srfn.surfaceRefinement(geomsurface,1,'p','U')
        # srfn.surfaceRefinement(geomsurface,1,'p','V')
        srfn.surfaceRefinement(geomsurface,2,'h','U')
        srfn.surfaceRefinement(geomsurface,2,'h','V')
    # End if

    # dirichletConditionsData = [[[1.0,0.0],[1.0,1.0],0.0],[[0.0,1.0],[1.0,1.0],0.0]]
    dirichletConditionsData = [[[1.0,0.0],[1.0,1.0],0.0],[[0.0,0.0],[1.0,0.0],0.0],
                               [[0.0,1.0],[1.0,1.0],0.0],[[0.0,0.0],[0.0,1.0],0.0]]
    # neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"tangent",flux],[[0.0,1.0],[1.0,1.0],"tangent",flux]]
    neumannConditionsData = None

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F,M = waveEq.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)

    # Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement_Reduced(K,F,enforcedDOF,enforcedValues)
    # dSolution = matEqnSol.solveReducedMatrixEquations(Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    T = 1.0
    dt = 0.025
    # uInitial = np.full((F.shape[0],1),0.0)
    d0,v0 = defineInitialConditions(geomsurface.P)
    # uTransient = timeIntg.hyperbolicExplicitScheme(M,K,F,d0,v0,dt,T,enforcedDOF,enforcedValues)
    uTransient = timeIntg.hyperbolicImplicitScheme(M,K,F,d0,v0,dt,T,enforcedDOF,enforcedValues)
    # print(uTransient)

    # Mmod,Kmod,Fmod,totalDofs = matEqnSol.dirichletBCEnforcement_Reduced(M,K,F,enforcedDOF,enforcedValues)
    # dSolution = matEqnSol.solveReducedMatrixEquations(Kmod,Fmod,totalDofs)

    # print(np.hstack((dtotal,dtotal2)))

    # post2D.postProcessing(phenomenon,geomsurface,surfacePreprocessing,dSolution,materialProperties)
    post2D.plotTransientField(geomsurface,surfacePreprocessing,materialProperties,uTransient,T,dt,False)
# End function

mainProgram()

# import cProfile
# import pstats
# profiler = cProfile.Profile()
# profiler.enable()
# mainProgram()
# profiler.disable()
# # stats = pstats.Stats(profiler).sort_stats('ncalls')
# stats = pstats.Stats(profiler).sort_stats('tottime')
# stats.print_stats()
