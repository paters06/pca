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
import src.phenomenon.diffusionSolver as diffSol
import src.phenomenon.timeIntegration as timeIntg
import src.matrixEquationSolver as matEqnSol
import src.postprocessor2D as post2D
import src.spline_functions.surfaceRefinements as srfn
import src.debug.debugScripts as dbg_scrpt

####################################################
################## MAIN PROBLEM ####################
####################################################

def mainProgram():
    #Data
    phenomenon = "Heat"
    Rmax = 1.0
    Rmin = 0.7
    kappa = 1 #Pa
    rho = 1.0
    source = 1.0 #kg/m3
    materialProperties = [kappa,rho,source]
    flux = 0.0 #Pa

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0.0,0.0],[1.0,0.0],[0.0,1.0],[1.0,1.0]])

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

    # dirichletConditionsData = [[[1.0,0.0],[1.0,1.0],0.0],[[0.0,1.0],[1.0,1.0],0.0]]
    dirichletConditionsData = [[[1.0,0.0],[1.0,1.0],0.0],[[0.0,0.0],[1.0,0.0],0.0],
                               [[0.0,1.0],[1.0,1.0],0.0],[[0.0,0.0],[0.0,1.0],100.0]]
    # neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"tangent",flux],[[0.0,1.0],[1.0,1.0],"tangent",flux]]
    neumannConditionsData = None

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F,M = diffSol.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)

    T = 1.0
    dt = 0.005
    uInitial = np.full((F.shape[0],1),0.0)
    uTransient = timeIntg.parabolicExplicitScheme(M,K,F,uInitial,dt,T,enforcedDOF,enforcedValues)
    # uTransient = timeIntg.parabolicImplicitScheme(M,K,F,uInitial,dt,T,enforcedDOF,enforcedValues)

    # print(np.hstack((geomsurface.P,uTransient[:,5,None])))

    # Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforcedDOF,enforcedValues)
    # dSolution = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    # print(np.hstack((geomsurface.P,uTransient[:,-1,None],dSolution)))

    # post2D.postProcessing(phenomenon,geomsurface,surfacePreprocessing,dSolution,materialProperties)
    post2D.plotTransientField(geomsurface,surfacePreprocessing,materialProperties,uTransient,T,dt,False)

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
