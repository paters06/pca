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
    L = 4
    H = 1.0
    kappa = 237.0
    rho = 903.0
    source = 0.0
    materialProperties = [kappa,rho,source]
    flux = 0.0

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0.0,0.0],[L,0.0],[0.0,H],[L,H]])

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

    # dirichletConditionsData = [[[1.0,0.0],[1.0,1.0],85.0]]
    # dirichletConditionsData = [[[0.0,0.0],[0.0,1.0],0.0],[[1.0,0.0],[1.0,1.0],85.0]]
    dirichletConditionsData = [[[0.0,0.0],[0.0,1.0],25.0],[[1.0,0.0],[1.0,1.0],85.0]]
    # dirichletConditionsData = [[[1.0,0.0],[1.0,1.0],0.0],[[0.0,0.0],[1.0,0.0],0.0],
    #                            [[0.0,1.0],[1.0,1.0],0.0],[[0.0,0.0],[0.0,1.0],100.0]]
    # neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"tangent",flux],[[0.0,1.0],[1.0,1.0],"tangent",flux]]
    neumannConditionsData = None

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F,M = diffSol.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)


    # Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement_Reduced(K,F,enforcedDOF,enforcedValues)
    # dtotal,D = matEqnSol.solveReducedMatrixEquations(phenomenon,Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    # Mmod,Kmod,Fmod,totalDofs = matEqnSol.dirichletBCEnforcement_Modified(M,K,F,enforcedDOF,enforcedValues)
    # dtotal,D = matEqnSol.solveModifiedMatrixEquations(phenomenon,Kmod,Fmod,totalDofs)

    # post2D.postProcessing(phenomenon,geomsurface,D,dtotal,surfacePreprocessing,materialProperties)

    # T = 5.0
    # dt = 0.01
    # uInitial = np.full((F.shape[0],1),25.0)
    # # uTransient = timeIntg.explicitScheme(M,K,F,uInitial,dt,T,enforcedDOF,enforcedValues)
    # uTransient = timeIntg.implicitScheme(M,K,F,uInitial,dt,T,enforcedDOF,enforcedValues)
    # post2D.plotTransientField(phenomenon,geomsurface,surfacePreprocessing,materialProperties,uTransient,T,dt,True)

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
