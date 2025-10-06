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
import src.preprocessorPlate as preP
import src.phenomenon.reissnerMidlinPlateSolver as reisMidSol
import src.matrixEquationSolver as matEqnSol
import src.postprocessorPlate as postP
import src.spline_functions.surfaceRefinements as srfn
import src.debug.debugScripts as dbg_scrpt

####################################################
################## MAIN PROBLEM ####################
####################################################

def mainProgram():
    #Data
    phenomenon = "Elasticity"
    A = 1650
    B = 1650
    t = 30
    E = 200e3 #MPa
    nu = 0.3
    rho = 0.0 #kg/m3
    load = -1 #Pa
    materialProperties = [E,nu,rho,load,t]

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0,0],[0.5*A,0],[A,0],[0,0.5*B],
                  	 [0.5*A,0.5*B],[A,0.5*B],[0,B],[0.5*A,B],[A,B]])

    winit = np.ones((Pinit.shape[0],1))

    #Isogeometric routines
    Uinit = np.array([0,0,0,1,1,1])
    Vinit = np.array([0,0,0,1,1,1])

    pinit = 2
    qinit = 2

    geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)

    doRefinement = 'Y'

    if doRefinement == 'Y':
        # srfn.surfaceRefinement(geomsurface,1,'p','U')
        # srfn.surfaceRefinement(geomsurface,1,'p','V')
        srfn.surfaceRefinement(geomsurface,2,'h','U')
        srfn.surfaceRefinement(geomsurface,2,'h','V')
    # End function

    side_1 = [[0.0,0.0],[0.0,1.0],0.0,"C"]
    side_2 = [[0.0,1.0],[1.0,1.0],0.0,"C"]
    side_3 = [[1.0,1.0],[1.0,0.0],0.0,"C"]
    side_4 = [[0.0,0.0],[1.0,0.0],0.0,"C"]
    # dirichletConditionsData = [side_1]
    dirichletConditionsData = [side_1,side_2,side_3,side_4]
    neumannConditionsData = None

    surfacePreprocessing,boundaryPreprocessing,enforcedCtrlPts,enforcedDOF,enforcedValues = \
    preP.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = preP.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)
    
    # preP.plotGeometry(phenomenon,geomsurface,enforcedCtrlPts,boundaryPreprocessing)

    K,F,M = reisMidSol.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)

    Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforcedDOF,enforcedValues)

    dSolution = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    # postP.postProcessing(phenomenon,geomsurface,surfacePreprocessing,dSolution,materialProperties)

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
