# Python libraries
import numpy as np
import time
import math

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
import src.phenomenon.schrodingerEquationSolver as schroSol
import src.matrixEquationSolver as matEqnSol
import src.postprocessor2D as post2D
import src.spline_functions.surfaceRefinements as srfn
import src.debug.debugScripts as dbg_scrpt

####################################################
################## MAIN PROBLEM ####################
####################################################

def mainProgram():
    #Data
    phenomenon = "Schrodinger"
    hbar = 1.0
    m = 1.0
    L = 5.0

    x = np.linspace(0,L)
    y = np.linspace(0,L)
    X,Y = np.meshgrid(x,y)

    def potentialField(x):
        # return 0*x[0]
        return -(1.6e-19)**2*1.0/(2*8.854e-12*math.sqrt((x[0] - 2.5)**2 + (x[1] - 2.5)**2))
    # End function

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0.0,0.0],[L,0.0],[0.0,L],[L,L]])
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
        srfn.surfaceRefinement(geomsurface,4,'h','U')
        srfn.surfaceRefinement(geomsurface,4,'h','V')

    dirichletConditionsData = [[[1.0,0.0],[1.0,1.0],0.0],[[0.0,0.0],[1.0,0.0],0.0],
                               [[0.0,1.0],[1.0,1.0],0.0],[[0.0,0.0],[0.0,1.0],0.0]]
    # neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"tangent",flux],[[0.0,1.0],[1.0,1.0],"tangent",flux]]
    neumannConditionsData = None

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F,M = schroSol.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                      potentialField,boundaryPreprocessing)

    Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement_Reduced(M,K,F,enforcedDOF,enforcedValues)
    # dSolution = matEqnSol.solveReducedMatrixEquations(Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    eigenSolution = matEqnSol.solveReducedEigenvalueProblem(Mred,Kred,6,totalDofs,enforcedDOF,enforcedValues)

    # print(np.hstack((geomsurface.P,eigenSolution)))

    post2D.postProcessing(phenomenon,geomsurface,surfacePreprocessing,eigenSolution)
# End function

start = time.time()
mainProgram()
end = time.time()
print("Elapsed time {:.3f} seconds".format(end - start))

# import cProfile
# import pstats
# profiler = cProfile.Profile()
# profiler.enable()
# mainProgram()
# profiler.disable()
# # stats = pstats.Stats(profiler).sort_stats('ncalls')
# stats = pstats.Stats(profiler).sort_stats('tottime')
# stats.print_stats()
