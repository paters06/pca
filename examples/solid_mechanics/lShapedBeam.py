# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

#######################################################################
# DO NOT REMOVE THIS SEGMENT
import os
import sys

# Get current working directory
dir1 = os.getcwd()
# Insert .. command to go to the upper directory
dir2 = dir1 + '/..'
# Change directory
os.chdir(dir2)
# Get the new current working directory
dir3 = os.getcwd()
# Setting the package directory path for the modules execution
sys.path.append(dir3)
#######################################################################

# Local project
import src.plotting.plottingScripts as plts
import src.preprocessor2D as pre2D
import src.phenomenon.linearElastoStaticsSolver as linElastStat
import src.matrixEquationSolver as matEqnSol
import src.postprocessor2D as post2D
import src.spline_functions.surfaceRefinements as srfn
import src.debug.debugScripts as dbg_scrpt

####################################################
################## MAIN PROBLEM ####################
####################################################

def mainProgram():
    #Data
    L = 1.0
    H = 0.2
    E = 2e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = 1 #Pa
    # uDirichlet = [1,4,7]
    # uAxis = [1,0,1,0,1,0]

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0.2,0.0],
                      [0.2,0.4],
                      [0.6,0.4],
                      [0.0,0.0],
                      [0.0,0.6],
                      [0.6,0.6]])

    winit = np.array([[1],[1],[1],[1],[1],[1]])

    #Isogeometric routines
    Uinit = np.array([0,0,0.5,1,1])
    Vinit = np.array([0,0,1,1])

    pinit = 1
    qinit = 1

    doRefinement = 'Y'

    if doRefinement == 'Y':
#        reflist = ['h','h','h','h','p','p']
#        dirlist = ['V','U','V','U','U','V']
#        reflist = ['h','h','h','h','h','h','h','h','h','h']
#        dirlist = ['V','U','U','V','U','V','U','V','U','V']
        Uinp,Vinp,pinp,qinp,Pinp,winp = srfn.surfaceRefinement(reflist,dirlist,Uinit,Vinit,pinit,qinit,Pinit,winit)
    else:
        Uinp = Uinit
        Vinp = Vinit
        pinp = pinit
        qinp = qinit
        Pinp = Pinit
        winp = winit

#    for pi in Pinp:
#        print(pi)

    displacementConditions = [[[0.0,0.0],[0.2,0.0],"C",0.0]]
    neumannConditions = [[[1.0,0.0],[1.0,1.0],"normal",tv]]

    parametricNodes,nodesInElement,surfacePreprocessing = pre2D.parametricGrid(Uinp,Vinp,pinp,qinp)
    loadElements,boundaryPreprocessing = pre2D.loadPreprocessing(parametricNodes,nodesInElement,neumannConditions,Uinp,Vinp,pinp,qinp)
    dirichletBCList = pre2D.dirichletBCPreprocessingOnFaces(Pinp,displacementConditions)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

#    pre2D.plotGeometry(Uinp,Vinp,pinp,qinp,Pinp,winp,dirichletBCList,neumannConditions,boundaryPreprocessing)

    K,F = linElastStat.assemblyWeakForm(Uinp,Vinp,winp,pinp,qinp,Pinp,surfacePreprocessing,numericalquadrature,materialProperties,boundaryPreprocessing,neumannConditions)

    Kred,Fred,removedDofs,totalDofs = matEqnSol.boundaryConditionsEnforcement(K,F,dirichletBCList)

    dtotal,D = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)

    post2D.postProcessing(Uinp,Vinp,pinp,qinp,Pinp,D,winp,dtotal,surfacePreprocessing,materialProperties)

mainProgram()

#import cProfile
#import pstats
#profiler = cProfile.Profile()
#profiler.enable()
#mainProgram()
#profiler.disable()
## stats = pstats.Stats(profiler).sort_stats('ncalls')
#stats = pstats.Stats(profiler).sort_stats('tottime')
#stats.print_stats()
