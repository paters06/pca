# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

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
import src.plottingScripts as plts
import src.nurbs as rbs
import src.preprocessor2D as pre2D
import src.linearElastoStaticsSolver as linElastStat
import src.matrixEquationSolver as matEqnSol
import src.postprocessor2D as post2D
import src.surfaceRefinements as srfn
import src.debugScripts as dbg_scrpt

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
    tv = -1 #Pa
    # uDirichlet = [1,4,7]
    # uAxis = [1,0,1,0,1,0]

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0,0],
                  	 [0.5*L,0],
                  	 [L,0],
                  	 [0,0.5*H],
                  	 [0.5*L,0.5*H],
                  	 [L,0.5*H],
                  	 [0,H],
                  	 [0.5*L,H],
                  	 [L,H]])

    winit = np.array([[1],[1],[1],[1],[1],[1],[1],[1],[1]])

    #Isogeometric routines
    Uinit = np.array([0,0,0,1,1,1])
    Vinit = np.array([0,0,0,1,1,1])
    # Uinit = np.array([0,0,0.5,1,1])
    # Vinit = np.array([0,0,0.5,1,1])

    pinit = 2
    qinit = 2

    geomsurface = rbs.NURBSSurface(Uinit,Vinit,pinit,qinit,Pinit,winit)

    doRefinement = 'Y'

    if doRefinement == 'Y':
        reflist = ['h','h','h','h']
        dirlist = ['U','V','U','V']
        srfn.surfaceRefinement(geomsurface,reflist,dirlist)

    displacementConditions = [[[0.0,0.0],[0.0,H],"C",0.0]]
    neumannConditions = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList = \
    pre2D.problemPreprocessing(geomsurface,displacementConditions,neumannConditions)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # pre2D.plotGeometry(geomsurface,dirichletBCList,neumannConditions,boundaryPreprocessing)

    # K,F = linElastStat.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
    #                                     materialProperties,boundaryPreprocessing,neumannConditions)

    # Kred,Fred,removedDofs,totalDofs = matEqnSol.boundaryConditionsEnforcement(K,F,dirichletBCList)

    # dtotal,D = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)

    # post2D.postProcessing(geomsurface,D,dtotal,surfacePreprocessing,materialProperties)

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
