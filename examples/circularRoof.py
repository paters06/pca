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
import src.nurbs as rbs
import src.preprocessor2D as pre2D
import src.linearElastoStaticsSolver as linElastStat
import src.matrixEquationSolver as matEqnSol
import src.postprocessor2D as post2D
import src.surfaceRefinements as srfn
import src.debugScripts as dbg_scrpt

####################################################
################### MAIN PROBLEM ###################
####################################################

def mainProgram():
    #Data
    phenomenon = "Elasticity"
    E = 1e5 #Pa
    nu = 0.31
    rho = 10.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = 10 #Pa
    Rmax = 1.0
    Rmin = 0.5
    numGaussPoints = 4

    Pinit = np.array([[Rmin,0],[Rmin,Rmin],[0,Rmin],[-Rmin,Rmin],[-Rmin,0],
                      [Rmax,0],[Rmax,Rmax],[0,Rmax],[-Rmax,Rmax],[-Rmax,0]])
    winit = np.array([[1],[0.5*np.sqrt(2)],[1],[0.5*np.sqrt(2)],[1],
                      [1],[0.5*np.sqrt(2)],[1],[0.5*np.sqrt(2)],[1]])

    #Isogeometric routines
    Uinit = np.array([0,0,0,0.5,0.5,1,1,1])
    Vinit = np.array([0,0,1,1])

    pinit = 2
    qinit = 1

    geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)

    doRefinement = 'Y'

    if doRefinement == 'Y':
        # srfn.surfaceRefinement(geomsurface,1,'p','U')
        # srfn.surfaceRefinement(geomsurface,1,'p','V')
        srfn.surfaceRefinement(geomsurface,1,'h','U')
        srfn.surfaceRefinement(geomsurface,1,'h','V')

    dirichletConditionsData = [[[0.0,0.0],[0.0,1.0],0.0,"C"],[[1.0,0.0],[1.0,1.0],0.0,"C"]]
    # neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"normal",tv]]
    neumannConditionsData = None

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # print(geomsurface.P)
    # print(np.vstack((enforcedDOF,enforcedValues)))

    pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F,M = linElastStat.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)

    Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement_Reduced(K,F,enforcedDOF,enforcedValues)
    dtotal,D = matEqnSol.solveReducedMatrixEquations(phenomenon,Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    # post2D.postProcessing(phenomenon,geomsurface,D,dtotal,surfacePreprocessing,materialProperties)

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
