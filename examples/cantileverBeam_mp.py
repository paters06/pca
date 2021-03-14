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
import src.plottingScripts as plts
import src.multipatchPreprocessor2D as multipatchpre2D
import src.linearElastoStaticsSolver as linElastStat
import src.matrixEquationSolver as matEqnSol
import src.multipatchPostprocessor2D as multipatchpost2D
import src.surfaceRefinements as srfn
import src.debugScripts as dbg_scrpt

####################################################
################## MAIN PROBLEM ####################
####################################################

def mainProgram():
    #Data
    
    fullP = np.array([[0,0],[0.4,0],[0.0,0.2],[0.4,0.2],
                  [1.0,0.0],[1.0,0.2]])

    fullw = np.array([[1],[1],[1],[1],[1],[1]])

    idcontrolpoints = [[0,1,2,3],[1,4,3,5]]
    localcontrolpoints = [fullP[idcontrolpoints[0],:],fullP[idcontrolpoints[1],:]]

    multip = [1,1]
    multiq = [1,1]

#    multiU = [np.array([0,0,0.5,0.5]),np.array([0.5,0.5,1,1])]
    multiU = [np.array([0,0,1,1]),np.array([0,0,1,1])]
    multiV = [np.array([0,0,1,1]),np.array([0,0,1,1])]
    
    E = 2e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = -1 #Pa

    numGaussPoints = 4
    numericalquadrature = multipatchpre2D.numericalIntegrationPreprocessing(numGaussPoints)

    localRefinement = 'Y'
    patchesToRefine = [0,1]
#    reflist = [['h','h','h','h'],['h','h']]
#    dirlist = [['U','V','U','V'],['U','V']]
#    reflist = [['k','k','h','h'],['k','k']]
#    dirlist = [['U','V','U','V'],['U','V']]

    if localRefinement == 'Y':
        multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints = \
        srfn.localPatchRefinement(patchesToRefine,reflist,dirlist,multiU,multiV,multip,multiq,fullP,fullw,\
        idcontrolpoints,localcontrolpoints)

    for mU in multiU:
        print(mU)

    #disp_i = [id patch,[startpt,endpt],restriction,value]
    displacementConditions = [[0,[0.0,0.0],[0.0,0.2],"C",0.0]]
    #neumann_i = [id patch,[startpt,endpt],type_load,value]
    neumannConditions = [[1,[1.0,0.0],[1.0,1.0],"tangent",tv]]
    
    fullParametricNodes,fullNodesInElement,fullSurfacePreprocessing = multipatchpre2D.parametricGrid(multiU,multiV,multip,multiq)
    boundaryPreprocessing = multipatchpre2D.loadPreprocessing(fullParametricNodes,fullNodesInElement,neumannConditions,\
                            multiU,multiV,multip,multiq)
    dirichletBCList = multipatchpre2D.dirichletBCPreprocessingOnFaces(fullP,idcontrolpoints,displacementConditions)
    
#    multipatchpre2D.plotMultipatchGeometry(multiU,multiV,multip,multiq,fullP,fullw,\
#                                           idcontrolpoints,dirichletBCList,boundaryPreprocessing)
    
    K,F = linElastStat.assemblyMultipatchWeakForm(multiU,multiV,fullw,multip,multiq,fullP,idcontrolpoints, \
          fullSurfacePreprocessing,numericalquadrature,materialProperties,boundaryPreprocessing)

    Kred,Fred,removedDofs,totalDofs = matEqnSol.boundaryConditionsEnforcement(K,F,dirichletBCList)

    dtotal,D = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)

    multipatchpost2D.postProcessing(multiU,multiV,multip,multiq,fullP,D,fullw,dtotal, \
                                    idcontrolpoints,fullSurfacePreprocessing,materialProperties)

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
