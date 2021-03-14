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
    
    Ra = 0.5
    Rb = 1.0
    crack = 0.1
    
    # km = cos(pi/4)
    km = 0.5*np.sqrt(2)
    
    # Complete list of control points
    fullP = np.array([[0,Ra],[Ra,Ra],[Ra,0],
                      [0,Ra+crack],[Ra+crack,Ra+crack],[Ra+crack,0],
                      [0,Rb],[Rb,Rb],[Rb,0]])

    fullw = np.array([[1],[km],[1],[1],[km],[1],[1],[km],[1]])

    idcontrolpoints = [[0,1,2,3,4,5],[3,4,5,6,7,8]]
    localcontrolpoints = [fullP[idcontrolpoints[0],:],fullP[idcontrolpoints[1],:]]

    multip = [2,2]
    multiq = [1,1]

    multiU = [np.array([0,0,0,1,1,1]),np.array([0,0,0,1,1,1])]
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
    reflist = [['h','h'],['h','h']]
    dirlist = [['U','V'],['U','V']]
#    reflist = [['k','k','h','h'],['k','k','h','h']]
#    dirlist = [['U','V','U','V'],['U','V','U','V']]

    if localRefinement == 'Y':
        multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints = \
        srfn.localPatchRefinement(patchesToRefine,reflist,dirlist,multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints)

#    for pi in fullP:
#        print(pi)

    #disp_i = [id patch,startpt,endpt,restriction,value]
    displacementConditions = [[0,[Ra,0.0],[Ra+crack,0.0],"S",0.0],[1,[Ra+crack,0.0],[Rb,0.0],"S",0.0],
                              [1,[0.0,Ra+crack],[0.0,Rb],"S",0.0]]
    #neumann_i = [id patch,startpt,endpt,type_load,value]
    neumannConditions = [[0,[0.0,0.0],[1.0,0.0],"normal",tv]]
    
    fullSurfacePreprocessing,boundaryPreprocessing,dirichletBCList = \
    multipatchpre2D.problemPreprocessing(multiU,multiV,multip,multiq,fullP,idcontrolpoints,displacementConditions,neumannConditions)
    
#    multipatchpre2D.plotMultipatchGeometry(multiU,multiV,multip,multiq,fullP,fullw,\
#                                           idcontrolpoints,dirichletBCList,boundaryPreprocessing)
    
#    K,F = linElastStat.assemblyMultipatchWeakForm(multiU,multiV,fullw,multip,multiq,fullP,idcontrolpoints, \
#          fullSurfacePreprocessing,numericalquadrature,materialProperties,boundaryPreprocessing)

#    Kred,Fred,removedDofs,totalDofs = matEqnSol.boundaryConditionsEnforcement(K,F,dirichletBCList)

#    dtotal,D = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)

#    multipatchpost2D.postProcessing(multiU,multiV,multip,multiq,fullP,D,fullw,dtotal, \
#                                    idcontrolpoints,fullSurfacePreprocessing,materialProperties)

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
