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
import src.plottingScripts as plts
import src.nurbs as rbs
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

    fullP = np.array([[0,0],[0.2,0],[0,0.6],[0.2,0.4],
                  [0.6,0.4],[0.6,0.6]])

    fullw = np.array([[1],[1],[1],[1],[1],[1]])

    idcontrolpoints = [[0,1,2,3],[3,4,2,5]]
    localcontrolpoints = [fullP[idcontrolpoints[0],:],fullP[idcontrolpoints[1],:]]

    multip = [1,1]
    multiq = [1,1]

    # multiU = [np.array([0,0,0.5,0.5]),np.array([0.5,0.5,1,1])]
    multiU = [np.array([0,0,1,1]),np.array([0,0,1,1])]
    multiV = [np.array([0,0,1,1]),np.array([0,0,1,1])]

    geomsurface = rbs.MultiPatchNURBSSurface(multiU,multiV,multip,multiq,\
                                             fullP,fullw,idcontrolpoints)

    E = 2e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = 1 #Pa

    numGaussPoints = 4
    numericalquadrature = multipatchpre2D.numericalIntegrationPreprocessing(numGaussPoints)

    localRefinement = 'Y'
    patchesToRefine = [0,1]
    # reflist = [['h','h','p','p','h','h'],['h','h','p','p']]
    # dirlist = [['U','V','U','V','U','V'],['U','V','U','V']]
    reflist = [['k','k','h','h'],['k','k','h','h']]
    dirlist = [['U','V','U','V'],['U','V','U','V']]

    if localRefinement == 'Y':
        localcontrolpoints = srfn.localPatchRefinement(geomsurface,patchesToRefine,reflist,dirlist,localcontrolpoints)

    # for pi in fullP:
    #     print(pi)

    #disp_i = [id patch,[startpt,endpt],restriction,value]
    displacementConditions = [[0,[0.0,0.0],[0.2,0.0],"C",0.0]]
    #neumann_i = [id patch,[startpt,endpt],type_load,value]
    neumannConditions = [[1,[1.0,0.0],[1.0,1.0],"normal",tv]]

    fullSurfacePreprocessing,boundaryPreprocessing,dirichletBCList = \
    multipatchpre2D.problemPreprocessing(geomsurface,displacementConditions,neumannConditions)

    # multipatchpre2D.plotMultipatchGeometry(geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F = linElastStat.assemblyMultipatchWeakForm(geomsurface,fullSurfacePreprocessing, \
          numericalquadrature,materialProperties,boundaryPreprocessing)

    Kred,Fred,removedDofs,totalDofs = matEqnSol.boundaryConditionsEnforcement(K,F,dirichletBCList)

    dtotal,D = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)

    # multipatchpost2D.postProcessing(geomsurface,D,dtotal, \
    #                                 fullSurfacePreprocessing,materialProperties)

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
