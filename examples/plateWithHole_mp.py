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

    Ra = 2.0
    Rb = 4.0
    # ka = cos(pi/2 - pi/8)
    ka = 0.5*np.sqrt(2 - np.sqrt(2))
    # kb = cos(pi/4)
    kb = 0.5*np.sqrt(2)
    # kw = cos(pi/8)
    kw = 0.5*np.sqrt(2 + np.sqrt(2))

    # Complete list of control points
    fullP = np.array([[0,Ra],[Ra*ka,Ra],[Ra*kb,Ra*kb],
                      [0,Rb],[Rb*ka,Rb],[Rb,Rb],
                      [Ra,Ra*ka],[Ra,0],[Rb,Rb*ka],
                      [Rb,0]])

    fullw = np.array([[1],[kw],[1],[1],[kw],[1],[kw],[1],[kw],[1]])

    idcontrolpoints = [[0,1,2,3,4,5],[2,6,7,5,8,9]]
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
    reflist = [['p','p','h','h'],['p','p','h','h']]
    dirlist = [['U','V','U','V'],['U','V','U','V']]
#    reflist = [['k','k','h','h'],['k','k','h','h']]
#    dirlist = [['U','V','U','V'],['U','V','U','V']]

    if localRefinement == 'Y':
        multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints = \
        srfn.localPatchRefinement(patchesToRefine,reflist,dirlist,multiU,multiV,multip,multiq,fullP,fullw,idcontrolpoints,localcontrolpoints)

#    for pi in fullP:
#        print(pi)

    #disp_i = [id patch,startpt,endpt,restriction,value]
    displacementConditions = [[0,[0.0,Ra],[0.0,Rb],"S",0.0],[1,[Ra,0.0],[Rb,0.0],"S",0.0]]
    #neumann_i = [id patch,startpt,endpt,type_load,value]
    neumannConditions = [[1,[0.0,1.0],[1.0,1.0],"normal",tv]]

    fullSurfacePreprocessing,boundaryPreprocessing,dirichletBCList = \
    multipatchpre2D.problemPreprocessing(multiU,multiV,multip,multiq,fullP,idcontrolpoints,displacementConditions,neumannConditions)

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
