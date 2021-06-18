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
    phenomenon = "Elasticity"
    E = 2e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = -1 #Pa
    numGaussPoints = 4

    controlPointsPatch1 = np.array([[0,0],[0.4,0],[0.0,0.2],[0.4,0.2]])
    weightsPatch1 = np.ones((controlPointsPatch1.shape[0],1))

    controlPointsPatch2 = np.array([[0.4,0.0],[1.0,0.0],[0.4,0.2],[1.0,0.2]])
    weightsPatch2 = np.ones((controlPointsPatch2.shape[0],1))

    multiP = [controlPointsPatch1,controlPointsPatch2]
    multiw = [weightsPatch1,weightsPatch2]

    multip = [1,1]
    multiq = [1,1]

    # multiU = [np.array([0,0,0.5,0.5]),np.array([0.5,0.5,1,1])]
    multiU = [np.array([0,0,1,1]),np.array([0,0,1,1])]
    multiV = [np.array([0,0,1,1]),np.array([0,0,1,1])]

    geomsurface = rbs.MultiPatchNURBSSurface(multiU,multiV,multip,multiq,\
                                             multiP,multiw)

    localRefinement = 'N'
    patchesToRefine = [0,1]
    numreflist = [2,1]
    reflist = [['h'],['h']]
    dirlist = [['U','V'],['U','V']]

    if localRefinement == 'Y':
        srfn.localPatchRefinement(geomsurface,patchesToRefine,numreflist,reflist,dirlist)
    # End if

    #disp_i = [startpt,endpt,value,restriction]
    dirichletData_0 = [[0.0,0.0],[0.0,1.0],0.0,"C"]
    dirichletConditionsData = [dirichletData_0,None]
    #neumann_i = [startpt,endpt,type_load,value]
    neumannData_1 = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]
    neumannConditionsData = [None,neumannData_1]

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
    multipatchpre2D.multiPatchProblemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # multipatchpre2D.plotMultiPatchGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    Ktotal,Ftotal,Mtotal = linElastStat.assemblyMultipatchWeakForm(geomsurface,surfacePreprocessing,\
            numericalquadrature,materialProperties,boundaryPreprocessing)

    Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement_Reduced(Ktotal,Ftotal,enforcedDOF,enforcedValues)

    dtotal,D = matEqnSol.solveReducedMatrixEquations(phenomenon,Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    multipatchpost2D.postProcessing(geomsurface,D,dtotal,surfacePreprocessing,materialProperties)

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
