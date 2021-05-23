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
    phenomenon = "Elasticity"
    E = 1e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = -10 #Pa
    L = 4.0
    R = 1.0
    numGaussPoints = 4

    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[-R,0],[-R,R*(np.sqrt(2)-1)],[R*(1-np.sqrt(2)),R],[0,R],[-2.5,0],[-2.5,0.75],
                  [-0.75,2.5],[0,2.5],[-L,0],[-L,L],[-L,L],[0,L]])

    winit = np.array([[1],[0.5*(1+(1/np.sqrt(2)))],[0.5*(1+(1/np.sqrt(2)))],[1],[1],[1],
                  [1],[1],[1],[1],[1],[1]])

    #Isogeometric routines
    Uinit = np.array([0,0,0,0.5,1,1,1])
    Vinit = np.array([0,0,0,1,1,1])

    pinit = 2
    qinit = 2

    geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)

    doRefinement = 'Y'

    if doRefinement == 'Y':
        # srfn.surfaceRefinement(geomsurface,1,'p','U')
        # srfn.surfaceRefinement(geomsurface,1,'p','V')
        srfn.surfaceRefinement(geomsurface,3,'h','U')
        srfn.surfaceRefinement(geomsurface,3,'h','V')

    # displacementConditions = [[[-R,0.0],[-L,0.0],"S",0.0],[[0.0,R],[0.0,L],"S",0.0]]
    # neumannConditions = [[[0.0,1.0],[0.5,1.0],"normal",tv]]

    dirichletConditionsData = [[[0.0,0.0],[0.0,1.0],0.0,"S"],[[1.0,0.0],[1.0,1.0],0.0,"S"]]
    neumannConditionsData = [[[0.0,1.0],[0.5,1.0],"normal",tv]]

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F = linElastStat.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing,neumannConditionsData)

    Kred,Fred,totalDofs,removedDofs,dofValues = matEqnSol.dirichletBCEnforcement(phenomenon,K,F,dirichletBCList)

    dtotal,D = matEqnSol.solveMatrixEquations(phenomenon,Kred,Fred,totalDofs,removedDofs,dofValues)

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
