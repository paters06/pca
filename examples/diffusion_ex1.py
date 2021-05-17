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
import src.diffusionSolver as diffSol
import src.matrixEquationSolver as matEqnSol
import src.postprocessor2D as post2D
import src.surfaceRefinements as srfn
import src.debugScripts as dbg_scrpt

####################################################
################## MAIN PROBLEM ####################
####################################################

def mainProgram():
    #Data
    phenomenon = "Heat"
    kappa = 385 #Pa
    source = 0.0 #kg/m3
    materialProperties = [kappa,source]
    flux = 0.0 #Pa

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0.0,0.0],[1.0,0.0],[0.0,1.0],[1.0,1.0]])

    # Pinit=np.array([[0,0],[3.75,0],[7.5,0],[11.25,0],[15,0],[18.75,0],[22.5,0],
    #                 [26.25,0],[30,0],[0,3],[3.75,3],[7.5,3],[11.25,3],[15,3],
    #                 [18.75,3],[22.5,3],[26.25,3],[30,3],[0,6],[3.75,6],[7.5,6],
    #                 [11.25,6],[15,6],[18.75,6],[22.5,6],[26.25,6],[30,6]])

    winit = np.ones((Pinit.shape[0],1))

    #Isogeometric routines
    Uinit = np.array([0,0,1,1])
    Vinit = np.array([0,0,1,1])
    # Uinit = np.array([0,0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1,1])
    # Vinit = np.array([0,0,0.5,1,1])

    pinit = 1
    qinit = 1

    geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)

    doRefinement = 'Y'

    if doRefinement == 'Y':
        # reflist = ['h','h','h','h']
        # dirlist = ['U','V','U','V']
        reflist = ['p','p','h','h']
        dirlist = ['U','V','U','V']
        # reflist = ['h','h','h','h','h','h','h','h','h','h']
        # dirlist = ['U','V','U','V','U','V','U','V','U','V']
        # reflist = ['p','p','h','h','h','h','h','h','h','h']
        # dirlist = ['U','V','U','V','U','V','U','V','U','V']
        srfn.surfaceRefinement(geomsurface,reflist,dirlist)

    displacementConditionsData = [[[0.0,0.0],[0.0,1.0],"C",0.0],[[1.0,0.0],[1.0,1.0],"C",100.0]]
    neumannConditionsData = [[[0.0,1.0],[1.0,1.0],"tangent",flux]]
    # neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"tangent",5],[[0.0,1.0],[1.0,1.0],"tangent",10]]

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList = \
    pre2D.problemPreprocessing(geomsurface,displacementConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F = diffSol.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)

    Kred,Fred,removedDofs,totalDofs = matEqnSol.boundaryConditionsEnforcement(K,F,dirichletBCList)

    # dtotal,D = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)

    # post2D.postProcessing(geomsurface,D,dtotal,surfacePreprocessing,materialProperties)

mainProgram()

# import cProfile
# import pstats
# profiler = cProfile.Profile()
# profiler.enable()
# mainProgram()
# profiler.disable()
# # stats = pstats.Stats(profiler).sort_stats('ncalls')
# stats = pstats.Stats(profiler).sort_stats('tottime')
# stats.print_stats()
