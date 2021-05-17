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
    Rmax = 1.0
    Rmin = 0.7
    kappa = 385 #Pa
    source = 0.0 #kg/m3
    materialProperties = [kappa,source]
    flux = 0.0 #Pa

    numGaussPoints = 4
    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[Rmin,0],[Rmin,Rmin],[0,Rmin],[Rmax,0],[Rmax,Rmax],[0,Rmax]])

    winit = np.array([[1],[0.5*np.sqrt(2)],[1],[1],[0.5*np.sqrt(2)],[1]])

    #Isogeometric routines
    Uinit = np.array([0,0,0,1,1,1])
    Vinit = np.array([0,0,1,1])

    pinit = 2
    qinit = 1

    geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)

    doRefinement = 'N'

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

    displacementConditionsData = [[[0.0,Rmin],[0.0,Rmax],"S",0.0],[[Rmin,0.0],[Rmax,0.0],"S",0.0]]
    neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"normal",tv]]

    # neumannConditionsData = [[[0.0,1.0],[1.0,1.0],"tangent",flux]]
    # neumannConditionsData = [[[0.0,0.0],[1.0,0.0],"tangent",flux],[[0.0,1.0],[1.0,1.0],"tangent",flux]]

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,displacementConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F = diffSol.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)

    Kred,Fred,removedDofs,totalDofs = matEqnSol.dirichletBCEnforcement(phenomenon,K,F,dirichletBCList)

    dtotal,D = matEqnSol.solveMatrixEquations(phenomenon,Kred,Fred,totalDofs,removedDofs,dirichletBCList)

    post2D.postProcessing(phenomenon,geomsurface,D,dtotal,surfacePreprocessing,materialProperties)

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
