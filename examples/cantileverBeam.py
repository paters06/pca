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
from src.profiling_script import profiling_script
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

    Pinit = np.array([[0,0],[0.5*L,0],[L,0],[0,0.5*H],
                  	 [0.5*L,0.5*H],[L,0.5*H],[0,H],[0.5*L,H],[L,H]])

    winit = np.array([[1],[1],[1],[1],[1],[1],[1],[1],[1]])

    # Pinit=np.array([[0,0],[3.75,0],[7.5,0],[11.25,0],[15,0],[18.75,0],[22.5,0],
    #                 [26.25,0],[30,0],[0,3],[3.75,3],[7.5,3],[11.25,3],[15,3],
    #                 [18.75,3],[22.5,3],[26.25,3],[30,3],[0,6],[3.75,6],[7.5,6],
    #                 [11.25,6],[15,6],[18.75,6],[22.5,6],[26.25,6],[30,6]])
    #
    # winit = np.ones((Pinit.shape[0],1))

    gridsize = [9,3]

    #Isogeometric routines
    Uinit = np.array([0,0,0,1,1,1])
    Vinit = np.array([0,0,0,1,1,1])
    # Uinit = np.array([0,0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1,1])
    # Vinit = np.array([0,0,0.5,1,1])

    pinit = 2
    qinit = 2

    geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,U=Uinit,V=Vinit)
    # geomsurface = rbs.NURBSSurface(Pinit,winit,pinit,qinit,gridsize=[9,3])

    doRefinement = 'N'

    if doRefinement == 'Y':
        # srfn.surfaceRefinement(geomsurface,1,'p','U')
        # srfn.surfaceRefinement(geomsurface,1,'p','V')
        srfn.surfaceRefinement(geomsurface,1,'h','U')
        srfn.surfaceRefinement(geomsurface,1,'h','V')

    dirichletConditionsData = [[[0.0,0.0],[0.0,1.0],0.0,"C"]]
    neumannConditionsData = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]

    surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
    pre2D.problemPreprocessing(phenomenon,geomsurface,dirichletConditionsData,neumannConditionsData)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)

    # dbg_scrpt.calculateArea(geomsurface,surfacePreprocessing,numericalquadrature)

    # pre2D.plotGeometry(phenomenon,geomsurface,dirichletBCList,boundaryPreprocessing)

    K,F,M = linElastStat.assemblyWeakForm(geomsurface,surfacePreprocessing,numericalquadrature,\
                                        materialProperties,boundaryPreprocessing)

    Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,enforcedDOF,enforcedValues)

    dSolution = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,enforcedDOF,enforcedValues)

    # post2D.postProcessing(phenomenon,geomsurface,surfacePreprocessing,dSolution,materialProperties)

if __name__ == '__main__':
    # mainProgram()
    profiling_script(mainProgram)
