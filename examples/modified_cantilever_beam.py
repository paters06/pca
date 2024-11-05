# Python libraries
import numpy as np
# import matplotlib.pyplot as plt

#######################################################################
# DO NOT REMOVE THIS SEGMENT
from os.path import dirname
import sys

sys.path.append(dirname(dirname(__file__)))

#######################################################################

# Local project
import src.spline_functions.nurbs as rbs
import src.spline_functions.surfaceRefinements as srfn
from src.numerical_model import MultiPatchNumericalModel


def main_program():
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

    controlPointsPatch3 = np.array([[-0.2,-0.2],[0.0,-0.2],[-0.2,0.4],[0.0,0.4]])
    weightsPatch3 = np.ones((controlPointsPatch2.shape[0],1))

    # controlPointsPatch4 = np.array([[-0.2,0.1],[0.0,0.1],[-0.2,-0.2],[0.0,-0.2]])
    # weightsPatch4 = np.ones((controlPointsPatch2.shape[0],1))

    multiP = [controlPointsPatch1,controlPointsPatch2,controlPointsPatch3]
    multiw = [weightsPatch1,weightsPatch2,weightsPatch3]

    multip = [1,1,1]
    multiq = [1,1,1]

    # multiU = [np.array([0,0,0.5,0.5]),np.array([0.5,0.5,1,1])]
    multiU = [np.array([0,0,1,1]),np.array([0,0,1,1]),np.array([0,0,1,1])]
    multiV = [np.array([0,0,1,1]),np.array([0,0,1,1]),np.array([0,0,1,1])]

    geomsurface = rbs.MultiPatchNURBSSurface(multiU,multiV,multip,multiq,\
                                             multiP,multiw)

    localRefinement = True
    patchesToRefine = [0,1,2]
    numreflist = [2,2,2]
    # reflist = [['h'],['h'],['p'],['p'],['h'],['h']]
    # dirlist = [['U','V'],['U','V'],['U','V'],['U','V'],['U','V'],['U','V']]
    reflist = [['k'],['k'],['k']]
    dirlist = [['U','V'],['U','V'],['U','V']]

    if localRefinement:
        srfn.localPatchRefinement(geomsurface,patchesToRefine,numreflist,reflist,dirlist)

    #disp_i = [startpt,endpt,value,restriction]
    dirichletData_0 = [[0.0,0.0],[0.0,1.0],0.0,"C"]
    dirichletData_3 = [[0.0,0.0],[0.0,1.0],0.0,"C"]
    dirichletData_4 = [[0.0,0.0],[0.0,1.0],0.0,"C"]
    dirichletConditionsData = [None,None,dirichletData_3]
    #neumann_i = [startpt,endpt,type_load,value]
    neumannData_1 = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]
    neumannConditionsData = [None,neumannData_1,None]

    x_range = [0.0, 1.0]
    y_range = [1.0, 1.0]
    id_patches = [0, 1]

    cantilever_beam_mp = MultiPatchNumericalModel(phenomenon,geomsurface,
                                        dirichletConditionsData,neumannConditionsData,
                                        numGaussPoints,materialProperties,
                                        x_range, y_range, id_patches)
    cantilever_beam_mp.select_stage('Preprocessing')

if __name__ == '__main__':
    main_program()
    # profiling_script(mainProgram)
