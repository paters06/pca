# Python libraries
import numpy as np
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
from src.profiling_script import profiling_script
import src.nurbs as rbs
import src.surfaceRefinements as srfn
from src.numerical_model import NumericalModel


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

    multiP = [controlPointsPatch1,controlPointsPatch2]
    multiw = [weightsPatch1,weightsPatch2]

    multip = [1,1]
    multiq = [1,1]

    # multiU = [np.array([0,0,0.5,0.5]),np.array([0.5,0.5,1,1])]
    multiU = [np.array([0,0,1,1]),np.array([0,0,1,1])]
    multiV = [np.array([0,0,1,1]),np.array([0,0,1,1])]

    geomsurface = rbs.MultiPatchNURBSSurface(multiU,multiV,multip,multiq,\
                                             multiP,multiw)

    localRefinement = True
    patchesToRefine = [0,1]
    numreflist = [2,2]
    # reflist = [['h'],['h'],['p'],['p'],['h'],['h']]
    reflist = [['h'],['h']]
    # dirlist = [['U','V'],['U','V'],['U','V'],['U','V'],['U','V'],['U','V']]
    dirlist = [['U','V'],['U','V']]

    if localRefinement:
        srfn.localPatchRefinement(geomsurface,patchesToRefine,numreflist,reflist,dirlist)

    #disp_i = [startpt,endpt,value,restriction]
    dirichletData_0 = [[0.0,0.0],[0.0,1.0],0.0,"C"]
    dirichletConditionsData = [dirichletData_0,None]
    #neumann_i = [startpt,endpt,type_load,value]
    neumannData_1 = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]
    neumannConditionsData = [None,neumannData_1]

    plate_with_hole_mp = NumericalModel(phenomenon,geomsurface,
                                        dirichletConditionsData,neumannConditionsData,
                                        numGaussPoints,materialProperties)
    plate_with_hole_mp.select_stage('Preprocessing')

if __name__ == '__main__':
    main_program()
    # profiling_script(mainProgram)
