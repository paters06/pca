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
import src.surfaceRefinements as srfn
# import src.debugScripts as dbg_scrpt
from src.numerical_model import NumericalModel


def main_program():
    # Data
    phenomenon = "Elasticity"
    E = 2e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = -1 #Pa
    numGaussPoints = 4

    Ra = 2.0
    Rb = 4.0
    # ka = cos(pi/2 - pi/8)
    ka = 0.5*np.sqrt(2 - np.sqrt(2))
    # kb = cos(pi/4)
    kb = 0.5*np.sqrt(2)
    # kw = cos(pi/8)
    kw = 0.5*np.sqrt(2 + np.sqrt(2))

    controlPointsPatch1 = np.array([[0,Ra],[Ra*ka,Ra],[Ra*kb,Ra*kb],
                                    [0,Rb],[Rb*ka,Rb],[Rb,Rb]])
    weightsPatch1 = np.array([[1],[kw],[1],[1],[kw],[1]])

    controlPointsPatch2 = np.array([[Ra*kb,Ra*kb],[Ra,Ra*ka],[Ra,0],[Rb,Rb],[Rb,Rb*ka],
                                    [Rb,0]])
    weightsPatch2 = np.array([[1],[kw],[1],[1],[kw],[1]])

    multiP = [controlPointsPatch1,controlPointsPatch2]
    multiw = [weightsPatch1,weightsPatch2]

    multip = [2,2]
    multiq = [1,1]

    multiU = [np.array([0,0,0,1,1,1]),np.array([0,0,0,1,1,1])]
    multiV = [np.array([0,0,1,1]),np.array([0,0,1,1])]

    geomsurface = rbs.MultiPatchNURBSSurface(multiU,multiV,multip,multiq,\
                                             multiP,multiw)

    localRefinement = False
    patchesToRefine = [0,1]
    numreflist = [2,2]
    reflist = [['h'],['h']]
    dirlist = [['U','V'],['U','V']]

    if localRefinement:
        srfn.localPatchRefinement(geomsurface,patchesToRefine,numreflist,reflist,dirlist)

    # disp_i = [startpt,endpt,value,restriction]
    dirichletData_0 = [[0.0,0.0],[0.0,1.0],0.0,"S"]
    dirichletData_1 = [[1.0,0.0],[1.0,1.0],0.0,"S"]
    dirichletConditionsData = [dirichletData_0,dirichletData_1]
    #neumann_i = [startpt,endpt,type_load,value]
    neumannData_1 = [[[0.0,1.0],[1.0,1.0],"normal",tv]]
    neumannConditionsData = [None,neumannData_1]

    plate_with_hole_mp = NumericalModel(phenomenon,geomsurface,
                                        dirichletConditionsData,neumannConditionsData,
                                        numGaussPoints,materialProperties)
    plate_with_hole_mp.select_stage('Preprocessing')


if __name__ == '__main__':
    main_program()
    # profiling_script(mainProgram)
