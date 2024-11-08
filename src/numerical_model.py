import src.preprocessor.preprocessor2D as pre2D
import src.postprocessor.postprocessor2D as post2D
import src.multipatch.multipatchPreprocessor2D as multipatchpre2D
import src.phenomenon.linearElastoStaticsSolver as linElastStat
import src.matrixEquationSolver as matEqnSol
import src.multipatch.multipatchPostprocessor2D as multipatchpost2D
# import src.nurbs as rbs
from src.spline_functions.nurbs import MultiPatchNURBSSurface
from src.spline_functions.nurbs import NURBSSurface


# import numpy as np

class NumericalModel:
    def __init__(self, phenomenon: str, geomsurface: NURBSSurface,
                 dirichletConditionsData, neumannConditionsData,
                 numGaussPoints: int, materialProperties: list) -> None:
        self.phenomenon = phenomenon
        self.geomsurface = geomsurface
        self.dirichlet_conditions = dirichletConditionsData
        self.neumann_conditions = neumannConditionsData
        self.gauss_points_number = numGaussPoints
        self.material_properties = materialProperties

    def preprocessing(self):
        surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
        pre2D.problemPreprocessing(self.phenomenon,self.geomsurface,
                                   self.dirichlet_conditions,
                                   self.neumann_conditions)
        
        self.surface_preprocessing = surfacePreprocessing
        self.boundary_preprocessing = boundaryPreprocessing
        self.enforced_DOF = enforcedDOF
        self.enforced_values = enforcedValues
        
        self.numericalquadrature = pre2D.numericalIntegrationPreprocessing(self.gauss_points_number)

        pre2D.plotGeometry(self.phenomenon,self.geomsurface,dirichletBCList,boundaryPreprocessing)

    def analysis(self):
        K,F,M = linElastStat.assemblyWeakForm(self.geomsurface,self.surface_preprocessing,
                                              self.numericalquadrature,
                                              self.material_properties,
                                              self.boundary_preprocessing)

        Mred,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(M,K,F,self.enforced_DOF,self.enforced_values)

        self.d_solution = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,self.enforced_DOF,self.enforced_values)

    def postprocessing(self):
        post2D.postProcessing(self.phenomenon,self.geomsurface,
                              self.surface_preprocessing,
                              self.d_solution,self.material_properties)

    def select_stage(self, stage: str):
        if stage == 'Preprocessing':
            self.preprocessing()
        elif stage == 'Analysis':
            self.preprocessing()
            self.analysis()
        elif stage == 'Postprocessing':
            self.preprocessing()
            self.analysis()
            self.postprocessing()
        else:
            print('Error')


class MultiPatchNumericalModel:
    def __init__(self, phenomenon: str, geomsurface: MultiPatchNURBSSurface,
                 dirichletConditionsData, neumannConditionsData,
                 numGaussPoints: int, materialProperties: list,
                 x_pt, y_pt, id_patches):
        self.phenomenon = phenomenon
        self.geomsurface = geomsurface
        self.dirichlet_conditions = dirichletConditionsData
        self.neumann_conditions = neumannConditionsData
        self.gauss_points_number = numGaussPoints
        self.material_properties = materialProperties
        self.path_start = x_pt
        self.path_end = y_pt
        self.patches_on_path = id_patches

    def preprocessing(self):
        surfacePreprocessing,boundaryPreprocessing,dirichletBCList,enforcedDOF,enforcedValues = \
        multipatchpre2D.multiPatchProblemPreprocessing(self.phenomenon,self.geomsurface,
                                                       self.dirichlet_conditions,
                                                       self.neumann_conditions)
        
        self.surface_preprocessing = surfacePreprocessing
        self.boundary_preprocessing = boundaryPreprocessing
        self.enforced_DOF = enforcedDOF
        self.enforced_values = enforcedValues

        self.numericalquadrature = pre2D.numericalIntegrationPreprocessing(self.gauss_points_number)

        multipatchpre2D.plotMultiPatchGeometry(self.phenomenon,self.geomsurface,
                                               dirichletBCList,boundaryPreprocessing,
                                               self.path_start, self.path_end, self.patches_on_path)

    def analysis(self):
        Ktotal,Ftotal,Mtotal = linElastStat.assemblyMultipatchWeakForm(self.geomsurface,
                                                                       self.surface_preprocessing,
                                                                       self.numericalquadrature,
                                                                       self.material_properties,
                                                                       self.boundary_preprocessing)

        Mtotal,Kred,Fred,totalDofs = matEqnSol.dirichletBCEnforcement(Mtotal,Ktotal,Ftotal,
                                                                      self.enforced_DOF,self.enforced_values)

        self.d_solution = matEqnSol.solveMatrixEquations(Kred,Fred,totalDofs,self.enforced_DOF,self.enforced_values)

    def postprocessing(self):
        multipatchpost2D.postProcessing(self.phenomenon,self.geomsurface,
                                        self.surface_preprocessing,self.d_solution,
                                        self.material_properties)
    
    def path_postprocessing(self):
        
        multipatchpost2D.pathPostProcessing(self.phenomenon,self.geomsurface,
                                            self.surface_preprocessing,self.d_solution,
                                            self.material_properties)

    def select_stage(self, stage: str):
        if stage == 'Preprocessing':
            self.preprocessing()
        elif stage == 'Analysis':
            self.preprocessing()
            self.analysis()
        elif stage == 'Postprocessing':
            self.preprocessing()
            self.analysis()
            self.postprocessing()
        elif stage == 'Path_postprocessing':
            self.preprocessing()
            self.analysis()
            self.path_postprocessing()
        else:
            print('Error')
