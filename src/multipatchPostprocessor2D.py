# Python libraries
import math
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs
from src.linearElastoStaticsSolver import elasticMatrix
import src.plottingScripts as plts

################ POSTPROCESSING ####################

def find_match_rows(A_mat: np.ndarray, B_mat: np.ndarray) -> np.ndarray:
    """
    A is the sample matrix
    B is the global matrix

    This functions return the indices of the rows of the B matrix
    that corresponds to each row of the A matrix
    """
    match_indices_list = []

    for i in range(0, A_mat.shape[0]):
        for j in range(0, B_mat.shape[0]):
            equal_row = np.where(A_mat[i,:] == B_mat[j,:])[0]
            if len(equal_row) == 2:
                match_indices_list.append(j)
    
    match_indices = np.array(match_indices_list, dtype=int)
    return match_indices

def defineNumberOfEvaluationPoints(numelems):
    if numelems < 5:
        numpoints = 11
    elif numelems >= 5 and numelems < 10:
        numpoints = 9
    elif numelems >= 10 and numelems < 20:
        numpoints = 7
    else:
        numpoints = 5

    return numpoints

class SolutionField:
    """
    A class for the multipatch solution field
    """
    def __init__(self,multisurface,dtot, surfaceprep):
        """
        Constructor of the class
        """
        multiU,multiV,multip,multiq,multiP,multiw,globalPatchIndices = multisurface.retrieveSurfaceInformation()
        self.multiU = multiU
        self.multiV = multiV
        self.multip = multip
        self.multiq = multiq
        self.multiP = multiP
        self.multiw = multiw
        self.globalPatchIndices = globalPatchIndices
        self.fullP = multisurface.fullP
        # self.Dsol = D
        self.dtot = dtot
        self.Dsol = self._create_D_matrix_from_d_vector(dtot)
        self._calculate_field_matrix_size(surfaceprep)

    def _calculate_field_matrix_size(self, surfaceprep):
        numpatches = len(self.multiU)

        totalpoints = 0

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.multiP[ipatch]
            wi = self.multiw[ipatch]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu = mu - pi - 1
            nv = mv - qi - 1

            elementcorners = surfaceprep[ipatch][2]

            numElems = len(elementcorners)
            numelemsu = len(np.unique(Ui)) - 1
            numelemsv = len(np.unique(Vi)) - 1

            numpoints = defineNumberOfEvaluationPoints(numElems)

            totalpoints += (numpoints*numelemsu*numpoints*numelemsv)

        # self.totalpoints = totalpoints

        self.fullcpts = np.zeros((totalpoints, 2))
        self.fullupts = np.zeros((totalpoints, 2))
        self.fullsigmapts = np.zeros((totalpoints, 3))

    def _create_D_matrix_from_d_vector(self, d_vec):
        num_dofs = d_vec.shape[0]
        num_ctrl_pts = int(num_dofs/2)
        Dsol = np.zeros((num_ctrl_pts, 2))
        Dsol[:,None,0] = d_vec[0::2]
        Dsol[:,None,1] = d_vec[1::2]
        return Dsol

    def displacement_field(self, surfaceprep):
        """
        Computing the displacement field for each patch
        """
        numpatches = len(self.multiU)

        id_point = 0

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.multiP[ipatch]
            wi = self.multiw[ipatch]

            ids_patch = find_match_rows(Pi, self.fullP)

            Di = self.Dsol[self.globalPatchIndices[ipatch],:]
            Pi2 = self.fullP[self.globalPatchIndices[ipatch],:]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu = mu - pi - 1
            nv = mv - qi - 1

            # Extraction of surface preprocessing
            nonzeroctrlpts = surfaceprep[ipatch][0]
            surfacespan = surfaceprep[ipatch][1]
            elementcorners = surfaceprep[ipatch][2]

            numElems = len(elementcorners)
            numelemsu = len(np.unique(Ui)) - 1
            numelemsv = len(np.unique(Vi)) - 1

            numpoints = defineNumberOfEvaluationPoints(numElems)

            Pwl = rbs.weightedControlPoints(Pi,wi)
            Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            # Geometric coordinates
            cpts = np.zeros((1,Pi.shape[1]))

            # Displacements
            upts = np.zeros((1,Pi.shape[1]))

            for ielem in range(0,numElems):
                # Extracting the indices of the non-zero control points
                idR = nonzeroctrlpts[ielem]
                # Extracting the indices for the location of the parametric element
                uspan = surfacespan[ielem][0]
                vspan = surfacespan[ielem][1]
                # Extracting the corners of the parametric element
                apt = elementcorners[ielem][0]
                cpt = elementcorners[ielem][1]

                urank = np.linspace(apt[0],cpt[0],numpoints)
                vrank = np.linspace(apt[1],cpt[1],numpoints)

                jv = ielem//numelemsu
                iu = ielem%numelemsu

                for j in range(0,numpoints):
                    for i in range(0,numpoints):
                        R = rbs.bivariateRationalFunction(mu,mv,pi,qi,uspan,vspan,urank[i],vrank[j],Ui,Vi,Pwi)

                        upts = R@Di[idR,:]
                        cpts = R@Pi[idR,:]

                        self.fullcpts[id_point,:] = cpts[0]
                        self.fullupts[id_point,:] = upts[0]
                        id_point += 1

    def stress_field(self,matprop,surfaceprep):
        """
        Computing the stress field in list form
        """
        numpatches = len(self.multiU)

        id_point = 0

        # Definition of the material matrix
        E = matprop[0]
        nu = matprop[1]
        rho = matprop[2]
        dMat = elasticMatrix(E,nu)

        paramgrad = np.zeros((2,2))

        for ipatch in range(0,numpatches):
            Ui = self.multiU[ipatch]
            Vi = self.multiV[ipatch]

            pi = self.multip[ipatch]
            qi = self.multiq[ipatch]

            Pi = self.multiP[ipatch]
            wi = self.multiw[ipatch]
            Di = self.Dsol[self.globalPatchIndices[ipatch],:]

            mu = len(Ui) - 1
            mv = len(Vi) - 1
            nu = mu - pi - 1
            nv = mv - qi - 1

            # Extraction of surface preprocessing
            nonzeroctrlpts = surfaceprep[ipatch][0]
            surfacespan = surfaceprep[ipatch][1]
            elementcorners = surfaceprep[ipatch][2]

            numElems = len(elementcorners)
            numelemsu = len(np.unique(Ui)) - 1
            numelemsv = len(np.unique(Vi)) - 1

            numpoints = defineNumberOfEvaluationPoints(numElems)

            Pwl = rbs.weightedControlPoints(Pi,wi)
            Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)

            sigma = np.zeros((1,3))

            # Global degrees of freedom
            globalDOF = np.zeros(2*len(self.globalPatchIndices[ipatch]),dtype=int)
            dof0 = 2*np.array(self.globalPatchIndices[ipatch])
            dof1 = dof0 + 1
            globalDOF[0::2] = dof0
            globalDOF[1::2] = dof1

            for ielem in range(0,numElems):
                # Extracting the indices of the non-zero control points
                idR = nonzeroctrlpts[ielem]
                # Extracting the indices for the location of the parametric element
                uspan = surfacespan[ielem][0]
                vspan = surfacespan[ielem][1]
                # Extracting the corners of the parametric element
                apt = elementcorners[ielem][0]
                cpt = elementcorners[ielem][1]

                urank = np.linspace(apt[0],cpt[0],numpoints)
                vrank = np.linspace(apt[1],cpt[1],numpoints)

                # Patch degrees of freedom
                patchDOF = np.zeros(2*len(idR),dtype=int)
                dof0 = 2*np.array(idR)
                dof1 = dof0 + 1
                patchDOF[0::2] = dof0
                patchDOF[1::2] = dof1

                for j in range(0,numpoints):
                    for i in range(0,numpoints):
                        xpcoor = urank[i]
                        ypcoor = vrank[j]

                        biRatGrad = rbs.bivariateRationalGradient(mu,mv,pi,qi,uspan,vspan,xpcoor,ypcoor,Ui,Vi,Pwi)

                        jac = (biRatGrad[1:3,:]@Pi[idR,:]).T
                        detJac = np.linalg.det(jac)

                        if abs(detJac) > 1e-5:
                            invJac = np.linalg.inv(jac)
                            dN2 = biRatGrad[1:3,:]
                            dN2dxi = invJac.T@dN2

                            numpts = dN2dxi.shape[1]
                            bmat = np.zeros((3,2*numpts))
                            #dNx
                            bmat[0,0::2] = dN2dxi[0,:]
                            bmat[2,0::2] = dN2dxi[1,:]
                            #dNy
                            bmat[1,1::2] = dN2dxi[1,:]
                            bmat[2,1::2] = dN2dxi[0,:]

                            sigma = dMat@(bmat@self.dtot[globalDOF[patchDOF],:])
                        else:
                            print("Singularity")
                        self.fullsigmapts[id_point,:] = sigma.flatten()
                        id_point += 1
    
    def showExtremaValues(self):
        print("Displacements")
        print("UX ==> Max: {:.5f} m. Min: {:.5f} m".format(np.max(self.fullupts[:,0]),np.min(self.fullupts[:,0])))
        print("UY ==> Max: {:.5f} m. Min: {:.5f} m".format(np.max(self.fullupts[:,1]),np.min(self.fullupts[:,1])))
        print("Stresses")
        print("SXX ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.fullsigmapts[:,0]),np.min(self.fullsigmapts[:,0])))
        print("SYY ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.fullsigmapts[:,1]),np.min(self.fullsigmapts[:,1])))
        print("SXY ==> Max: {:.3f} Pa. Min: {:.3f} Pa".format(np.max(self.fullsigmapts[:,2]),np.min(self.fullsigmapts[:,2])))

################ MAIN POSTPROCESSING FUNCTION ####################

def postProcessing(phenomenon,multisurface,surfaceprep,dtot,matprop):
    solfield = SolutionField(multisurface,dtot,surfaceprep)

    solfield.displacement_field(surfaceprep)
    solfield.stress_field(matprop,surfaceprep)

    solfield.showExtremaValues()

    # plts.plot_multipatch_grid(solfield.fullcpts)
    # plts.plot_multipatch_field(solfield.fullcpts,solfield.fullupts,0,["Ux Displacement Field","[m]"])
    plts.plot_multipatch_field(solfield.fullcpts,solfield.fullsigmapts,0,["Sx Stress Field","[Pa]"])
