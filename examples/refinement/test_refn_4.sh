tput reset
gfortran -c -fcheck=all -Wall -Wno-uninitialized -g -fbounds-check ../../src/spline_functions/utils.f90
gfortran -c -fcheck=all -Wall -g -fbounds-check ../../src/spline_functions/curve_refinement.f90
gfortran -c -fcheck=all -Wall -g ../../src/spline_functions/bspline_basis_functions.f90
gfortran -c -fcheck=all -Wall -g ../../src/spline_functions/nurbs_curve_module.f90

gfortran -c -fcheck=all -Wall -Wno-uninitialized -g -fbounds-check refinement_example_4.f90 curve_refinement.o utils.o bspline_basis_functions.o nurbs_curve_module.o
gfortran -fcheck=all -Wall -o refinement_example_4 refinement_example_4.o curve_refinement.o utils.o bspline_basis_functions.o nurbs_curve_module.o

rm *.o
rm *.mod
./refinement_example_4
