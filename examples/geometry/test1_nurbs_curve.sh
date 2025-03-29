tput reset
gfortran -c -g ../../src/spline_functions/utils.f90
gfortran -c -g ../../src/spline_functions/bspline_basis_functions.f90
gfortran -c -g ../../src/spline_functions/nurbs_curve_module.f90 bspline_basis_functions.o
gfortran -c -g -fbounds-check ../../src/spline_functions/nurbs_curve_example.f90 nurbs_curve_module.o utils.o
gfortran -o nurbs_curve nurbs_curve_example.o nurbs_curve_module.o bspline_basis_functions.o utils.o
rm *.o
rm *.mod
./nurbs_curve
