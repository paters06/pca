tput reset
gfortran -c -g ../../src/spline_functions/utils.f90
gfortran -c -g ../../src/spline_functions/bspline_basis_functions.f90
gfortran -c -g ../../src/phenomenon/isogeometric_beam_1D.f90

# gfortran -c -g ../../src/spline_functions/nurbs_curve_module.f90 bspline_basis_functions.o
# gfortran -c -g -fbounds-check iga_1D.f90 nurbs_curve_module.o utils.o
gfortran -c -g -fbounds-check iga_1D.f90 utils.o isogeometric_beam_1D.o bspline_basis_functions.o
# gfortran -o nurbs_curve nurbs_curve_example.o nurbs_curve_module.o bspline_basis_functions.o utils.o
gfortran -L/usr/lib/x86_64-linux-gnu/ -o iga_1D iga_1D.o utils.o isogeometric_beam_1D.o bspline_basis_functions.o -llapack -lblas
rm *.o
rm *.mod
# cd ../../examples/geometry/
./iga_1D
