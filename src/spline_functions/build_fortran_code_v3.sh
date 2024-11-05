gfortran -c -o bspline_basis_functions.o bspline_basis_functions.f90
gfortran -c -o test_fortran_scripts_v2.o test_fortran_scripts_v2.f90
gfortran -o test test_fortran_scripts_v2.o bspline_basis_functions.o