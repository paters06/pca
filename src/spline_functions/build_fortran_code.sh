gfortran -c -o one_basis_function.o one_basis_function.f90 
gfortran -c -o test_fortran_scripts.o test_fortran_scripts.f90
gfortran -c -o find_span.o find_span.f90
gfortran -c -o basis_functions.o basis_functions.f90
gfortran -o test one_basis_function.o find_span.o test_fortran_scripts.o basis_functions.o