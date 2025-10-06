tput reset
gfortran -c -fcheck=all -Wall -g ../../src/spline_functions/utils.f90
gfortran -c -fcheck=all -Wall -g -fbounds-check ../../src/spline_functions/curve_refinement.f90
gfortran -c -fcheck=all -Wall -g ../../src/spline_functions/bspline_basis_functions.f90
gfortran -c -fcheck=all -Wall -g ../../src/spline_functions/nurbs_curve_module.f90
# gfortran -c -fcheck=all -Wall -g ../../src/phenomenon/reference_solutions.f90
# gfortran -c -fcheck=all -Wall -g -fbounds-check ../../src/phenomenon/euler_bernoulli_beam.f90

gfortran -c -fcheck=all -Wall -g -fbounds-check refinement_example_3.f90 curve_refinement.o utils.o bspline_basis_functions.o nurbs_curve_module.o
gfortran -fcheck=all -Wall -o refinement_example_3 refinement_example_3.o curve_refinement.o utils.o bspline_basis_functions.o nurbs_curve_module.o

# gfortran -c -g ../../src/spline_functions/nurbs_curve_module.f90 bspline_basis_functions.o
# gfortran -c -g -fbounds-check beam_example_1.f90 nurbs_curve_module.o utils.o
# gfortran -c -fcheck=all -Wall -g -fbounds-check beam_example_1.f90 utils.o euler_bernoulli_beam.o nurbs_curve_module.o bspline_basis_functions.o reference_solutions.o
# gfortran -o nurbs_curve nurbs_curve_example.o nurbs_curve_module.o bspline_basis_functions.o utils.o
# gfortran -L/usr/lib/x86_64-linux-gnu/ -fcheck=all -Wall -o beam_example_1 beam_example_1.o utils.o euler_bernoulli_beam.o nurbs_curve_module.o bspline_basis_functions.o reference_solutions.o -llapack -lblas
rm *.o
rm *.mod
# cd ../../examples/geometry/
./refinement_example_3
