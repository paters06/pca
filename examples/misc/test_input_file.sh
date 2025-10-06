tput reset
gfortran -c -g ../../src/spline_functions/utils.f90

gfortran -c -g -fbounds-check read_input_file.f90 utils.o
gfortran -L/usr/lib/x86_64-linux-gnu/ -o read_input_file read_input_file.o utils.o
rm *.o
rm *.mod
# cd ../../examples/geometry/
./read_input_file
