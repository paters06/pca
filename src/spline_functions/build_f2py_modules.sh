# python3 -m numpy.f2py -c -m bspline_basis_functions_pmod bspline_basis_functions.f90
python3 -m numpy.f2py -c -m nurbs_curve_pymod nurbs_curve_module.f90 bspline_basis_functions.f90