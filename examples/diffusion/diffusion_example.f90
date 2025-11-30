program diffusion_example
    use utils
    use nurbs_surface_module
    use surface_refinement
    use input_output
    use diffusion_solver
    use derived_types
    implicit none

    character(:), allocatable :: file_name
    character(len=50) :: file_name_inp
    character(len=50), dimension(:), allocatable :: line_array
    character(len=1), dimension(:,:), allocatable :: ref_list
    integer :: num_gauss_pts
    real :: kappa
    real, dimension(:,:), allocatable :: Kmat, Fvec, Kred, Usol
    real, dimension(:), allocatable :: Fred
    integer, dimension(:), allocatable :: id_disp, remainder_dofs
    real, dimension(:), allocatable :: u_pres
    type(boundary_condition), dimension(:), allocatable :: bc_array
    type(nurbs_surface) :: input_nurbs_surface, refined_surface

    call get_command_argument(1, file_name_inp)
    file_name = file_name_inp

    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, input_nurbs_surface, ref_list)
    call convert_data_to_solver(line_array, num_gauss_pts, kappa, bc_array)

    ! call print_string_matrix(ref_list)

    call surface_spline_refinement(input_nurbs_surface, ref_list, refined_surface)
    call assess_surface_refinement(input_nurbs_surface, refined_surface)

    call get_boundary_conditions_dof(refined_surface, bc_array, id_disp, u_pres)

    call assemble_weak_form(refined_surface, num_gauss_pts, kappa, Kmat, Fvec)
    call matrix_reduction(Kmat, Fvec, u_pres, id_disp, Kred, Fred, remainder_dofs)
    call solve_matrix_equations(Kred, Fred, remainder_dofs, id_disp, u_pres, Usol)
    call compute_postprocessing_solutions(Usol, refined_surface, file_name)
end program diffusion_example