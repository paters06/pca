program isogeom_preprocessor
    use utils
    use nurbs_surface_module
    use surface_refinement
    use input_output
    use derived_types
    implicit none

    character(:), allocatable :: file_name
    character(len=50) :: file_name_inp
    character(len=50), dimension(:), allocatable :: line_array
    real, dimension(:,:), allocatable :: sbpts
    character(len=1), dimension(:,:), allocatable :: ref_list
    integer :: num_gauss_pts
    real :: kappa
    integer, dimension(:), allocatable :: id_disp
    real, dimension(:), allocatable :: u_pres
    type(boundary_condition), dimension(:), allocatable :: bc_array
    character(:), allocatable :: file_output_1, file_output_2, file_output_3
    type(nurbs_surface) :: input_nurbs_surface, refined_surface

    call get_command_argument(1, file_name_inp)
    file_name = file_name_inp

    ! Information import
    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, input_nurbs_surface, ref_list)
    call convert_data_to_solver(line_array, num_gauss_pts, kappa, bc_array)

    call print_nurbs_surface_info(input_nurbs_surface)
    call print_string_matrix(ref_list)

    call surface_spline_refinement(input_nurbs_surface, ref_list, refined_surface)
    call print_nurbs_surface_info(refined_surface)
    call assess_surface_refinement(input_nurbs_surface, refined_surface)

    call create_surface_boundary(refined_surface, sbpts)
    call get_boundary_conditions_dof(refined_surface, bc_array, id_disp, u_pres)

    file_output_1 = 'boundary_points_'//file_name
    file_output_2 = 'control_points_'//file_name
    file_output_3 = 'dirichlet_points_'//file_name
    call export_matrix(sbpts, file_output_1)
    call export_matrix(refined_surface%control_points, file_output_2)
    call export_matrix(refined_surface%control_points(id_disp+1,:), file_output_3)
end program isogeom_preprocessor