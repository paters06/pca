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

    real, dimension(:,:), allocatable :: ctrl_pts, P_pts, w_pts
    real, dimension(:,:), allocatable :: Pref_pts, wref_pts
    real, dimension(:,:), allocatable :: spts_1, spts_2, sbpts
    real, dimension(:), allocatable :: UP, VP, Uref, Vref

    integer :: p, q, pref, qref
    integer :: size_1, size_2, size_vec_1, size_vec_2, num_points
    
    character(len=1), dimension(:,:), allocatable :: ref_list

    integer :: num_gauss_pts
    real :: kappa
    integer, dimension(:), allocatable :: id_disp
    real, dimension(:), allocatable :: u_pres

    type(boundary_condition), dimension(:), allocatable :: bc_array

    character(:), allocatable :: file_output_1, file_output_2, file_output_3

    call get_command_argument(1, file_name_inp)
    file_name = file_name_inp

    ! Information import
    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, p, q, UP, VP, ctrl_pts, ref_list)
    call convert_data_to_solver(line_array, num_gauss_pts, kappa, bc_array)

    size_1 = size(ctrl_pts, 1)
    size_2 = size(ctrl_pts, 2) - 1
    size_vec_1 = size(UP)
    size_vec_2 = size(VP)

    allocate(P_pts(0:size_1-1,0:size_2-1))
    allocate(w_pts(0:size_1-1,0))

    P_pts = ctrl_pts(1:size_1,1:size_2)
    w_pts = reshape(ctrl_pts(1:size_1,size_2+1), (/size_1, 1/))

    ! call print_matrix(P_pts)
    ! call print_matrix(w_pts)
    ! call print_string_matrix(ref_list)

    ! Surface refinement
    call surface_spline_refinement(p, q, P_pts, w_pts, UP, VP, ref_list, pref, qref, Pref_pts, wref_pts, Uref, Vref)

    call print_matrix(Pref_pts)
    call print_row_vector(Uref)
    call print_row_vector(Vref)

    num_points = 41
    call create_surface(num_points, p, q, P_pts, w_pts, UP, VP, spts_1)
    call create_surface(num_points, pref, qref, Pref_pts, wref_pts, Uref, Vref, spts_2)

    call assess_surface_refinement(spts_1, spts_2)

    call create_surface_boundary(pref, qref, Pref_pts, wref_pts, Uref, Vref, sbpts)
    
    call get_boundary_conditions_dof(pref, qref, Uref, Vref, bc_array, id_disp, u_pres)

    file_output_1 = 'boundary_points_'//file_name
    file_output_2 = 'control_points_'//file_name
    file_output_3 = 'dirichlet_points_'//file_name
    call export_matrix(sbpts, file_output_1)
    call export_matrix(Pref_pts, file_output_2)
    call export_matrix(Pref_pts(id_disp+1,:), file_output_3)
end program isogeom_preprocessor