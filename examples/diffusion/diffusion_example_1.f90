program diffusion_example_1
    use utils
    use nurbs_surface
    use surface_refinement
    use input_output
    use diffusion_solver
    implicit none

    character(:), allocatable :: file_name
    character(len=50) :: file_name_inp
    character(len=50), dimension(:), allocatable :: line_array

    real, dimension(:,:), allocatable :: ctrl_pts, P_pts, w_pts
    real, dimension(:,:), allocatable :: Pref_pts, wref_pts
    real, dimension(:,:), allocatable :: spts_1, spts_2
    real, dimension(:), allocatable :: UP, VP, Uref, Vref

    integer :: p, q, pref, qref
    integer :: size_1, size_2, size_vec_1, size_vec_2, num_points
    
    character(len=1), dimension(:,:), allocatable :: ref_list

    integer :: num_gauss_pts
    real :: kappa
    real, dimension(:,:), allocatable :: Kmat, Fvec, Kred, Usol
    real, dimension(:), allocatable :: Fred
    integer, dimension(:), allocatable :: id_disp, remainder_dofs
    real, dimension(:), allocatable :: u_pres

    ! Information import
    ! file_name = "input_file_diffusion.txt"
    write (*,*) "Write the name of the input file"
    write (*,*) "for the IGA solver"
    read (*,*) file_name_inp

    file_name = file_name_inp

    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, p, q, UP, VP, ctrl_pts, ref_list)
    call convert_data_to_solver(line_array, num_gauss_pts, kappa, id_disp, u_pres)

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

    ! IGA implementation
    call assemble_weak_form(pref, qref, Uref, Vref, Pref_pts, wref_pts, num_gauss_pts, kappa, Kmat, Fvec)
    call matrix_reduction(Kmat, Fvec, u_pres, id_disp, Kred, Fred, remainder_dofs)
    call solve_matrix_equations(Kred, Fred, remainder_dofs, id_disp, u_pres, Usol)
    call compute_postprocessing_solutions(pref, qref, Usol, Pref_pts, wref_pts, Uref, Vref, file_name)
end program diffusion_example_1