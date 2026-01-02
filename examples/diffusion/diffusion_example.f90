program diffusion_example
    use utils
    use utils_solver
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
    integer :: num_gauss_pts, refn_flag, interf_flag
    real :: kappa
    real, dimension(:,:), allocatable :: Kmat, Fvec, Kred, Usol
    real, dimension(:), allocatable :: Fred
    integer, dimension(:), allocatable :: id_disp, remainder_dofs, id_patches
    real, dimension(:,:), allocatable :: sbpts, ctrl_pts_pres
    real, dimension(:), allocatable :: u_pres
    type(boundary_condition), dimension(:), allocatable :: bc_array
    real, dimension(:,:), allocatable :: p_ctrl_pts
    type(nurbs_surface), dimension(:), allocatable :: input_patches
    type(interface_line), dimension(:), allocatable :: interf_var
    integer :: i_patch

    call get_command_argument(1, file_name_inp)
    file_name = file_name_inp

    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, input_patches, refn_flag, ref_list, interf_flag, interf_var)
    call convert_data_to_solver(line_array, num_gauss_pts, kappa, bc_array, id_patches)

    do i_patch = 1, size(id_patches)
        if (input_patches(i_patch)%refn_flag == 1) then
            call print_string_matrix(input_patches(i_patch)%refn_input)
            call surface_spline_refinement(input_patches(i_patch))
        else
            print "(A, I3)", "No refinement for patch #", i_patch
        end if
        call print_nurbs_surface_info(input_patches(i_patch))
    end do

    call create_patch_boundaries(input_patches, sbpts, p_ctrl_pts)
    call get_boundary_conditions_dof(input_patches, bc_array, id_patches, id_disp, u_pres, ctrl_pts_pres)

    call assemble_weak_form(input_patches, id_patches, num_gauss_pts, kappa, Kmat, Fvec)
    call matrix_reduction(Kmat, Fvec, u_pres, id_disp, Kred, Fred, remainder_dofs)
    call solve_matrix_equations(Kred, Fred, remainder_dofs, id_disp, u_pres, Usol)
    ! call compute_postprocessing_solutions(Usol, input_nurbs_surface, file_name)
end program diffusion_example