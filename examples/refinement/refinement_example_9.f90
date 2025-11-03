program refinement_example_9
    use utils
    use nurbs_surface
    use nurbs_curve
    use surface_refinement
    use input_output
    implicit none

    character(:), allocatable :: file_name
    character(len=50), dimension(:), allocatable :: line_array

    real, dimension(:,:), allocatable :: ctrl_pts, P_pts, w_pts
    real, dimension(:,:), allocatable :: Pref_pts, wref_pts
    real, dimension(:,:), allocatable :: spts_1, spts_2
    real, dimension(:), allocatable :: UP, VP, Uref, Vref

    integer :: p, q, pref, qref
    integer :: size_1, size_2, size_vec_1, size_vec_2
    
    character(len=1), dimension(:,:), allocatable :: ref_list

    ! character :: dir
    ! integer :: r, s, nu, nv

    file_name = "input_file_surface.txt"

    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, p, q, UP, VP, ctrl_pts, ref_list)

    size_1 = size(ctrl_pts, 1)
    size_2 = size(ctrl_pts, 2) - 1
    size_vec_1 = size(UP)
    size_vec_2 = size(VP)

    allocate(P_pts(0:size_1-1,0:size_2-1))
    allocate(w_pts(0:size_1-1,0))

    P_pts = ctrl_pts(1:size_1,1:size_2)
    w_pts = reshape(ctrl_pts(1:size_1,size_2+1), (/size_1, 1/))

    call print_matrix(P_pts)
    call print_matrix(w_pts)
    call print_string_matrix(ref_list)

    call surface_spline_refinement(p, q, P_pts, w_pts, UP, VP, ref_list, pref, qref, Pref_pts, wref_pts, Uref, Vref)

    call print_matrix(Pref_pts)
    call print_row_vector(Uref)
    call print_row_vector(Vref)

    call create_surface(p, q, P_pts, w_pts, UP, VP, spts_1)
    call create_surface(pref, qref, Pref_pts, wref_pts, Uref, Vref, spts_2)

    call assess_surface_refinement(spts_1, spts_2)
end program refinement_example_9