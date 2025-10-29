program refinement_example_7
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
    ! real, dimension(:,:,:), allocatable :: Pw_net, Pwref_net

    integer :: p, q, pref, qref
    integer :: size_1, size_2, size_vec_1, size_vec_2
    
    character(len=1), dimension(:), allocatable :: ref_list

    character :: dir
    integer :: r, s, nu, nv

    file_name = "input_file_surface.txt"

    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, p, q, UP, VP, ctrl_pts, ref_list)

    size_1 = size(ctrl_pts, 1)
    size_2 = size(ctrl_pts, 2) - 1
    size_vec_1 = size(UP)
    size_vec_2 = size(VP)

    ! size_1 = 4
    ! size_2 = 3
    ! size_vec_1 = 4
    ! size_vec_2 = 4

    allocate(P_pts(0:size_1-1,0:size_2-1))
    allocate(w_pts(0:size_1-1,0))

    ! P_pts = reshape((/0., 0.2, 0. , 0.2, &
    !                   0., 0.0, 0.4, 0.4 , &
    !                   0., 0. , 0. , 0. /), (/size_1, size_2/))

    ! w_pts = reshape((/1.,1.,1.,1./), (/size_1, 1/))

    P_pts = ctrl_pts(1:size_1,1:size_2)
    w_pts = reshape(ctrl_pts(1:size_1,size_2+1), (/size_1, 1/))

    ! p = 1
    ! q = 1

    ! U_knot = (/0.,0.,1.,1./)
    ! V_knot = (/0.,0.,1.,1./)

    call print_matrix(P_pts)
    call print_matrix(w_pts)
    print *, ref_list

    dir = "U"

    r = size(UP) - 1
    s = size(VP) - 1
    nu = r - p - 1
    nv = s - q - 1

    ! call surface_knot_insertion(p, U_knot, q, V_knot, Pw_net, uv, dir, nq, UQ, mq, VQ, Qw_net)
    ! call surface_knot_refinement(p, UP, q, VP, Pw_net, X_array, dir, Uref, Vref, Pwref_net)
    call surface_h_refinement(p, q, P_pts, w_pts, UP, VP, dir, pref, Pref_pts, wref_pts, Uref, Vref)

    call print_matrix(Pref_pts)
    call print_row_vector(Uref)
    call print_row_vector(Vref)

    call create_surface(p, q, P_pts, w_pts, UP, VP, spts_1)
    call create_surface(p, q, Pref_pts, wref_pts, Uref, Vref, spts_2)

    call assess_surface_refinement(spts_1, spts_2)
end program refinement_example_7