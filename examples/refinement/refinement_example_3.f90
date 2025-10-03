program refinement_example_3
    use utils
    use curve_refinement
    use nurbs_curve
    implicit none

    real, dimension(:,:), allocatable :: P_pts, w_pts, Ph_pts, wh_pts, cpts_1, cpts_2
    real, dimension(:), allocatable :: U_knot, Ubar
    integer :: pinp, ph

    integer :: size_1, size_2, size_vec

    character(len=1), dimension(:), allocatable :: ref_list

    size_1 = 3
    size_2 = 2
    size_vec = 6

    allocate(P_pts(0:size_1-1,0:size_2-1))
    allocate(w_pts(0:size_1-1,0))
    allocate(U_knot(0:size_vec-1))
    P_pts = reshape((/(/1,1,0,0,1,1/)/),(/size_1,size_2/))
    w_pts = reshape((/1.,sqrt(2.),1./), (/size_1, 1/))
    pinp = 2
    U_knot = (/0,0,0,1,1,1/)
    
    ! ref_list = (/'h','h'/)
    ref_list = (/'h', 'p', 'h'/)

    call print_matrix(P_pts)

    ! call h_refinement(pinp, P_pts, w_pts, U_knot, ph, Ph_pts, wh_pts, Ubar)
    ! call p_refinement(pinp, P_pts, w_pts, U_knot, ph, Ph_pts, wh_pts, Ubar)
    ! call k_refinement(pinp, P_pts, w_pts, U_knot, ph, Ph_pts, wh_pts, Ubar)

    call spline_refinement(pinp, P_pts, w_pts, U_knot, ref_list, ph, Ph_pts, wh_pts, Ubar)

    call print_row_vector(Ubar)
    call print_matrix(Ph_pts)
    call print_matrix(wh_pts)

    call create_curve(P_pts, w_pts, U_knot, pinp, cpts_1)
    call create_curve(Ph_pts, wh_pts, Ubar, ph, cpts_2)

    call assess_refinement(cpts_1, cpts_2)

end program refinement_example_3