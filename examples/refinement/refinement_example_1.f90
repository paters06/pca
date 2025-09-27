program refinement_example_1
    use utils
    use curve_refinement
    use nurbs_curve
    implicit none

    real, dimension(:,:), allocatable :: P_pts, w_pts, Pw_pts, Qw_pts, cpts_1, cpts_2
    real, dimension(:,:), allocatable :: Q_pts, wq_pts
    real, dimension(:), allocatable :: U_knot, X_array, Ubar
    integer :: pinp
    real :: u

    allocate(P_pts(3,2))
    P_pts = reshape((/(/1,1,0,0,1,1/)/),(/size(P_pts,1),size(P_pts,2)/))
    w_pts = reshape((/1.,sqrt(2.),1./), (/size(P_pts,1), 1/))
    pinp = 2
    U_knot = (/0,0,0,1,1,1/)

    call print_matrix(P_pts)

    X_array = (/0.5/)
    u = 0.5

    call weighted_control_points(P_pts, w_pts, Pw_pts)

    ! call knot_refinement(pinp, U_knot, X_array, Pw_pts, Qw_pts, Ubar)
    call knot_insertion(pinp, U_knot, Pw_pts, u, Ubar, Qw_pts)

    call geometric_control_points(Qw_pts, Q_pts, wq_pts)

    call print_matrix(Qw_pts)
    call print_matrix(Q_pts)
    call print_matrix(wq_pts)

    call create_curve(P_pts, w_pts, U_knot, pinp, cpts_1)
    call create_curve(Q_pts, wq_pts, Ubar, pinp, cpts_2)

    call assess_refinement(cpts_1, cpts_2)

end program refinement_example_1