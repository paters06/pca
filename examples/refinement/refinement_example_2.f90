program refinement_example_1
    use utils
    use curve_refinement
    use nurbs_curve
    implicit none

    real, dimension(:,:), allocatable :: P_pts, w_pts, Pw_pts, Qw_pts, cpts_1, cpts_2
    real, dimension(:,:), allocatable :: Q_pts, wq_pts
    real, dimension(:), allocatable :: U_knot, X_array, Ubar
    integer :: pinp, t
    real :: u
    integer :: size_1, size_2, size_vec

    size_1 = 3
    size_2 = 2
    size_vec = 6

    ! allocate(P_pts(3,2))
    allocate(P_pts(0:size_1-1,0:size_2-1))
    allocate(w_pts(0:size_1-1,0))
    allocate(U_knot(0:size_vec-1))
    P_pts = reshape((/(/1,1,0,0,1,1/)/),(/size_1,size_2/))
    ! P_pts = (/1,1,0,0,1,1/)
    w_pts = reshape((/1.,sqrt(2.),1./), (/size_1, 1/))
    pinp = 2
    U_knot = (/0,0,0,1,1,1/)

    call print_matrix(P_pts)

    ! X_array = (/0.5/)
    X_array = (/0.5, 0.75/)
    u = 0.5
    t = 1

    call weighted_control_points(P_pts, w_pts, Pw_pts)

    call degree_elevation(pinp, U_knot, Pw_pts, t, Qw_pts, Ubar)

    call geometric_control_points(Qw_pts, Q_pts, wq_pts)

    call print_row_vector(Ubar)
    call print_matrix(Qw_pts)
    call print_matrix(Q_pts)
    call print_matrix(wq_pts)

    call create_curve(P_pts, w_pts, U_knot, pinp, cpts_1)
    call create_curve(Q_pts, wq_pts, Ubar, pinp+t, cpts_2)

    call assess_refinement(cpts_1, cpts_2)

end program refinement_example_1