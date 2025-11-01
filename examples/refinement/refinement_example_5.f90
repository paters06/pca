program refinement_example_5
    use utils
    use nurbs_surface
    use nurbs_curve
    use surface_refinement
    implicit none

    integer :: p, q
    integer :: size_1, size_2, size_vec_1, size_vec_2
    real, dimension(:,:), allocatable :: P_pts, Pw, Qw, Q_pts
    real, dimension(:,:), allocatable :: w_pts, wq_pts
    real, dimension(:,:,:), allocatable :: Pw_net, Qw_net
    real, dimension(:), allocatable :: U_knot, V_knot, UQ, VQ
    real, dimension(:,:), allocatable :: spts, spts_2
    ! character(:), allocatable :: file_name

    real :: uv
    character :: dir
    integer :: r, s, nu, nv, nq, mq

    size_1 = 4
    size_2 = 3
    size_vec_1 = 4
    size_vec_2 = 4

    allocate(P_pts(0:size_1-1,0:size_2-1))
    allocate(w_pts(0:size_1-1,0))
    allocate(U_knot(0:size_vec_1-1))
    allocate(V_knot(0:size_vec_1-1))

    P_pts = reshape((/0., 0.2, 0. , 0.2, &
                      0., 0.0, 0.4, 0.4 , &
                      0., 0. , 0. , 0. /), (/size_1, size_2/))

    w_pts = reshape((/1.,1.,1.,1./), (/size_1, 1/))

    p = 1
    q = 1

    U_knot = (/0.,0.,1.,1./)
    V_knot = (/0.,0.,1.,1./)

    call print_matrix(P_pts)

    uv = 0.5
    dir = "V"

    r = size(U_knot) - 1
    s = size(V_knot) - 1
    nu = r - p - 1
    nv = s - q - 1

    call weighted_control_points(P_pts, w_pts, Pw)
    call create_control_net(nu,nv,Pw,Pw_net)
    
    call surface_knot_insertion(p, U_knot, q, V_knot, Pw_net, uv, dir, nq, UQ, mq, VQ, Qw_net)

    call create_control_list(Qw_net, Qw)
    call geometric_control_points(Qw, Q_pts, wq_pts)

    call print_matrix(Q_pts)
    call print_row_vector(UQ)
    call print_row_vector(VQ)

    call create_surface(p, q, P_pts, w_pts, U_knot, V_knot, spts)
    call create_surface(p, q, Q_pts, wq_pts, UQ, VQ, spts_2)

    call assess_surface_refinement(spts, spts_2)
end program refinement_example_5