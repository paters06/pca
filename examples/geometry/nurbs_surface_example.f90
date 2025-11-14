program nurbs_surface_example
    use utils
    use nurbs_surface_module
    implicit none

    integer :: p, q
    integer :: size_1, size_2, size_vec_1, size_vec_2, num_points
    real, dimension(:,:), allocatable :: P_pts
    real, dimension(:,:), allocatable :: w_pts
    real, dimension(:), allocatable :: U_knot, V_knot
    real, dimension(:,:), allocatable :: spts
    real, dimension(:,:,:,:), allocatable :: dspts
    ! character(:), allocatable :: file_name

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
    call print_row_vector(U_knot)
    call print_row_vector(V_knot)

    num_points = 41
    call create_surface(num_points, p, q, P_pts, w_pts, U_knot, V_knot, spts)
    call create_tangent_surface(p, q, P_pts, w_pts, U_knot, V_knot, dspts)

    ! call print_matrix(spts)
end program nurbs_surface_example