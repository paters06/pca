program iga_1D
    use utils
    use euler_bernoulli_beam
    use reference_solutions
    implicit none

    real :: L, E, u0, tx, I, t, b
    integer :: numGaussPoints, p, num_control_points, num_knots, id_load

    real, dimension(:), allocatable :: Uinit
    real, dimension(:), allocatable :: U_reduced
    real, dimension(:), allocatable :: deflection_ctrl_pts
    real, dimension(:,:), allocatable :: control_points, weight_ctrl_pts, Kmat, Fvec
    real, dimension(:,:), allocatable :: Usol, solution_points, solution_weights, cpts
    integer, dimension(:), allocatable :: id_disp
    ! integer, dimension(5,2) :: array

    L = 10.0
    E = 100.
    t = 1.
    b = 1.
    I = (1/12.)*b*t**3
    u0 = 0.0
    tx = 9000

    id_disp = (/1,2/)

    numGaussPoints = 3

    ! Uinit = (/0.,0.,0.,0.,1.,1.,1.,1./)*L
    ! Uinit = (/0.,0.,0.,0.,0.5,1.,1.,1.,1./)*L
    ! Uinit = (/0.,0.,0.,0.,0.25,0.5,0.75,1.,1.,1.,1./)*L
    Uinit = (/0.,0.,0.,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.,1.,1./)*L
    p = 3 
    ! Uinit = (/0., 0., 0., 0.25, 0.5, 0.75, 1., 1., 1./)*L
    ! p = 2
    ! Uinit = (/0., 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1./)*L
    ! p = 1
    ! Uinit = (/0., 0., 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1., 1./)*L
    ! p = 2
    ! Uinit = (/0.,0.,0.,0.5,1.,1.,1./)*L
    ! p = 2

    num_knots = size(Uinit)
    num_control_points = size(Uinit) - p - 1

    U_reduced = Uinit(p+1:num_knots-p)

    id_load = num_control_points

    allocate(control_points(num_control_points, 2))
    allocate(solution_points(num_control_points, 2))
    allocate(solution_weights(num_control_points,1))
    allocate(weight_ctrl_pts(num_control_points,1))

    call linspace_real(0., L, num_control_points, 2, control_points)
    weight_ctrl_pts(:,1) = 1.0
    ! control_points=reshape((/1.,1.,0.,0.,1.,1./),shape(array))
    ! weight_ctrl_pts(:,1) = (/1.,1.,2./)
    ! control_points=reshape((/0., 1., 3., 4., 5., 0., 1., 2., 1., -1./),shape(array))
    ! weight_ctrl_pts(:,1) = (/1., 4., 1., 1., 1./)

    ! call print_row_vector(U_reduced)
    ! call print_matrix(control_points)

    call assemble_weak_form(control_points, weight_ctrl_pts, Uinit, U_reduced, p, numGaussPoints, E, I, tx, id_load, Kmat, Fvec)

    call solve_matrix_equations(Kmat, Fvec, id_disp, Usol)

    allocate(deflection_ctrl_pts(num_control_points))
    deflection_ctrl_pts = Usol(1:2*num_control_points:2,1)

    solution_points(:,1) = control_points(:,1)
    solution_points(:,2) = deflection_ctrl_pts
    solution_weights(:,1) = 1.0

    ! call print_matrix(solution_points)

    call compute_field_solution(Uinit, p, solution_points, solution_weights, cpts)

    call print_matrix(cpts)
end program iga_1D