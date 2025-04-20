program iga_1D
    use utils
    use isogeometric_beam_1D
    use reference_solutions
    implicit none

    real :: L, E, u0, tx, I
    integer :: numGaussPoints, p, num_control_points, num_knots, id_disp, id_load

    real, dimension(:), allocatable :: Uinit
    real, dimension(:), allocatable :: U_reduced
    real, dimension(:,:), allocatable :: control_points, weight_ctrl_pts, Kmat, Fvec
    real, dimension(:,:), allocatable :: Usol, solution_points, solution_weights, cpts, theo_vals
    

    L = 10.0
    E = 1000
    I = 1.
    u0 = 0.0
    tx = 25.0

    id_disp = 1

    numGaussPoints = 3

    ! Uinit = (/0.,0.,0.5,1.,1./)
    ! Uinit = (/0., 0., 0.25, 0.5, 0.75, 1., 1./)
    ! Uinit = (/0., 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1./)
    ! p = 1
    Uinit = (/0., 0., 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1., 1./)
    p = 2

    num_knots = size(Uinit)
    num_control_points = size(Uinit) - p - 1

    U_reduced = Uinit(p+1:num_knots-p)
    call print_row_vector(U_reduced)

    id_load = num_control_points

    allocate(control_points(num_control_points, 2))
    allocate(solution_points(num_control_points,2))
    allocate(solution_weights(num_control_points,1))
    allocate(weight_ctrl_pts(num_control_points,1))

    call linspace(0., L, num_control_points, 2, control_points)
    weight_ctrl_pts(:,1) = 1.0

    ! call print_matrix(control_points)

    call assemble_weak_form(control_points, weight_ctrl_pts, Uinit, U_reduced, p, numGaussPoints, E, I, tx, id_load, Kmat, Fvec)

    call solve_matrix_equations(Kmat, Fvec, id_disp, Usol)

    solution_points(:,1) = control_points(:,1)
    solution_points(:,2) = Usol(:,1)
    solution_weights(:,1) = 1.0

    call compute_field_solution(Uinit, p, solution_points, solution_weights, cpts)

    ! call print_matrix(cpts)
end program iga_1D