program iga_1D
    use utils
    use isogeometric_beam_1D
    implicit none

    real :: L, E, u0, tx, I
    integer :: numGaussPoints, p, num_control_points, num_knots, id_disp, id_load

    real, dimension(5) :: Uinit
    real, dimension(3) :: U_reduced
    real, dimension(:,:), allocatable :: control_points

    L = 10.0
    E = 1000
    I = 1.
    u0 = 0.0
    tx = 25.0

    id_disp = 1
    id_load = 3

    numGaussPoints = 3

    Uinit = (/0.,0.,0.5,1.,1./)
    p = 1

    num_knots = size(Uinit)
    num_control_points = size(Uinit) - p - 1

    U_reduced = Uinit(2:num_knots-1)

    allocate(control_points(num_control_points, 2))

    call linspace(0., L, num_control_points, 2, control_points)

    ! call print_matrix(control_points)

    call assemble_weak_form(Uinit, U_reduced, p, numGaussPoints, E, I, tx, id_load, id_disp)

end program iga_1D