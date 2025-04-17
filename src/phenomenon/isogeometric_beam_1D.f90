module isogeometric_beam_1D
    implicit none
contains
    subroutine gauss_legendre_nodes_weights(num_gauss_pts, gauss_nodes, gauss_weights)
        ! Info got from https://pomax.github.io/bezierinfo/legendre-gauss.html
        integer, intent(in) :: num_gauss_pts
        real, dimension(:), intent(out), allocatable :: gauss_nodes, gauss_weights

        allocate(gauss_nodes(num_gauss_pts))
        allocate(gauss_weights(num_gauss_pts))

        if (num_gauss_pts == 2) then
            gauss_nodes = (/-0.5774, 0.5774/)
            gauss_weights = (/1.0, 1.0/)
        else if (num_gauss_pts == 3) then
            gauss_nodes = (/-0.7746, 0.0, 0.7746/)
            gauss_weights = (/0.5556, 0.8889, 0.5556/)
        else if (num_gauss_pts == 4) then
            gauss_nodes = (/-0.8611, -0.3400, 0.3400, 0.8611/)
            gauss_weights = (/0.3479, 0.6521, 0.6521, 0.3479/)
        else if (num_gauss_pts == 5) then
            gauss_nodes = (/-0.9062, -0.5385, 0.0, 0.5385, 0.9062/)
            gauss_weights = (/0.2369, 0.4786, 0.5689, 0.4786, 0.2369/)
        else
            print *, "Value error"
        end if
    end subroutine gauss_legendre_nodes_weights
    
    subroutine calculate_integration_points(gauss_pts, ua, ub, intg_pts)
        real, intent(in) :: ua, ub
        real, dimension(:), intent(in) :: gauss_pts
        real, dimension(:), intent(out) :: intg_pts

        integer :: i, num_gauss_pts

        num_gauss_pts = size(gauss_pts)

        do i = 1, num_gauss_pts
            intg_pts(i) = 0.5*(ub - ua)*gauss_pts(i) + 0.5*(ub + ua)
        end do

    end subroutine calculate_integration_points
    
    subroutine element_matrix_stiffness(U_array, p, num_gauss_pts, intg_pts, young, inertia, Kmat_i, global_dofs_k, global_dofs_l)
        ! id_K_elem_i is the ith integration point
        ! element_nonzero_K_values is the number of nonzero values in an element stiffness matrix
        use bspline_basis_functions

        real, dimension(:), intent(in) :: U_array, intg_pts
        integer, intent(in) :: p, num_gauss_pts
        real, intent(in) :: young, inertia
        real, dimension(:), intent(out), allocatable :: Kmat_i
        integer, dimension(:), intent(out), allocatable :: global_dofs_k, global_dofs_l

        integer :: span, m, j, l, k, id_K_elem_i, element_nonzero_K_values, d
        real :: u
        real, dimension(:,:), allocatable :: nders

        element_nonzero_K_values = num_gauss_pts*(p+1)*(p+1)

        allocate(Kmat_i(element_nonzero_K_values))
        allocate(global_dofs_k(element_nonzero_K_values))
        allocate(global_dofs_l(element_nonzero_K_values))

        id_K_elem_i = 1

        d = 1
        allocate(nders(d+1,p+1))

        m = size(U_array)

        integration_loop: do j = 1, num_gauss_pts
            u = intg_pts(j)
            call find_span(m, p, u, U_array, span)
            call der_basis_functions(span, u, p, d, m, U_array, nders)
            row_loop: do l = 0, p
                column_loop: do k = 0, p
                    ! print "(a, I3, a, I3, a, I3, a, I3, a, I3)", "ele = ", i, " gauss_id", j," l = ", l, " k = ", k
                    global_dofs_k(id_K_elem_i) = span-p+k+1
                    global_dofs_l(id_K_elem_i) = span-p+l+1
                    Kmat_i(id_K_elem_i) = young*inertia*nders(p+1, k+1)*nders(p+1, l+1)
                    ! Kmat(span-p+k+1, span-p+l+1) = Kmat(span-p+k+1, span-p+l+1) + young*inertia*nders(p+1, k+1)*nders(p+1, l+1)
                    id_K_elem_i = id_K_elem_i + 1
                end do column_loop
            end do row_loop
        end do integration_loop
    end subroutine element_matrix_stiffness

    subroutine element_force_vector()
    end subroutine element_force_vector

    subroutine reduce_matrices(num_control_points, id_disp, Kmat, Fvec, Kred, Fred)
        use utils
        integer, intent(in) :: num_control_points, id_disp
        real, dimension(:,:), intent(in) :: Kmat, Fvec
        real, dimension(:,:), intent(out) :: Kred, Fred

        integer :: i_red, j_red, i, j

        j_red = 1
        reduction_loop_row: do j = 1, num_control_points
            if (j .ne. id_disp) then
                i_red = 1 
                reduction_loop_column: do i = 1, num_control_points
                    if (i .ne. id_disp) then
                        Kred(i_red, j_red) =  Kmat(i,j)
                        i_red = i_red + 1
                    end if
                end do reduction_loop_column
                
                Fred(j_red, 1) =  Fvec(j,1)
                j_red = j_red + 1
            end if
        end do reduction_loop_row

        call print_matrix(Kred)
        call print_column_vector(Fred)
    end subroutine reduce_matrices
    
    subroutine assemble_weak_form(U_array, U_reduced, p, num_gauss_pts, young, inertia, load, id_load, id_disp)
        ! U_array is the original knot vector. dimension(n)
        ! U_reduced is the knot vector with the repeated 0s and 1s. dimension(:)
        ! p is the degree of the spline. Integer
        ! num_gauss_pts is the number of integration points. Integer
        ! young is the young modulus. Real
        ! inertia is the second moment of area. Real
        ! load is the point load at the cantilever end. Real
        ! id_load is the dof for the load condition. Integer
        ! id_disp is the dof for the displacement condition. Integer
        ! d is the 1st derivative of the splines. Integer
        use bspline_basis_functions
        use utils

        external :: sgesv
        real, dimension(:), intent(in) :: U_reduced, U_array
        integer, intent(in) :: num_gauss_pts, p, id_load, id_disp
        real, intent(in) :: young, inertia, load

        integer :: num_unique_knots, i, j, num_control_points, num_dofs_red
        integer :: rc
        integer, dimension(:), allocatable :: pivot
        real, dimension(:), allocatable :: intg_pts, gauss_nodes, gauss_weights, Kmat_i
        real, dimension(:,:), allocatable :: Kmat, Fvec, Kred, Fred, Ured
        integer, dimension(:), allocatable :: global_dofs_k, global_dofs_l

        num_unique_knots = size(U_reduced)

        call gauss_legendre_nodes_weights(num_gauss_pts, gauss_nodes, gauss_weights)

        allocate(intg_pts(num_gauss_pts))

        num_control_points = size(U_array) - p - 1

        allocate(Kmat(num_control_points, num_control_points))
        allocate(Kmat_i(num_gauss_pts*(p+1)*(p+1)))
        allocate(Fvec(num_control_points,1))

        element_loop: do i = 1, num_unique_knots - 1
            call calculate_integration_points(gauss_nodes, U_reduced(i), U_reduced(i+1), intg_pts)
            call element_matrix_stiffness(U_array, p, num_gauss_pts, intg_pts, young, inertia, Kmat_i, global_dofs_k, global_dofs_l)
            do j = 1, size(Kmat_i,1)
                Kmat(global_dofs_k(j), global_dofs_l(j)) = Kmat(global_dofs_k(j), global_dofs_l(j)) + Kmat_i(j)
            end do
        end do element_loop
        ! call print_matrix(Kmat)

        Fvec(id_load,1) = load

        ! call print_column_vector(Fvec)

        num_dofs_red = num_control_points-1

        allocate(Kred(num_dofs_red, num_dofs_red))
        allocate(Fred(num_dofs_red,1))
        allocate(pivot(num_dofs_red))

        call reduce_matrices(num_control_points, id_disp, Kmat, Fvec, Kred, Fred)
        
        Ured = Fred

        call sgesv(num_dofs_red, 1, Kred, num_dofs_red, pivot, Ured, num_dofs_red, rc)

        ! call print_column_vector(Ured)

        ! call element_matrix_stiffness(Ke)
        ! call element_force_vector(Fe)
    end subroutine assemble_weak_form
end module isogeometric_beam_1D