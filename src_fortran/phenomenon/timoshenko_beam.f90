module euler_bernoulli_beam
    implicit none
contains
    subroutine gauss_legendre_nodes_weights(num_gauss_pts, gauss_nodes, gauss_weights)
        ! Info from https://pomax.github.io/bezierinfo/legendre-gauss.html
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
        real, dimension(:), intent(out), allocatable :: intg_pts

        integer :: i, num_gauss_pts

        num_gauss_pts = size(gauss_pts)

        allocate(intg_pts(num_gauss_pts))

        do i = 1, num_gauss_pts
            intg_pts(i) = 0.5*(ub - ua)*gauss_pts(i) + 0.5*(ub + ua)
        end do
    end subroutine calculate_integration_points
    
    subroutine calculate_parametric_gradient(apt, cpt, detJac2)
        real, intent(in) :: apt, cpt
        real, intent(out) :: detJac2
        detJac2 = 0.5*(cpt - apt)
    end subroutine calculate_parametric_gradient
    
    subroutine calculate_1D_jacobian(Pw_array, U_array, p, d, u, jacob)
        use nurbs_curve
        use utils
        integer, intent(in) :: p, d
        real, intent(in) :: u
        real, dimension(:), intent(in) :: U_array
        real, dimension(:,:), intent(in) :: Pw_array
        real, dimension(2), intent(out) :: jacob

        real, dimension(:,:), allocatable :: CK
        real, dimension(2) :: tang_vec

        allocate(CK(d+1,2))

        call rat_curve_derivs(Pw_array, U_array, p, d, u, CK)
        call compute_tangent_vector(Pw_array, U_array, p, u, tang_vec)

        ! jacob = tang_vec
        jacob = CK(2,:)
    end subroutine calculate_1D_jacobian

    subroutine calculate_dR2_dx2(Pw_array, U_array, p, d, u, nders, dR2_dx2)
        ! Compute the second derivative of the shape functions
        ! with respect to the physical coordinates
        use nurbs_curve
        use utils
        real, dimension(:,:), intent(in) :: Pw_array, nders
        real, dimension(:), intent(in) :: U_array
        integer, intent(in) :: p, d
        real, intent(in) :: u
        real, dimension(p+1), intent(out) :: dR2_dx2
        real, dimension(:,:), allocatable :: CK
        real, dimension(p+1) :: dR_dx, dR_dxi, dR2_dxi2
        real, dimension(2) :: dx_dxi, dx2_dxi2
        real :: dx_dxi_val, dx2_dxi2_val

        allocate(CK(d+1,2))

        call rat_curve_derivs(Pw_array, U_array, p, d, u, CK)
        dx_dxi = CK(2,:)
        dx2_dxi2 = CK(3,:)

        dR_dxi = nders(2,:)
        dR2_dxi2 = nders(3,:)
        
        ! call print_row_vector(dR_dxi)
        ! call print_row_vector(dx_dxi)

        dx_dxi_val = norm2(dx_dxi)
        dx2_dxi2_val = norm2(dx2_dxi2)

        dR_dx = (1/dx_dxi_val)*dR_dxi

        dR2_dx2 = (1/(dx_dxi_val**2))*(dR2_dxi2 - dx2_dxi2_val*dR_dx)
    end subroutine calculate_dR2_dx2

    subroutine element_matrix_stiffness(Pw_array, U_array, p, apt, cpt, intg_pts, gauss_weights, &
                                        Kmat_i, global_dofs_k, global_dofs_l, length_i)
        ! id_K_elem_i is the ith integration point
        ! element_nonzero_K_values is the number of nonzero values in an element stiffness matrix
        use bspline_basis_functions
        use utils

        real, dimension(:), intent(in) :: U_array, intg_pts, gauss_weights
        real, dimension(:,:), intent(in) :: Pw_array
        real, intent(in) :: apt, cpt
        integer, intent(in) :: p
        real, dimension(:), allocatable :: Bmat_i
        real, dimension(:), intent(out), allocatable :: Kmat_i
        integer, dimension(:), intent(out), allocatable :: global_dofs_k, global_dofs_l
        real, intent(out) :: length_i

        integer :: span, m, j, l, k
        integer :: id_K_elem_i, element_nonzero_B_values, d, num_gauss_pts
        real :: u, detJac2, wJac
        real, dimension(:,:), allocatable :: nders
        real, dimension(p+1) :: dR2_dx2
        real, dimension(2) :: jacob

        num_gauss_pts = size(intg_pts)
        
        element_nonzero_B_values = (p+1)*(p+1)

        allocate(Kmat_i(2*element_nonzero_B_values))
        allocate(Bmat_i(2*element_nonzero_B_values))
        allocate(global_dofs_k(2*element_nonzero_B_values))
        allocate(global_dofs_l(2*element_nonzero_B_values))

        d = 1
        allocate(nders(d+1,p+1))

        m = size(U_array)

        length_i = 0.
        Kmat_i = 0.

        integration_loop: do j = 1, num_gauss_pts
            u = intg_pts(j)
            call find_span(m, p, u, U_array, span)
            call der_basis_functions(span, u, p, d, m, U_array, nders)
            call calculate_1D_jacobian(Pw_array, U_array, p, d, u, jacob)
            call calculate_parametric_gradient(apt, cpt, detJac2)
            call calculate_dR2_dx2(Pw_array, U_array, p, d, u, nders, dR2_dx2)
            ! call print_row_vector(dR2_dx2)
            wJac = norm2(jacob)*detJac2*gauss_weights(j)
            Bmat_i = 0.

            id_K_elem_i = 1
            row_loop: do l = 0, p
                column_loop: do k = 0, p
                    ! global_dofs_k(id_K_elem_i) = span-p+k+1
                    ! global_dofs_l(id_K_elem_i) = span-p+l+1
                    ! Bmat_i(id_K_elem_i) = dR2_dx2(k+1)*dR2_dx2(l+1)
                    global_dofs_k(2*id_K_elem_i-1) = 2*(span-p+k+1) - 1
                    global_dofs_k(2*id_K_elem_i) = 2*(span-p+k+1)
                    global_dofs_l(2*id_K_elem_i-1) = 2*(span-p+l+1) - 1
                    global_dofs_l(2*id_K_elem_i) = 2*(span-p+l+1)
                    B_be(2*id_K_elem_i) = dR_dx(k+1)*dR_dx(l+1)
                    B_se(2*id_K_elem_i - 1) = dR_dx(k+1)*dR_dx(l+1)
                    B_se(2*id_K_elem_i) = dR_dx(k+1)*dR_dx(l+1)
                    id_K_elem_i = id_K_elem_i + 1
                end do column_loop
            end do row_loop
            Kmat_i = Kmat_i + Bmat_i*wJac
            length_i = length_i + wJac
        end do integration_loop
        ! call print_row_vector(Kmat_i)
        ! call print_row_vector(Bmat_i)
        ! call print_row_vector_intg(global_dofs_k)
        ! call print_row_vector_intg(global_dofs_l)
    end subroutine element_matrix_stiffness

    subroutine element_force_vector()
    end subroutine element_force_vector

    subroutine assemble_weak_form(P_array, w_array, U_array, U_reduced, p, num_gauss_pts, young, inertia, load, id_load, Kmat, Fvec)
        ! U_array is the original knot vector. dimension(n)
        ! U_reduced is the knot vector with the repeated 0s and 1s. dimension(:)
        ! p is the degree of the spline. Integer
        ! num_gauss_pts is the number of integration points. Integer
        ! young is the young modulus. Real
        ! inertia is the second moment of area. Real
        ! load is the point load at the cantilever end. Real
        ! id_load is the dof for the load condition. Integer
        ! d is the 1st derivative of the splines. Integer
        use bspline_basis_functions
        use nurbs_curve
        use utils

        external :: sgesv
        real, dimension(:), intent(in) :: U_reduced, U_array
        real, dimension(:,:), intent(in) :: P_array, w_array
        integer, intent(in) :: num_gauss_pts, p, id_load
        real, intent(in) :: young, inertia, load
        real, dimension(:,:), intent(out), allocatable :: Kmat, Fvec

        integer :: num_unique_knots, i, j, num_control_points
        
        real, dimension(:), allocatable :: intg_pts, gauss_nodes, gauss_weights, Kmat_i
        real, dimension(:,:), allocatable :: Pw_array
        integer, dimension(:), allocatable :: global_dofs_k, global_dofs_l
        real :: total_length, length_i

        ! call print_matrix(P_array)
        allocate(Pw_array(size(P_array,1),size(P_array)+1))
        call weighted_control_points(P_array, w_array, Pw_array)
        ! call print_matrix(Pw_array)

        num_unique_knots = size(U_reduced)

        call gauss_legendre_nodes_weights(num_gauss_pts, gauss_nodes, gauss_weights)

        allocate(intg_pts(num_gauss_pts))

        num_control_points = size(U_array) - p - 1

        allocate(Kmat(2*num_control_points, 2*num_control_points))
        allocate(Kmat_i(2*num_gauss_pts*(p+1)*(p+1)))
        allocate(Fvec(2*num_control_points,1))

        total_length = 0.

        element_loop: do i = 1, num_unique_knots - 1
            call calculate_integration_points(gauss_nodes, U_reduced(i), U_reduced(i+1), intg_pts)
            call element_matrix_stiffness(Pw_array, U_array, p, U_reduced(i), U_reduced(i+1), &
                                          intg_pts, gauss_weights, &
                                          Kmat_i, global_dofs_k, global_dofs_l, length_i)
            ! print "(a, f9.3)", "length_i = ", length_i
            ! call print_row_vector_intg(global_dofs_k)
            ! call print_row_vector_intg(global_dofs_l)
            do j = 1, size(Kmat_i,1)
                Kmat(global_dofs_k(j), global_dofs_l(j)) = Kmat(global_dofs_k(j), global_dofs_l(j)) + &
                                                           young*inertia*Kmat_i(j)
            end do
            ! total_length = total_length + length_i
        end do element_loop
        ! call print_matrix(Kmat)

        Fvec(2*id_load-1,1) = -load
        ! call print_column_vector(Fvec)

        ! call element_force_vector(Fe)
        ! print "(a, f9.3)", "Total length = ", total_length
    end subroutine assemble_weak_form

    subroutine reduce_matrices(num_control_points, id_disp, Kmat, Fvec, Kred, Fred, remaining_dofs)
        use utils
        integer, intent(in) :: num_control_points
        real, dimension(:,:), intent(in) :: Kmat, Fvec
        real, dimension(:,:), intent(out) :: Kred, Fred
        integer, dimension(:), allocatable, intent(out) :: remaining_dofs
        integer, dimension(:), intent(in) :: id_disp
        integer, dimension(:), allocatable :: dofs_indices
        logical, dimension(:), allocatable :: temp_mask_array

        integer :: k_i, k_j, num_reduced_dofs

        call linspace_intg(1, num_control_points, 1, dofs_indices)
        ! call print_row_vector_intg(id_disp)
        ! call print_row_vector_intg(dofs_indices)

        num_reduced_dofs = num_control_points - size(id_disp)

        allocate(remaining_dofs(num_reduced_dofs))

        k_j = 1
        control_point_reduction_loop: do k_i = 1, num_control_points
            temp_mask_array = dofs_indices(k_i) /= id_disp
            if (all(temp_mask_array)) then
                remaining_dofs(k_j) = k_i
                k_j = k_j + 1
            end if
        end do control_point_reduction_loop

        ! call print_row_vector_intg(remaining_dofs)
        ! call print_matrix(Kmat)
        Kred = Kmat(remaining_dofs, remaining_dofs)
        Fred = Fvec(remaining_dofs, :)
        ! call print_matrix(Kred)

        ! call print_matrix(Kred)
        ! call print_column_vector(Fred)
    end subroutine reduce_matrices

    subroutine solve_matrix_equations(Kmat, Fvec, id_disp, Usol)
        use utils
        ! id_disp is the dof for the displacement condition. Integer

        real, dimension(:,:), intent(in) :: Kmat, Fvec
        real, dimension(:,:), intent(out), allocatable :: Usol
        integer, dimension(:), intent(in) :: id_disp

        real, dimension(:,:), allocatable :: Kred, Fred, Ured
        integer, dimension(:), allocatable :: pivot, remaining_dofs

        integer :: num_control_points, num_dofs_red, rc
        
        num_control_points = size(Kmat,1)
        num_dofs_red = num_control_points - size(id_disp)

        allocate(Kred(num_dofs_red, num_dofs_red))
        allocate(Fred(num_dofs_red,1))
        allocate(pivot(num_dofs_red))
        allocate(Usol(num_control_points,1))

        call reduce_matrices(num_control_points, id_disp, Kmat, Fvec, Kred, Fred, remaining_dofs)
        
        Ured = Fred

        call sgesv(num_dofs_red, 1, Kred, num_dofs_red, pivot, Ured, num_dofs_red, rc)

        Usol(remaining_dofs, :) = Ured
        Usol(id_disp, :) = 0.

        call print_column_vector(Usol)
    end subroutine solve_matrix_equations

    subroutine compute_field_solution(U_array, p, P_array, w_array, cpts)
        use nurbs_curve
        use utils
        real, dimension(:), intent(in) :: U_array
        real, dimension(:,:), intent(in) :: P_array, w_array
        integer, intent(in) :: p
        real, dimension(:,:), intent(out), allocatable :: cpts

        call create_curve(P_array, w_array, U_array, p, cpts)
    end subroutine compute_field_solution
end module euler_bernoulli_beam