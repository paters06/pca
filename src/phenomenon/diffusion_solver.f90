module diffusion_solver
    use utils
    implicit none
contains

    function conductivity_matrix(kappa) result(kappa_mat)
        real, intent(in) :: kappa
        real, dimension(0:1,0:1) :: kappa_mat
        kappa_mat = 0.0
        kappa_mat(0,0) = kappa
        kappa_mat(1,1) = kappa
    end function conductivity_matrix

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

    subroutine compute_parametric_coords(ua, va, ub, vb, gauss_pts, gauss_weights, param_pts, param_weights)
        real, intent(in) :: ua, ub, va, vb
        real, dimension(:), intent(in) :: gauss_pts, gauss_weights
        real, dimension(:,:), allocatable, intent(out) :: param_pts
        real, dimension(:), allocatable, intent(out) :: param_weights

        integer :: i, j, i_param, num_gauss_pts

        num_gauss_pts = size(gauss_pts)

        allocate(param_pts(num_gauss_pts**2,0:1))
        allocate(param_weights(num_gauss_pts**2))

        param_pts = 0.0
        i_param = 1

        outer_loop: do i = 1, num_gauss_pts
            inner_loop: do j = 1, num_gauss_pts
                param_pts(i_param,0) = 0.5*(ub - ua)*gauss_pts(i) + 0.5*(ub + ua)
                param_pts(i_param,1) = 0.5*(vb - va)*gauss_pts(j) + 0.5*(vb + va)
                param_weights(i_param) = gauss_weights(i)*gauss_weights(j)
                i_param = i_param + 1
            end do inner_loop
        end do outer_loop
    end subroutine compute_parametric_coords

    subroutine derivative_shape_function_parametric(p, q, UP, VP, n_der, Pw_net, u, v, dR_dxi, non_zero_indices, global_dofs)
        use bspline_basis_functions, only: der_basis_functions, find_span
        integer, intent(in) :: p, q, n_der
        real, intent(in) :: u, v
        real, dimension(0:), intent(in) :: UP, VP
        real, dimension(0:,0:,0:), intent(in) :: Pw_net
        real, dimension(0:(p+1)*(q+1)-1,0:2), intent(out) :: dR_dxi
        integer, dimension(0:(p+1)*(q+1)-1,0:1), intent(out) :: non_zero_indices
        integer, dimension(0:(p+1)*(q+1)-1), intent(out) :: global_dofs

        integer :: i, j, rp, sq, nu, nv, uspan, vspan
        integer :: iu, jv, local_num
        real :: sum_W

        real, dimension(0:n_der,0:p) :: dN_dxi
        real, dimension(0:n_der,0:q) :: dM_deta
        
        real, dimension(0:(p+1)*(q+1)-1,0:1) :: dA_dxi
        real, dimension(0:1) :: dW_dxi

        rp = size(UP) - 1
        sq = size(VP) - 1
        nu = rp - p - 1
        nv = sq - q - 1

        call find_span(rp, p, u, UP, uspan)
        call find_span(sq, q, v, VP, vspan)
        call der_basis_functions(uspan, u, p, n_der, rp, UP, dN_dxi)
        call der_basis_functions(vspan, v, q, n_der, sq, VP, dM_deta)

        local_num = 0
        dR_dxi = 0.
        dA_dxi = 0.

        outer_loop: do j = 0, q
            inner_loop: do i = 0, p
                iu = uspan - p + i
                jv = vspan - q + j
                non_zero_indices(local_num, 0) = iu
                non_zero_indices(local_num, 1) = jv
                global_dofs(local_num) = compute_global_dof(iu,jv,nu)
                ! print "(A, I3, A, I3, A, I3, A, I3)", "iu:", iu, " |jv:", jv, " |nu:", nu, " |dof:", global_dofs(local_num)
                dR_dxi(local_num, 0) = dN_dxi(0,i)*dM_deta(0,j)*Pw_net(iu,jv,3)
                dA_dxi(local_num, 0) = dN_dxi(n_der,i)*dM_deta(0,j)*Pw_net(iu,jv,3)
                dA_dxi(local_num, 1) = dN_dxi(0,i)*dM_deta(n_der,j)*Pw_net(iu,jv,3)
                local_num = local_num + 1
            end do inner_loop
        end do outer_loop

        sum_W = sum(dR_dxi(:,0))
        dW_dxi = sum(dA_dxi,dim=1)

        dR_dxi(:,0) = dR_dxi(:,0)/sum_W
        dR_dxi(:,1) = (dA_dxi(:,0) - dW_dxi(0)*dR_dxi(:,0))/sum_W
        dR_dxi(:,2) = (dA_dxi(:,1) - dW_dxi(1)*dR_dxi(:,0))/sum_W
    end subroutine derivative_shape_function_parametric

    function gradient_parametric_physical(p, q, non_zero_indices, dR_dxi, Pw_net) result(dx_dxi)
        integer, intent(in) :: p, q
        integer, dimension(0:(p+1)*(q+1)-1,0:1), intent(in) :: non_zero_indices
        real, dimension(0:(p+1)*(q+1)-1,0:2), intent(in) :: dR_dxi
        real, dimension(0:,0:,0:), intent(in) :: Pw_net
        real, dimension(0:1,0:1) :: dx_dxi
        integer :: i, aa, bb, local_num, iu, jv

        ! dR_dxi(:,0) is the trivariate shape function
        ! dR_dxi(:,1) and dR_dxi(:,1) are the derivatives
        ! w.r.t U and V direction

        dx_dxi = 0.0

        local_num = size(non_zero_indices,1)

        non_zero_loop: do i = 0, local_num - 1
            iu = non_zero_indices(i, 0)
            jv = non_zero_indices(i, 1)
            aa_loop: do aa = 0, 1
                bb_loop: do bb = 0, 1
                    dx_dxi(aa,bb) = dx_dxi(aa,bb) + Pw_net(iu,jv,aa)*dR_dxi(i,bb+1)
                end do bb_loop
            end do aa_loop
        end do non_zero_loop
    end function

    function gradient_parent_parametric(ua, va, ub, vb, dx_dxi) result (J_mat)
        real, intent(in) :: ua, ub, va, vb
        real, dimension(0:1,0:1), intent(in) :: dx_dxi
        real, dimension(0:1,0:1) :: J_mat, dxi_dtildexi
        integer :: aa, bb, cc

        J_mat = 0.0
        dxi_dtildexi = 0.0

        dxi_dtildexi(0,0) = 0.5*(ub - ua)
        dxi_dtildexi(1,1) = 0.5*(vb - va)

        do aa = 0, 1
            do bb = 0, 1
                do cc = 0, 1
                    J_mat(aa,bb) = J_mat(aa,bb) + dx_dxi(aa,cc)*dxi_dtildexi(cc,bb)
                end do
            end do
        end do
    end function gradient_parent_parametric

    function derivative_basis_function_physical(p, q, dxi_dx, dR_dxi) result(dR_dx)
        integer, intent(in) :: p, q
        real, dimension(0:1,0:1), intent(in) :: dxi_dx
        real, dimension(0:(p+1)*(q+1)-1,0:2), intent(in) :: dR_dxi
        integer :: aa, bb, local_num, i
        real, dimension(0:(p+1)*(q+1)-1,0:1) :: dR_dx

        local_num = size(dR_dxi, 1) - 1

        dR_dx = 0.

        do i = 0, local_num
            aa_loop: do aa = 0, 1
                bb_loop: do bb = 0, 1
                    dR_dx(i,aa) = dR_dx(i,aa) + dxi_dx(bb,aa)*dR_dxi(i,bb+1)
                end do bb_loop
            end do aa_loop
        end do
    end function derivative_basis_function_physical

    subroutine element_stiffness_matrix(p, q, UP, VP, ua, va, ub, vb, n_der, kappa, Pw_net, &
                                        parametric_coordinates, parametric_weights, &
                                        K_local, global_dofs_loc)
        ! Adapted from Algorithm 1. Appendix 3.A Shape function routine from
        ! Isogeometric Analysis, Cottrell, Hughes and Bazilevs (2009)
        ! --------------------------------------------------------------------
        ! N_arr, M_arr, L_arr: Array of univariate B-spline basis functions
        ! dN_dxi, dM_deta, dL_dzeta: Univariate B-spline function derivatives
        !                            w.r.t appropiate parametric coordinates
        ! n_der: n-th derivative to compute
        ! p, q: degrees of the B-splines in the U and V direction
        use bspline_basis_functions, only: der_basis_functions, find_span
        integer, intent(in) :: p, q, n_der
        real, intent(in) :: ua, ub, va, vb, kappa
        real, dimension(0:), intent(in) :: UP, VP, parametric_weights
        real, dimension(0:,0:,0:), intent(in) :: Pw_net
        real, dimension(0:,0:), intent(in) :: parametric_coordinates
        real, dimension(0:(p+1)*(q+1)-1,0:(p+1)*(q+1)-1), intent(out) :: K_local
        integer, dimension(:,:), allocatable, intent(out) :: global_dofs_loc

        integer :: i_coor, num_param_coords, aa, bb, num_non_zero, i_global_loop
        real :: u, v, J_det, J_mod
        real, dimension(0:(p+1)*(q+1)-1,0:2) :: dR_dxi
        real, dimension(0:(p+1)*(q+1)-1,0:1) :: dR_dx
        real, dimension(0:1,0:1) :: dx_dxi, dxi_dx, J_mat, kappa_mat
        integer, dimension(0:(p+1)*(q+1)-1,0:1) :: non_zero_indices
        integer, dimension(0:(p+1)*(q+1)-1) :: global_dofs

        J_mat = 0.0
        J_mod = 0.0
        K_local = 0.0
        kappa_mat = conductivity_matrix(kappa)

        num_param_coords = size(parametric_coordinates,1)
        num_non_zero = (p+1)*(q+1)

        allocate(global_dofs_loc(0:(num_non_zero**2)-1,0:1))
        global_dofs_loc = 0

        param_coor_loop: do i_coor = 0, num_param_coords-1
            u = parametric_coordinates(i_coor, 0)
            v = parametric_coordinates(i_coor, 1)

            call derivative_shape_function_parametric(p, q, UP, VP, n_der, &
                                                      Pw_net, u, v, dR_dxi, &
                                                      non_zero_indices, global_dofs)
            dx_dxi = gradient_parametric_physical(p, q, non_zero_indices, dR_dxi, Pw_net)
            dxi_dx = inverse_matrix_2(dx_dxi)

            J_mat = gradient_parent_parametric(ua, va, ub, vb, dx_dxi)
            J_det = J_mat(0,0)*J_mat(1,1) - J_mat(1,0)*J_mat(0,1)
            J_mod = J_det*parametric_weights(i_coor)

            dR_dx = derivative_basis_function_physical(p, q, dxi_dx, dR_dxi)

            i_global_loop = 0

            row_loop: do bb = 0, num_non_zero - 1
                column_loop: do aa = 0, num_non_zero - 1
                    global_dofs_loc(i_global_loop,0) = global_dofs(aa)
                    global_dofs_loc(i_global_loop,1) = global_dofs(bb)
                    
                    K_local(aa,bb) = K_local(aa,bb) + (dR_dx(aa,0)*kappa_mat(0,0)*dR_dx(bb,0) &
                                    + dR_dx(aa,1)*kappa_mat(1,1)*dR_dx(bb,1))*J_mod
                    
                    i_global_loop = i_global_loop + 1
                end do column_loop
            end do row_loop

            ! temp_prod = matmul(kappa_mat,transpose(dR_dx))
            ! temp_prod_2 = matmul(dR_dx, temp_prod)*J_mod
            ! K_local = K_local + temp_prod_2
        end do param_coor_loop
    end subroutine element_stiffness_matrix
    
    subroutine assemble_weak_form(p, q, UP, VP, P_pts, w_pts, num_gauss_pts, kappa, Kmat, Fvec)
        use nurbs_curve_module, only: weighted_control_points
        use nurbs_surface_module, only: create_control_net
        integer, intent(in) :: p, q, num_gauss_pts
        real, intent(in) :: kappa
        real, dimension(0:), intent(in) :: UP, VP
        real, dimension(0:,0:), intent(in) :: P_pts, w_pts
        real, dimension(:,:), allocatable, intent(out) :: Kmat, Fvec

        real, dimension(:), allocatable :: gauss_nodes, gauss_weights, param_intg_weights
        real, dimension(:,:), allocatable :: surface_elem_bounds, param_intg_points
        real, dimension(:,:), allocatable :: Pw_pts
        real, dimension(:,:,:), allocatable :: Pw_net
        real, dimension(0:(p+1)*(q+1)-1,0:(p+1)*(q+1)-1) :: K_local
        integer, dimension(:,:), allocatable :: global_dofs_loc

        integer :: n_der
        integer :: i_elem, num_elements, r, s, nu, nv, num_dofs, i_global_loop
        integer :: i_global, j_global, num_non_zero, aa, bb
        real :: ua, ub, va, vb

        r = size(UP) - 1
        s = size(VP) - 1
        nu = r - p - 1
        nv = s - q - 1

        call weighted_control_points(P_pts, w_pts, Pw_pts)
        call create_control_net(nu, nv, Pw_pts, Pw_net)

        ! Degree of derivatives to compute
        n_der = 1

        call gauss_legendre_nodes_weights(num_gauss_pts, gauss_nodes, gauss_weights)

        call calculate_surface_elements(UP, VP, surface_elem_bounds)
        num_elements = size(surface_elem_bounds,1)

        num_dofs = size(P_pts,1)
        allocate(Kmat(0:num_dofs-1,0:num_dofs-1))
        allocate(Fvec(0:num_dofs-1,1))
        Kmat = 0.
        Fvec = 0.
        num_non_zero = (p+1)*(q+1)

        element_loop: do i_elem = 0, num_elements - 1
            K_local = 0.
        
            ua = surface_elem_bounds(i_elem, 0)
            ub = surface_elem_bounds(i_elem, 1)
            va = surface_elem_bounds(i_elem, 2)
            vb = surface_elem_bounds(i_elem, 3)
            
            call compute_parametric_coords(ua, ub, va, vb, gauss_nodes, gauss_weights, param_intg_points, param_intg_weights)
            call element_stiffness_matrix(p, q, UP, VP, ua, ub, va, vb, n_der, kappa, &
                                          Pw_net, param_intg_points, param_intg_weights, K_local, global_dofs_loc)

            i_global_loop = 0

            row_loop: do bb = 0, num_non_zero - 1
                column_loop: do aa = 0, num_non_zero - 1
                    i_global = global_dofs_loc(i_global_loop,0)
                    j_global = global_dofs_loc(i_global_loop,1)
                    Kmat(i_global, j_global) = Kmat(i_global, j_global) + K_local(aa,bb)
                    i_global_loop = i_global_loop + 1
                end do column_loop
            end do row_loop
        end do element_loop
    end subroutine assemble_weak_form

    subroutine matrix_reduction(Kmat, Fvec, u_pres, id_disp, Kred, Fred, remainder_dofs)
        use utils
        real, dimension(0:), intent(in) :: u_pres
        integer, dimension(0:), intent(in) :: id_disp
        real, dimension(:,:), allocatable, intent(in) :: Kmat, Fvec
        real, dimension(:,:), allocatable, intent(out) :: Kred
        real, dimension(:), allocatable, intent(out) :: Fred
        integer, dimension(:), allocatable, intent(out) :: remainder_dofs

        integer :: num_dofs, i, num_reduced_dofs
        integer, dimension(:), allocatable :: full_dofs
        logical, dimension(:), allocatable :: logical_array

        num_dofs = size(Kmat, 1)

        allocate(full_dofs(0:num_dofs-1))
        allocate(logical_array(0:num_dofs-1))
        full_dofs = (/(i, i=0, num_dofs-1)/)

        logical_array = (/(.false., i = 0, num_dofs-1)/)

        do i = 0, num_dofs-1
            if (any(full_dofs(i) == id_disp)) then
                logical_array(i) = .true.
            end if
        end do

        remainder_dofs = pack(full_dofs,.not. logical_array)

        num_reduced_dofs = size(remainder_dofs)
        allocate(Kred(0:num_reduced_dofs-1,0:num_reduced_dofs-1))
        allocate(Fred(0:num_reduced_dofs-1))

        Kred = Kmat(remainder_dofs,remainder_dofs)
        Fred = Fvec(remainder_dofs,1)

        do i = 0, size(id_disp)-1
            Fred = Fred - Kmat(remainder_dofs,id_disp(i))*u_pres(i)
        end do
    end subroutine matrix_reduction

    subroutine solve_matrix_equations(Kred, Fred, remainder_dofs, id_disp, u_pres, Usol)
        real, dimension(:,:), intent(in) :: Kred
        real, dimension(:), intent(in) :: Fred
        integer, dimension(:), intent(in) :: id_disp, remainder_dofs
        real, dimension(:), intent(in) :: u_pres
        real, dimension(:,:), intent(out), allocatable :: Usol
        
        real, dimension(:,:), allocatable :: Ured
        integer, dimension(:), allocatable :: pivot

        integer :: num_dofs_full, num_dofs_red, rc, i

        num_dofs_red = size(Kred,1)
        num_dofs_full = num_dofs_red + size(id_disp)

        allocate(pivot(num_dofs_red))
        allocate(Usol(0:num_dofs_full-1,1))

        Usol = 0.

        Ured = reshape(Fred, shape=(/size(Fred),1/))

        call sgesv(num_dofs_red, 1, Kred, num_dofs_red, pivot, Ured, num_dofs_red, rc)

        do i = 0, size(remainder_dofs)-1
            Usol(remainder_dofs(i+1),1) = Ured(i+1,1)
        end do
        
        do i = 0, size(id_disp)-1
            Usol(id_disp(i+1),1) = u_pres(i+1) 
        end do

        ! Usol(remainder_dofs,1) = Ured
        ! Usol(id_disp) = u_pres

        call print_column_vector(Usol)
    end subroutine solve_matrix_equations

    subroutine compute_postprocessing_solutions(p, q, F_pts, P_pts, w_pts, UP, VP, file_name)
        use nurbs_surface_module
        use input_output, only: export_matrix
        integer, intent(in) :: p, q
        real, dimension(:,:), allocatable, intent(in) :: P_pts, w_pts, F_pts
        real, dimension(:,:), allocatable :: spts, fpts, post_pts
        real, dimension(:), allocatable, intent(in) :: UP, VP
        integer :: num_points, n_dim_1, n_dim_2
        character(:), allocatable, intent(in) :: file_name
        character(:), allocatable :: file_output

        num_points = 25
        call create_surface(num_points, p, q, P_pts, w_pts, UP, VP, spts)
        call create_n_surface(num_points, p, q, F_pts, UP, VP, fpts)

        n_dim_1 = size(spts,2)
        n_dim_2 = size(fpts,2)

        allocate(post_pts((num_points+1)**2,n_dim_1+n_dim_2))
        post_pts(:,1:n_dim_1) = spts
        post_pts(:,n_dim_1+1:) = fpts

        file_output = 'output_'//file_name
        call export_matrix(post_pts, file_output)
    end subroutine compute_postprocessing_solutions
end module diffusion_solver