module nurbs_surface_module
    implicit none
contains
    subroutine create_control_net(nu, nv, P_pts, P_net)
        use utils
        integer, intent(in) :: nu, nv
        real, dimension(:,:), allocatable, intent(in) :: P_pts
        real, dimension(:,:,:), allocatable, intent(out) :: P_net
        
        integer :: i

        allocate(P_net(0:nu,0:nv,0:3))
        P_net = 0.0

        do i = 0, 3
            P_net(:,:,i) = reshape(P_pts(:,i), (/nu+1,nv+1/))
        end do
    end subroutine create_control_net

    subroutine create_n_control_net(nu, nv, n_dim, P_pts, P_net)
        use utils
        integer, intent(in) :: nu, nv, n_dim
        real, dimension(:,:), allocatable, intent(in) :: P_pts
        real, dimension(:,:,:), allocatable, intent(out) :: P_net
        
        integer :: i

        allocate(P_net(0:nu,0:nv,0:n_dim))
        P_net = 0.0

        do i = 0, n_dim
            P_net(:,:,i) = reshape(P_pts(:,i), (/nu+1,nv+1/))
        end do
    end subroutine create_n_control_net

    subroutine create_control_list(P_net, P_pts)
        real, dimension(:,:,:), allocatable, intent(in) :: P_net
        real, dimension(:,:), allocatable, intent(out) :: P_pts
        
        integer :: i, nu, nv, num_ctrl_pts

        nu = size(P_net,1)
        nv = size(P_net,2)

        num_ctrl_pts = nu*nv

        allocate(P_pts(0:num_ctrl_pts-1,0:3))

        do i = 0, 3
            P_pts(:,i) = reshape(P_net(:,:,i), (/num_ctrl_pts/))
        end do
    end subroutine create_control_list

    subroutine surface_point(r, s, p, q, u, v, Pw_net, U_array, V_array, S_pt)
        ! Algorithm A4.3 from the NURBS Book
        ! Compute ponit on rational B-spline surface
        ! Input: n, p, U, m, q, V, Pw, u, v
        ! Output: S
        ! n+1 and m+1 are the number of control points in the U and V direction
        ! r+1 and s+1 are the number of elements of the U and V knot vectors
        ! p and q are the orders in the U and V knot vectors
        use utils
        use bspline_basis_functions, only: find_span, basis_function

        integer, intent(in) :: r, s, p, q
        real, intent(in) :: u, v
        real, dimension(0:r), intent(in) :: U_array
        real, dimension(0:s), intent(in) :: V_array
        real, dimension(:,:,:), allocatable, intent(in) :: Pw_net
        real, dimension(0:2), intent(out) :: S_pt

        integer :: uspan, vspan, l, k
        real, dimension(0:p) :: NU_arr
        real, dimension(0:q) :: NV_arr
        real, dimension(0:q,0:3) :: temp
        real, dimension(0:3) :: Sw

        temp = 0.0

        call find_span(r, p, u, U_array, uspan)
        call basis_function(uspan, u, r, p, U_array, NU_arr)
        call find_span(s, q, v, V_array, vspan)
        call basis_function(vspan, v, s, q, V_array, NV_arr)

        outer_loop: do l = 0, q
            temp(l,:) = 0.0
            inner_loop: do k = 0, p
                temp(l,:) = temp(l,:) + NU_arr(k)*Pw_net(uspan-p+k,vspan-q+l,:)
            end do inner_loop
        end do outer_loop

        S_pt = 0.0
        Sw = 0.0
        do l = 0, q
            Sw = Sw + NV_arr(l)*temp(l,:)
        end do

        S_pt = Sw(0:2)/Sw(3)
    end subroutine surface_point

    subroutine n_surface_point(r, s, p, q, u, v, n_dim, Pw_net, U_array, V_array, S_pt)
        ! Algorithm A4.3 from the NURBS Book
        ! Compute ponit on rational B-spline surface
        ! Input: n, p, U, m, q, V, Pw, u, v
        ! Output: S
        ! n+1 and m+1 are the number of control points in the U and V direction
        ! r+1 and s+1 are the number of elements of the U and V knot vectors
        ! p and q are the orders in the U and V knot vectors
        use utils
        use bspline_basis_functions, only: find_span, basis_function

        integer, intent(in) :: r, s, p, q, n_dim
        real, intent(in) :: u, v
        real, dimension(0:r), intent(in) :: U_array
        real, dimension(0:s), intent(in) :: V_array
        real, dimension(:,:,:), allocatable, intent(in) :: Pw_net
        real, dimension(0:n_dim-1), intent(out) :: S_pt

        integer :: uspan, vspan, l, k
        real, dimension(0:p) :: NU_arr
        real, dimension(0:q) :: NV_arr
        real, dimension(0:q,0:n_dim) :: temp
        real, dimension(0:n_dim) :: Sw

        temp = 0.0

        call find_span(r, p, u, U_array, uspan)
        call basis_function(uspan, u, r, p, U_array, NU_arr)
        call find_span(s, q, v, V_array, vspan)
        call basis_function(vspan, v, s, q, V_array, NV_arr)

        outer_loop: do l = 0, q
            temp(l,:) = 0.0
            inner_loop: do k = 0, p
                ! print *, Pw_net(uspan-p+k,vspan-q+l,:)
                temp(l,:) = temp(l,:) + NU_arr(k)*Pw_net(uspan-p+k,vspan-q+l,:)
            end do inner_loop
        end do outer_loop

        S_pt = 0.0
        Sw = 0.0
        do l = 0, q
            Sw = Sw + NV_arr(l)*temp(l,:)
        end do

        S_pt = Sw(0:n_dim-1)/Sw(n_dim)
    end subroutine n_surface_point

    subroutine create_surface(num_points, surf, spts)
        ! r+1 is the length of the U knot vector
        ! s+1 is the length of the V knot vector
        ! n+1 is the number of control points
        use utils
        use nurbs_curve_module, only: weighted_control_points
        use derived_types
        
        type(nurbs_surface), intent(in) :: surf
        integer, intent(in) :: num_points
        integer :: p, q
        real, dimension(:,:), allocatable :: P_pts, w_pts
        real, dimension(:), allocatable :: U_array, V_array
        real, intent(out), dimension(:,:), allocatable :: spts

        real :: U_min, U_max, V_min, V_max, u, v
        integer :: i, j, nu, nv, r, s, i_pt
        real, dimension(:,:), allocatable :: Pw
        real, dimension(:), allocatable :: S_pti
        real, dimension(:,:,:), allocatable :: Pw_net

        p = surf%p
        q = surf%q
        P_pts = surf%control_points
        w_pts = surf%weight_points
        U_array = surf%U_knot
        V_array = surf%V_knot

        r = size(U_array) - 1
        s = size(V_array) - 1
        nu = r - p - 1
        nv = s - q - 1
        U_max = maxval(U_array)
        U_min = minval(U_array)
        V_max = maxval(V_array)
        V_min = minval(V_array)

        call weighted_control_points(P_pts, w_pts, Pw)
        call create_control_net(nu,nv,Pw,Pw_net)
        
        allocate(spts((num_points+1)**2,3))
        allocate(S_pti(0:2))
        spts = 0.0
        S_pti = 0.0

        i_pt = 1

        do j = 1, num_points+1
            do i = 1, num_points+1
                u = ((U_max - U_min)/num_points)*(i-1) + U_min
                v = ((V_max - V_min)/num_points)*(j-1) + V_min
                call surface_point(r, s, p, q, u, v, Pw_net, U_array, V_array, S_pti)
                spts(i_pt,:) = S_pti(:)
                i_pt = i_pt + 1
            end do
        end do
    end subroutine create_surface

    subroutine create_n_surface(num_points, surf, F_pts, spts)
        ! r+1 is the length of the U knot vector
        ! s+1 is the length of the V knot vector
        ! n+1 is the number of control points
        ! F_pts are the control points of the Field of interest
        use utils
        use nurbs_curve_module, only: weighted_control_points
        use derived_types
        type(nurbs_surface), intent(in) :: surf
        integer :: p, q
        integer, intent(in) :: num_points
        real, intent(in), allocatable, dimension(:,:) :: F_pts
        real, dimension(:), allocatable :: U_array, V_array
        real, intent(out), dimension(:,:), allocatable :: spts

        real :: U_min, U_max, V_min, V_max, u, v
        integer :: i, j, nu, nv, r, s, n_dim, i_pt
        real, dimension(:), allocatable :: S_pti, wf_pts
        real, dimension(:,:), allocatable :: wf_col, Fw
        real, dimension(:,:,:), allocatable :: Pw_net

        p = surf%p
        q = surf%q
        U_array = surf%U_knot
        V_array = surf%V_knot

        ! num_points = 625
        r = size(U_array) - 1
        s = size(V_array) - 1
        nu = r - p - 1
        nv = s - q - 1
        U_max = maxval(U_array)
        U_min = minval(U_array)
        V_max = maxval(V_array)
        V_min = minval(V_array)

        ! The first dimension starts with index zero
        n_dim = size(F_pts,2) - 1
        
        allocate(wf_pts(size(F_pts,1)))
        wf_pts = 1.0

        wf_col = reshape(wf_pts, (/size(F_pts,1), 1/))

        call weighted_control_points(F_pts, wf_col, Fw)
        call create_n_control_net(nu,nv,n_dim+1,Fw,Pw_net)

        allocate(spts((num_points+1)**2,n_dim+1))
        allocate(S_pti(0:n_dim))
        spts = 0.0
        S_pti = 0.0

        i_pt = 1

        do j = 1, num_points+1
            do i = 1, num_points+1
                u = ((U_max - U_min)/num_points)*(i-1) + U_min
                v = ((V_max - V_min)/num_points)*(j-1) + V_min
                call n_surface_point(r, s, p, q, u, v, n_dim+1, Pw_net, U_array, V_array, S_pti)
                spts(i_pt,:) = S_pti(:)
                i_pt = i_pt + 1
            end do
        end do
    end subroutine create_n_surface

    subroutine bivariate_rational_function
    end subroutine bivariate_rational_function

    subroutine create_tangent_surface(p, q, P_pts, w_pts, U_array, V_array, spts)
        ! r+1 is the length of the U knot vector
        ! s+1 is the length of the V knot vector
        ! n+1 is the number of control points
        use utils
        use nurbs_curve_module, only: weighted_control_points
        
        integer, intent(in) :: p, q
        real, intent(in), dimension(:,:) :: P_pts
        real, intent(in), dimension(:,:) :: w_pts
        real, intent(in), dimension(:) :: U_array, V_array
        real, intent(out), dimension(:,:,:,:), allocatable :: spts

        real :: U_min, U_max, V_min, V_max, u, v
        integer :: num_points, i, nu, nv, r, s, d, k, l
        real, dimension(:,:), allocatable :: Pw
        real, dimension(:,:,:), allocatable :: SwKL, SKLi
        real, dimension(:,:,:), allocatable :: Aders
        real, dimension(:,:), allocatable :: wders
        real, dimension(:,:,:), allocatable :: Pw_net

        num_points = 41
        r = size(U_array) - 1
        s = size(V_array) - 1
        nu = r - p - 1
        nv = s - q - 1
        U_max = maxval(U_array)
        U_min = minval(U_array)
        V_max = maxval(V_array)
        V_min = minval(V_array)

        call weighted_control_points(P_pts, w_pts, Pw)

        call create_control_net(nu,nv,Pw,Pw_net)
        
        d = 1
        k = 1
        l = d - k

        allocate(spts(num_points+1,0:d,0:d,3))
        allocate(SwKL(0:d,0:d,0:3))
        allocate(SKLi(0:d,0:d,0:2))
        spts = 0.0
        SwKL = 0.0

        do i = 1, num_points+1
            u = ((U_max - U_min)/num_points)*(i-1) + U_min
            v = ((V_max - V_min)/num_points)*(i-1) + V_min
            call bspline_surface_gradient(r, p, U_array, s, q, V_array, Pw_net, u, v, d, SwKL)
            Aders = SwKL(:,:,0:2)
            wders = SwKL(:,:,3)
            call bivariate_rational_gradient(Aders, wders, d, SKLi)
            spts(i,:,:,:) = SKLi(:,:,:)
        end do
    end subroutine create_tangent_surface

    subroutine create_surface_boundary(surf, sbpts)
        ! r+1 is the length of the U knot vector
        ! s+1 is the length of the V knot vector
        ! n+1 is the number of control points
        use utils
        use nurbs_curve_module, only: weighted_control_points
        use derived_types
        
        type(nurbs_surface), intent(in) :: surf
        integer:: p, q
        real, dimension(:,:), allocatable :: P_pts
        real, dimension(:,:), allocatable :: w_pts
        real, dimension(:), allocatable :: U_array, V_array
        real, intent(out), dimension(:,:), allocatable :: sbpts

        real :: u, v
        integer :: j, nu, nv, r, s, num_boundary_points
        real, dimension(:,:), allocatable :: Pw, bound_param_pts
        real, dimension(:), allocatable :: S_pti
        real, dimension(:,:,:), allocatable :: Pw_net

        real, dimension(0:3) :: first_boundary, second_boundary
        real, dimension(0:3) :: third_boundary, fourth_boundary

        p = surf%p
        q = surf%q
        P_pts = surf%control_points
        w_pts = surf%weight_points
        U_array = surf%U_knot
        V_array = surf%V_knot

        r = size(U_array) - 1
        s = size(V_array) - 1
        nu = r - p - 1
        nv = s - q - 1

        call weighted_control_points(P_pts, w_pts, Pw)
        call create_control_net(nu,nv,Pw,Pw_net)

        first_boundary = (/0.0,0.0,1.0,0.0/)
        second_boundary = (/1.0,0.0,1.0,1.0/)
        third_boundary = (/1.0,1.0,0.0,1.0/)
        fourth_boundary = (/0.0,1.0,0.0,0.0/)

        bound_param_pts = compute_parametric_boundary_points(first_boundary, &
                          second_boundary, third_boundary, fourth_boundary)
        
        num_boundary_points = size(bound_param_pts,1)

        allocate(sbpts(0:num_boundary_points-1,0:2))
        allocate(S_pti(0:2))
        sbpts = 0.0
        S_pti = 0.0

        do j = 1, num_boundary_points
            u = bound_param_pts(j, 1)
            v = bound_param_pts(j, 2)
            call surface_point(r, s, p, q, u, v, Pw_net, U_array, V_array, S_pti)
            sbpts(j-1,:) = S_pti(:)
        end do
    end subroutine create_surface_boundary

    subroutine bspline_surface_gradient(r, p, U_array, s, q, V_array, Pw_net, u, v, d, SKL)
        ! Algorithm 3.6 from the NURBS Book
        ! Compute surface point and derivatives
        ! Input: n, p, U, m, q, V, P, u, v, d
        ! Output: SKL
        ! k is the k-th derivative with respect to u
        ! l is the l-th derivative with respect to v
        use bspline_basis_functions, only: find_span, der_basis_functions

        integer, intent(in) :: r, p, s, q, d
        real, intent(in), dimension(0:r) :: U_array
        real, intent(in), dimension(0:s) :: V_array
        real, intent(in) :: u, v
        real, dimension(0:,0:,0:), intent(in) :: Pw_net
        real, intent(out), dimension(:,:,:), allocatable :: SKL
        
        integer :: du, k, l, dv, uspan, vspan, dd, ri, si
        real, dimension(0:d,0:p) :: NU_ders
        real, dimension(0:d,0:q) :: NV_ders
        real, dimension(0:p,0:3) :: temp

        allocate(SKL(0:d,0:d,0:3))

        du = min(d, p)
        do k = p+1, d
            do l = 0, d-k
                SKL(k,l,:) = 0.0
            end do
        end do

        dv = min(d,q)
        do l = q+1, d
            do k = 0, d-l
                SKL(k,l,:) = 0.0
            end do
        end do

        call find_span(r, p, u, U_array, uspan)
        call der_basis_functions(uspan, u, p, du, r, U_array, NU_ders)
        call find_span(s, q, v, V_array, vspan)
        call der_basis_functions(vspan, v, q, dv, s, V_array, NV_ders)

        do k = 0, du
            do si = 0, q
                temp(si,:) = 0.0
                do ri = 0, p
                    temp(si,:) = temp(si,:) + NU_ders(k,ri)*Pw_net(uspan-p+ri,vspan-q+si,:)
                end do
            end do

            dd = min(d-k, dv)
            do l = 0, dd
                SKL(k,l,:) = 0.0
                do si = 0, q
                    SKL(k,l,:) =  SKL(k,l,:) + NV_ders(l,si)*temp(si,:)
                end do
            end do
        end do
    end subroutine bspline_surface_gradient
    
    subroutine bivariate_rational_gradient(Aders, wders, d, SKL)
        ! Algorithm A4.4 from the NURBS Book
        ! Compute S(u,v) derivatives
        ! from Sw(u,v) derivatives
        ! Input: Aders, dwers, d
        ! Output: SKL
        use nurbs_curve_module, only: binomial

        real, dimension(0:,0:,0:), intent(in) :: Aders
        real, dimension(0:,0:), intent(in) :: wders
        integer, intent(in) :: d
        real, intent(out), dimension(0:,0:,0:) :: SKL

        integer :: k, l, j, i
        real :: bc, bc2, bc3
        real, dimension(0:2) :: v, v2
        
        do k = 0, d
            do l = 0, d-k
                v = Aders(k,l,:)
                do j = 1, l
                    call binomial(l, j, bc)
                    v = v - bc*wders(0,j)*SKL(k,l-j,:)
                end do

                do i = 1, k
                    call binomial(k, i, bc)
                    v = v - bc*wders(i,0)*SKL(k-i,l,:)
                    v2 = 0.0
                    do j = 1, l
                        call binomial(l, j, bc2)
                        v2 = v2 + bc2*wders(i,j)*SKL(k-i,l-j,:)
                    end do
                    call binomial(k, i, bc3)
                    v = v - bc3*v2
                end do
                SKL(k,l,:) = v/wders(0,0)
            end do
        end do
    end subroutine bivariate_rational_gradient
end module nurbs_surface_module