module nurbs_curve
    implicit none
contains
    subroutine binomial(a, b, bc)
        integer, intent(in) :: a
        integer, intent(in) :: b
        real, intent(out) :: bc
        
        integer :: j

        bc = 1.0
        do j = 1, b
            bc = bc*((a+1-j)/j)
        end do
    end subroutine binomial
    
    subroutine weighted_control_points(P, w, Pw)
        ! Convert from real control points to homogeneous ones
        real, intent(in) :: P(:,:), w(:,:)
        real, dimension(:,:), intent(out), allocatable :: Pw(:,:)
        integer :: num_control_points, num_cols, j

        num_control_points = size(P,1)
        num_cols = size(P,2) + 1

        allocate(Pw(num_control_points,num_cols))
        Pw(:,1:num_cols-1) = P

        do j = 1,num_control_points
            Pw(j, 1:num_cols-1) = Pw(j, 1:num_cols-1)*w(j,1)
            Pw(j, num_cols) = w(j,1)
        end do
    end subroutine weighted_control_points

    subroutine geometric_control_points(Pw, P, w)
        ! Convert homogeneous control points to geometric ones
        real, intent(in) :: Pw(:,:)
        real, dimension(:,:), intent(out), allocatable :: P, w
        integer :: num_control_points, num_cols, j

        num_control_points = size(Pw,1)
        num_cols = size(Pw,2) - 1

        allocate(P(num_control_points,num_cols))
        allocate(w(num_control_points,1))

        do j = 1,num_control_points
            P(j,:) = Pw(j,1:num_cols)/Pw(j,num_cols+1)
            w(j,1) = Pw(j,num_cols+1)
        end do
    end subroutine geometric_control_points

    subroutine curve_point(n, p, m, U_arr, Pw, u, C)
        ! Compute point on rational B-spline curve
        ! Input: n, p, U, Pw, u
        ! n is the number of control points
        ! m is the number of knots
        ! N has p+1 nonzero elements
        ! Output: C
        use bspline_basis_functions

        integer, intent(in) :: n
        integer, intent(in) :: p
        integer, intent(in) :: m
        real, intent(in) :: u
        real, dimension(m), intent(in) :: U_arr
        real, dimension(n,3), intent(in) :: Pw
        real, dimension(p+1) :: N_arr
        real, dimension(2) :: Cw
        real :: w
        real, dimension(2), intent(out) :: C

        integer :: span, j

        call find_span(m, p, u, U_arr, span)
        call basis_function(span, u, m, p, U_arr, N_arr)

        Cw = 0.0
        w = 0.0

        do j = 0, p
            Cw = Cw + N_arr(j+1)*Pw(span-p+j+1,1:2)
            w = w + N_arr(j+1)*Pw(span-p+j+1,3)
        end do

        C = Cw/w
    end subroutine curve_point

    subroutine create_curve(P_array, w_array, U_array, p, cpts)
        use utils
        integer, intent(in) :: p
        real, intent(in), dimension(:,:) :: P_array
        real, intent(in), dimension(:,:) :: w_array
        real, intent(in), dimension(:) :: U_array
        real, intent(out), dimension(:,:), allocatable :: cpts
        real :: U_min, U_max, u
        integer :: num_points, i, n, m
        real, dimension(:,:), allocatable :: Pw
        real, dimension(:), allocatable :: Ci

        num_points = 41
        n = size(P_array, 1)
        m = size(U_array)
        U_max = maxval(U_array)
        U_min = minval(U_array)

        call weighted_control_points(P_array, w_array, Pw)

        allocate(cpts(num_points+1, 2))
        allocate(Ci(size(P_array,2)))

        do i = 1, num_points+1
            u = ((U_max - U_min)/num_points)*(i-1) + U_min
            call curve_point(n, p, m, U_array, Pw, u, Ci)
            cpts(i,1) = Ci(1)
            cpts(i,2) = Ci(2)
        end do
    end subroutine create_curve

    subroutine compute_curve_derivatives(Pw_array, U_array, p, d, u, Cki)
        ! Algorithm 3.2 from The NURBS Book
        ! Compute curve derivatives
        ! Input: n, p, U_array, Pw_array, u, d
        ! Input:
        ! n -> number of control points is n+1
        ! U_array -> knot vector
        ! Pw_array -> weighted control points
        ! p -> degree of the curve
        ! d -> dth derivative
        ! Output: CKi (d+1, p+1)
        use bspline_basis_functions
        use utils

        real, dimension(:,:), intent(in) :: Pw_array(:,:), U_array(:)
        real, intent(in) :: u
        integer, intent(in) :: p, d
        real, dimension(:,:), intent(out), allocatable :: CKi
        
        integer :: du, span, k, j, m
        real, dimension(:,:), allocatable :: nders

        allocate(Cki(d+1,3))

        du = min(d, p)
        do k = p+1, d
            Cki(k+1,:) = 0.0
        end do

        m = size(U_array,1)

        allocate(nders(d+1,p+1))

        call find_span(m, p, u, U_array, span)
        call der_basis_functions(span, u, p, d, m, U_array, nders)
        ! call print_matrix(nders)
        do k = 0, du
            CKi(k+1,:) = 0.0
            do j = 0, p
                Cki(k+1,:) = Cki(k+1,:) + nders(k+1,j+1)*Pw_array(span-p+j+1,:)
            end do
        end do

        ! call print_matrix(Cki)
    end subroutine compute_curve_derivatives

    subroutine compute_curve_derivatives_2(P_array, U_array, p, d, u, Cki)
        ! Algorithm 3.2 from The NURBS Book
        ! Compute curve derivatives
        ! Input: n, p, U_array, Pw_array, u, d
        ! Input:
        ! n -> number of control points is n+1
        ! U_array -> knot vector
        ! Pw_array -> weighted control points
        ! p -> degree of the curve
        ! d -> dth derivative
        ! Output: CKi (d+1, p+1)
        use bspline_basis_functions
        use utils

        real, dimension(:,:), intent(in) :: P_array(:,:), U_array(:)
        real, intent(in) :: u
        integer, intent(in) :: p, d
        real, dimension(:,:), intent(out), allocatable :: CKi
        
        integer :: du, span, k, j, m
        real, dimension(:,:), allocatable :: nders

        allocate(Cki(d+1,3))

        du = min(d, p)
        do k = p+1, d
            Cki(k+1,:) = 0.0
        end do

        m = size(U_array,1)

        allocate(nders(d+1,p+1))

        call find_span(m, p, u, U_array, span)
        call der_basis_functions(span, u, p, d, m, U_array, nders)
        ! call print_matrix(nders)
        do k = 0, du
            CKi(k+1,:) = 0.0
            do j = 0, p
                Cki(k+1,:) = Cki(k+1,:) + nders(k+1,j+1)*P_array(span-p+j+1,:)
            end do
        end do

        ! call print_matrix(Cki)
    end subroutine compute_curve_derivatives_2

    subroutine rat_curve_derivs(Pw_array, U_array, p, d, u, CK)
        ! Algorithm A4.2 from The NURBS Book
        ! Compute C(u) derivatives from Cw(u) derivatives
        ! Input: p, d
        ! Output: CK
        use bspline_basis_functions
        use utils
        
        ! integer, intent(in) :: n
        integer, intent(in) :: p, d
        real, intent(in) :: u
        real, dimension(:,:), intent(in) :: Pw_array
        real, dimension(:), intent(in) :: U_array
        real, dimension(d+1, 2), intent(out) :: CK
        ! real, dimension(:,:), allocatable :: Pw
        real, dimension(:,:), allocatable :: Cwders, Cders, w_ders
        
        integer :: ndims, n, k, i
        real :: bc

        real, dimension(d+1,2) :: Aders
        real, dimension(d+1,1) :: wders
        real, dimension(2) :: v

        n = size(Pw_array,1)
        ndims = size(Pw_array,2)

        allocate(Cwders(d+1,3))
        allocate(Cders(d+1,2))
        allocate(w_ders(d+1,1))

        call compute_curve_derivatives(Pw_array, U_array, p, d, u, Cwders)

        call geometric_control_points(Cwders, Cders, w_ders)

        CK(1,:) = Cders(1,:)

        Aders = Cwders(:,1:2)
        wders = reshape(Cwders(:,3), shape(wders))

        ! call print_matrix(Cwders)
        ! call print_matrix(Aders)
        ! call print_matrix(wders)

        do k = 0, d
            v = Aders(k+1,:)
            do i = 1, k
                if (k >= i) then
                    call binomial(k, i, bc)
                    v = v - bc*wders(i+1,1)*CK(k+1-i,:)
                end if
            end do
            CK(k+1,:) = v/wders(1,1)
        end do

        ! call print_matrix(CK)
    end subroutine rat_curve_derivs

    subroutine compute_tangent_vector(Pw_array, U_array, p, u, tang_vec)
        use utils
        integer, intent(in) :: p
        real, intent(in) :: u
        real, intent(in), dimension(:,:) :: Pw_array
        real, intent(in), dimension(:) :: U_array
        real, dimension(:,:), allocatable :: Cwders, Aders, wders
        real, dimension(2) :: Ci, Cders
        real, dimension(2), intent(out) :: tang_vec

        integer :: d

        d = 1

        allocate(Cwders(d+1,3))
        allocate(Aders(d+1,2))
        allocate(wders(d+1,1))

        call compute_curve_derivatives(Pw_array, U_array, p, d, u, Cwders)
        Aders = Cwders(:,1:2)
        wders = reshape(Cwders(:,3), shape(wders))

        Ci = Aders(1,:)/wders(1,1)

        Cders = (Aders(d+1,:) - wders(d+1,1)*Ci)/wders(1,1)
        tang_vec = Cders
        ! call print_row_vector(Cders)
    end subroutine

    subroutine compute_tangent_vector_2(P_array, U_array, p, u, tang_vec)
        use utils
        integer, intent(in) :: p
        real, intent(in) :: u
        real, intent(in), dimension(:,:) :: P_array
        real, intent(in), dimension(:) :: U_array
        real, dimension(:,:), allocatable :: Cders
        real, dimension(2), intent(out) :: tang_vec

        integer :: d

        d = 1

        allocate(Cders(d+1,2))

        call compute_curve_derivatives_2(P_array, U_array, p, d, u, Cders)
        tang_vec = Cders(d+1,:)
    end subroutine compute_tangent_vector_2

    subroutine create_tangent_curve(P_array, w_array, U_array, p, dcpts)
        use utils
        integer, intent(in) :: p
        real, intent(in), dimension(:,:) :: P_array
        real, intent(in), dimension(:,:) :: w_array
        real, intent(in), dimension(:) :: U_array
        real, intent(out), dimension(:,:), allocatable :: dcpts
        real, dimension(:,:), allocatable :: Cwders, Aders, wders
        real, dimension(:,:), allocatable :: Pw
        real, dimension(2) :: Ci, Cders

        integer :: num_points, m, n, d, i
        real :: U_min, U_max, u

        num_points = 41
        n = size(P_array, 1)
        m = size(U_array)
        U_max = maxval(U_array)
        U_min = minval(U_array)

        d = 1

        call weighted_control_points(P_array, w_array, Pw)

        allocate(dcpts(num_points+1, 2))

        allocate(Cwders(d+1,2))
        allocate(Aders(d+1,2))
        allocate(wders(d+1,2))

        do i = 1, num_points+1
            u = ((U_max - U_min)/num_points)*(i-1) + U_min
            ! call rat_curve_derivs(Pw, U_array, p, d, u, CKi)
            call compute_curve_derivatives(Pw, U_array, p, d, u, Cwders)
            Aders = Cwders(:,1:2)
            wders = reshape(Cwders(:,3), shape(wders))

            Ci = Aders(1,:)/wders(1,1)

            Cders = (Aders(d+1,:) - wders(d+1,1)*Ci)/wders(1,1)

            dcpts(i,1) = Cders(1)
            dcpts(i,2) = Cders(2)
            ! call print_matrix(Cki)
        end do
    end subroutine create_tangent_curve
end module nurbs_curve