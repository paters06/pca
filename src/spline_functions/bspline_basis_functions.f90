module bspline_basis_functions
    implicit none
contains
    subroutine find_span(n, p, u, U_arr, mid)
        ! Algorithm A2.1 from The NURBS Book
        ! Determine knot span index
        ! Input: n, p, u, U
        ! Return: the knot span index
        !
        implicit none
    
        integer, intent(in) :: n
        integer, intent(in) :: p
        real, intent(in) :: u
        real, dimension(n) :: U_arr
    
        integer, intent(out) :: mid

        !f2py intent(in) n, p, u, U_arr
        !f2py integer intent(hide), depend(U_arr) :: n = len(U_arr)
        !f2py intent(out) mid
    
        integer :: low, high
    
        integer :: ni, mid_i
    
        ni = n + 1
    
        ! Special case
        if (abs(u - U_arr(ni+1)) < 1e-5) then
            mid = n
        end if
    
        ! Starting binary search
        low = p
        high = ni + 1
        mid_i = (low + high)/2
    
        do while (u < U_arr(mid_i) .or. (u > U_arr(mid_i+1) .or. abs(u - U_arr(mid_i+1)) < 1e-5))
            if (u < U_arr(mid_i)) then
                high = mid_i
            else
                low = mid_i
            end if
    
            mid_i = (low + high)/2
        end do
    
        mid = mid_i - 1
    
    end subroutine find_span

    subroutine basis_function(i, u, m, p, U_arr, N_arr)
        ! Algorithm A2.2 from The NURBS Book
        ! Compute the nonvanishing basis functions
        ! Input: i,u,p,U
        ! Output: N
    
        implicit none
    
        integer, intent(in) :: i
        integer, intent(in) :: p
        integer, intent(in) :: m
        real, intent(in) :: u
        real, dimension(m), intent(in) :: U_arr
    
        real, dimension(p+1), intent(out) :: N_arr
        real, dimension(p+1) :: left_arr
        real, dimension(p+1) :: right_arr

        !f2py intent(in) i, u, p, U_arr
        !f2py intent(out) N_arr
        !f2py integer intent(hide), depend(U_arr) :: m = len(U_arr)
    
        integer :: j, r, ii
        real :: saved, temp
    
        N_arr(1) = 1.
    
        ii = i+1
    
        first_loop: do j = 1, p
            left_arr(j) = u - U_arr(ii+1-j)
            right_arr(j) = U_arr(ii+j) - u
            saved = 0.
            second_loop: do r = 0, j-1
                temp = N_arr(r+1)/(right_arr(r+1) + left_arr(j-r))
                N_arr(r+1) = saved + right_arr(r+1)*temp
                saved = left_arr(j-r)*temp
            end do second_loop
            N_arr(j+1) = saved
        end do first_loop
    
        ! do j = 1, p+1
        !     print "(a, I3, a, F5.3)", "N(i=", j+p, ")= ", N_arr(j)
        ! end do
    
    end subroutine basis_function

    subroutine one_basis_function(p, U_arr, U_arr_size, i, u, Nip)
        ! Algorithm A2.4 from The NURBS Book
        ! Compute the basis function Nip
        ! Input: p, m, U, i, u
        ! Output: Nip
        implicit none
    
        integer, intent(in) :: p
        integer, intent(in) :: i
        integer, intent(in) :: U_arr_size
        real, intent(in) :: u
    
        real, dimension(U_arr_size), intent(in) :: U_arr
        real, intent(out) :: Nip
    
        !f2py intent(in) p, i, u, U_arr_size
        !f2py intent(in) U_arr
        !f2py integer intent(hide), depend(U_arr) :: U_arr_size = len(U_arr)
        !f2py intent(out) Nip

        integer :: k, j, m
        integer :: ii
        real saved, Uleft, Uright, temp
        real, dimension(U_arr_size - 1 - p) :: N_arr
    
        ! print "(I3)", size(U_arr) - 1
    
        ii = i+1
    
        m = U_arr_size
        Nip = 0.
        N_arr = 0.0
    
        first_if: if ((ii == 1 .and. u == U_arr(1)) .or. (ii == m - p .and. u == U_arr(m))) then
            Nip = 1.0
        end if first_if
    
        second_if: if (u < U_arr(ii) .or. u >= U_arr(ii+p+1)) then
            Nip = 0.0
        end if second_if
    
        first_loop: do j = 0, p
            if (u >= U_arr(ii+j) .and. u < U_arr(ii+j+1)) then
                N_arr(j+1) = 1.0
            else
                N_arr(j+1) = 0.0
            end if
        end do first_loop
    
        second_loop: do k = 1, p
            if (N_arr(1) == 0.0) then
                saved = 0.0
            else
                saved = ((u - U_arr(ii))*N_arr(1)/(U_arr(ii+k) - U_arr(ii)))
            end if
            
            do j = 0, p-k
                Uleft = U_arr(ii+j+1)
                Uright = U_arr(ii+j+k+1)
                if (N_arr(j+1+1) == 0.0) then
                    N_arr(j+1) = saved
                    saved = 0.0
                else
                    temp = N_arr(j+1+1)/(Uright - Uleft)
                    N_arr(j+1) = saved + (Uright - u)*temp
                    saved = (u - Uleft)*temp
                end if
            end do
        end do second_loop
        
        Nip = N_arr(1)
    end subroutine one_basis_function

    subroutine der_basis_functions(i, u, p, n, m, U_arr, ders)
        ! Algorithm A2.3 from The NURBS Book
        ! Compute nonzero basis functions and their
        ! derivatives. First section is A2.2 modified
        ! to store functions and knot differences
        ! Input: i,u,p,n,U
        ! n is the nth derivative to compute
        ! Output: ders
        implicit none
    
        integer, intent(in) :: i
        integer, intent(in) :: p
        integer, intent(in) :: n
        integer, intent(in) :: m
        real, intent(in) :: u
        real, intent(in), dimension(m) :: U_arr
        real, intent(out), dimension(n+1,p+1) :: ders

        !f2py intent(in) i, u, p, n, m, U_arr
        !f2py intent(out) ders
        !f2py integer intent(hide), depend(U_arr) :: m = len(U_arr)
    
        real, dimension(p+1, p+1) :: ndu
        real, dimension(2, p+1) :: a
        real, dimension(p+1) :: left, right
    
        integer :: j, k, r
        integer :: j1, j2, rk, pk, s1, s2
        real :: d, saved, temp
    
        integer :: ii
    
        ii = i+1
        ndu = 0.0
        left = 0.0
        right = 0.0
        ders = 0.0
        a = 0.0
    
        ndu(1,1) = 1.0
    
        do j = 1, p
            left(j+1) = u - U_arr(ii+1-j)
            right(j+1) = U_arr(ii+j) - u
            saved = 0.0
            do r = 0, j-1
                ! Lower triangle
                ndu(j+1,r+1) = right(r+1+1) + left(j-r+1)
                temp = ndu(r+1,j-1+1)/ndu(j+1,r+1)
                ! Upper triangle
                ndu(r+1,j+1) = saved + right(r+1+1)*temp
                saved = left(j-r+1)*temp
            end do
            ndu(j+1,j+1) = saved
        end do
    
        ! Load the basis functions
        do j = 1, p+1
            ders(1,j) = ndu(j,p+1)
        end do
    
        ! This section computes the derivatives
        ! according to Eq (2.9) from The NURBS Book
        do r = 0, p
            ! Loop over function index
            s1 = 1
            s2 = 2
            ! Alternate rows in array a
            a(1,1) = 1.0
            ! Look to compute kth derivative
            do k = 1, n
                d = 0.0
                rk = r - k
                pk = p - k
    
                if (r >= k) then
                    a(s2,1) = a(s1,1)/ndu(pk+1+1,rk+1)
                    d = a(s2,1)*ndu(rk+1,pk+1)
                end if
    
                if (rk >= -1) then
                    j1 = 1
                else
                    j1 = -rk
                end if
    
                j1 = j1 + 1
    
                if ((r-1) <= pk) then
                    j2 = k - 1
                else
                    j2 = p - r
                end if
    
                j2 = j2 + 1
    
                do j = j1, j2
                    a(s2,j) = (a(s1,j) - a(s1,j-1))/ndu(pk+1,rk+j)
                    d = d + a(s2,j)*ndu(rk+j,pk)
                end do
    
                if (r <= pk) then
                    a(s2,k+1) = -a(s1,k)/ndu(pk+1+1,r+1)
                    d = d + a(s2,k+1)*ndu(r+1,pk+1)
                end if
    
                ders(k+1,r+1) = d
                ! Switch rows
                j = s1
                s1 = s2
                s2 = j
            end do
        end do
        
        ! Multiply through by the correct factors (Eq (2.9))
        r = p
        do k = 1, n
            do j = 0, p
                ders(k+1,j+1) = r*ders(k+1,j+1)
            end do
    
            r = (p-k)*r
        end do
    
        ! do j = 1, n+1
        !     print *, ders(j, 1:p+1)
        !     ! print "(a, I3, a, F5.3)", "dN(i=", j+p, ")= ", ders(j,1:p+1)
        ! end do
    
    end subroutine der_basis_functions

end module bspline_basis_functions