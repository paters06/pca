module bspline_basis_functions
    implicit none
contains
    subroutine find_span(m, p, u, U_arr, mid)
        ! Algorithm A2.1 from The NURBS Book
        ! Determine knot span index
        ! Input: m, p, u, U
        ! m+1 is the number of elements of the knot vector
        ! Return: the knot span index
        implicit none
    
        integer, intent(in) :: m
        integer, intent(in) :: p
        real, intent(in) :: u
        real, dimension(0:m) :: U_arr
        integer, intent(out) :: mid

        integer :: low, high
        integer :: n
    
        n = m - p - 1

        ! Special case
        if (abs(u - U_arr(n+1)) < 1e-5) then
            ! print *, "Special case"
            mid = n
        else
            ! Starting binary search
            low = p
            high = n + 1
            mid = (low + high)/2

            do while (u < U_arr(mid) .or. (u > U_arr(mid+1) .or. abs(u - U_arr(mid+1)) < 1e-5))
                if (u < U_arr(mid)) then
                    high = mid
                else
                    low = mid
                end if
        
                mid = (low + high)/2
            end do
        
        end if
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
        real, dimension(0:m), intent(in) :: U_arr
    
        real, dimension(0:p), intent(out) :: N_arr
        real, dimension(0:p) :: left_arr
        real, dimension(0:p) :: right_arr

        integer :: j, r
        real :: saved, temp
    
        N_arr(0) = 1.
    
        first_loop: do j = 1, p
            left_arr(j) = u - U_arr(i+1-j)
            right_arr(j) = U_arr(i+j) - u
            saved = 0.
            second_loop: do r = 0, j-1
                temp = N_arr(r)/(right_arr(r+1) + left_arr(j-r))
                N_arr(r) = saved + right_arr(r+1)*temp
                saved = left_arr(j-r)*temp
            end do second_loop
            N_arr(j) = saved
        end do first_loop
    
        ! do j = 0, p
        !     print "(a, I3, a, F5.3)", "N(i=", j+p, ")= ", N_arr(j)
        ! end do
    end subroutine basis_function

    subroutine one_basis_function(p, U_arr, m, i, u, Nip)
        ! Algorithm A2.4 from The NURBS Book
        ! Compute the basis function Nip
        ! Input: p, m, U, i, u
        ! Output: Nip
        implicit none
    
        integer, intent(in) :: p
        integer, intent(in) :: i
        integer, intent(in) :: m
        real, intent(in) :: u
    
        real, dimension(0:m), intent(in) :: U_arr
        real, intent(out) :: Nip
    
        integer :: k, j
        real saved, Uleft, Uright, temp
        real, dimension(0:m - 1 - p) :: N_arr
    
        Nip = 0.
        N_arr = 0.0
    
        first_if: if ((i == 0 .and. u == U_arr(0)) .or. (i == m - p - 1 .and. u == U_arr(m))) then
            Nip = 1.0
        end if first_if
    
        second_if: if (u < U_arr(i) .or. u >= U_arr(i+p+1)) then
            Nip = 0.0
        end if second_if
    
        ! Initialize zeroth-degree functions
        first_loop: do j = 0, p
            if (u >= U_arr(i+j) .and. u < U_arr(i+j+1)) then
                N_arr(j) = 1.0
            else
                N_arr(j) = 0.0
            end if
        end do first_loop
    
        ! Compute triangular table
        second_loop: do k = 1, p
            if (N_arr(0) == 0.0) then
                saved = 0.0
            else
                saved = ((u - U_arr(i))*N_arr(0)/(U_arr(i+k) - U_arr(i)))
            end if
            
            do j = 0, p-k
                Uleft = U_arr(i+j+1)
                Uright = U_arr(i+j+k+1)
                if (N_arr(j+1) == 0.0) then
                    N_arr(j) = saved
                    saved = 0.0
                else
                    temp = N_arr(j+1)/(Uright - Uleft)
                    N_arr(j) = saved + (Uright - u)*temp
                    saved = (u - Uleft)*temp
                end if
            end do
        end do second_loop
        
        Nip = N_arr(0)
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
        real, intent(in), dimension(0:m) :: U_arr
        real, intent(out), dimension(0:n,0:p) :: ders

        real, dimension(0:p, 0:p) :: ndu
        real, dimension(0:1, 0:p) :: a
        real, dimension(0:p) :: left, right
    
        integer :: j, k, r
        integer :: j1, j2, rk, pk, s1, s2
        real :: d, saved, temp
    
        ndu = 0.0
        left = 0.0
        right = 0.0
        ders = 0.0
        a = 0.0

        ndu(0,0) = 1.0

        do j = 1, p
            left(j) = u - U_arr(i+1-j)
            right(j) = U_arr(i+j) - u
            saved = 0.0
            do r = 0, j-1
                ! Lower triangle
                ndu(j,r) = right(r+1) + left(j-r)
                temp = ndu(r,j-1)/ndu(j,r)
                ! Upper triangle
                ndu(r,j) = saved + right(r+1)*temp
                saved = left(j-r)*temp
            end do
            ndu(j,j) = saved
        end do

        ! Load the basis functions
        do j = 0, p
            ders(0,j) = ndu(j,p)
        end do

        ! This section computes the derivatives
        ! according to Eq (2.9) from The NURBS Book
        do r = 0, p
            ! Loop over function index
            s1 = 0
            s2 = 1
            ! Alternate rows in array a
            a(0,0) = 1.0
            ! Look to compute kth derivative
            do k = 1, n
                d = 0.0
                rk = r - k
                pk = p - k
    
                if (r >= k) then
                    a(s2,0) = a(s1,0)/ndu(pk+1,rk)
                    d = a(s2,0)*ndu(rk,pk)
                end if
    
                if (rk >= -1) then
                    j1 = 1
                else
                    j1 = -rk
                end if
    
                if ((r-1) <= pk) then
                    j2 = k - 1
                else
                    j2 = p - r
                end if
    
                do j = j1, j2
                    a(s2,j) = (a(s1,j) - a(s1,j-1))/ndu(pk+1,rk+j)
                    d = d + a(s2,j)*ndu(rk+j,pk)
                end do
    
                if (r <= pk) then
                    a(s2,k) = -a(s1,k-1)/ndu(pk+1,r)
                    d = d + a(s2,k)*ndu(r,pk)
                end if
    
                ders(k,r) = d
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
                ders(k,j) = r*ders(k,j)
            end do
    
            r = (p-k)*r
        end do

        ! do j = 0, n
        !     print *, ders(j, 0:p)
        !     ! print "(a, I3, a, F5.3)", "dN(i=", j+p, ")= ", ders(0,0:p)
        ! end do
    end subroutine der_basis_functions
end module bspline_basis_functions