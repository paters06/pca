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