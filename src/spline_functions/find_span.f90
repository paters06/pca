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
