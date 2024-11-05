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

    do j = 1, p+1
        print "(a, I3, a, F5.3)", "N(i=", j+p, ")= ", N_arr(j)
    end do

end subroutine basis_function