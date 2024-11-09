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

!f2py intent(in) p, i, u, U_arr_size
!f2py intent(in) U_arr
!f2py integer intent(hide), depend(U_arr) :: U_arr_size = len(U_arr)
!f2py intent(out) Nip

    real, dimension(U_arr_size), intent(in) :: U_arr
    real, intent(out) :: Nip

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