program test_fortran_scripts
    use bspline_basis_functions
    implicit none

    character(len=6) :: real_out_format
    character(len=6) :: integer_out_format

    integer :: p, i, U_arr1_size, n, mid, ii, nd
    real :: u, Nip
    real, dimension(11) :: U_arr1
    real, dimension(3) :: N_arr
    real, dimension(2,3) :: ders

    real_out_format = "(F5.3)"
    integer_out_format = "(I3)"

    U_arr1 = [0,0,0,1,2,3,4,4,5,5,5]
    U_arr1_size = size(U_arr1)

    p = 2
    i = 4
    u = 5./2

    call one_basis_function(p, U_arr1, U_arr1_size, i, u, Nip)
    ! print "(a, I3, a, F5.3)", "N(i=", i, ")= ", Nip

    do i = 1, 5
        call one_basis_function(p, U_arr1, U_arr1_size, i, u, Nip)
        ! print "(a, I3, a, F5.3)", "N(i=", i, ")= ", Nip
    end do

    p = 2
    u = 5./2
    n = 7

    call find_span(n, p, u, U_arr1, mid)
    ! print integer_out_format, mid

    p = 2
    i = 4
    u = 5./2

    call basis_function(i, u, U_arr1_size, p, U_arr1, N_arr)

    nd = 1
    call der_basis_functions(i, u, p, nd, U_arr1_size, U_arr1, ders)

end program test_fortran_scripts