program misc
    implicit none
    INTEGER, DIMENSION(3, 3) :: array
    
    integer :: i, j

    array = reshape((/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /), shape(array))

    do i = 1, 3
        print *, array(i, 1:)
    end do

    print "(a, I3)", "a(2,3)=", array(2,3)
    print "(a, I3)", "a(3,2)=", array(3,2)

end program misc
