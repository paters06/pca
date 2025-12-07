module reference_solutions
    implicit none
contains
    subroutine cantilever_1D_point_load(xpts, L, E, I, load, theo_vals)
        use utils, only: print_matrix, export_matrix
        real, intent(in) :: L, E, I, load
        real, dimension(:,:), intent(in) :: xpts
        real, dimension(:,:), intent(out), allocatable :: theo_vals
        real :: y, a
        integer :: num_points, j
        character(:), allocatable :: file_name

        num_points = size(xpts, 1)

        allocate(theo_vals(num_points,2))

        do j = 1, num_points
            a = L - xpts(j, 1)
            y = ((-load)/(6.*E*I))*(2.*(l**3) - 3.*(l**2)*a + (a**3))
            theo_vals(j,1) = xpts(j,1)
            theo_vals(j,2) = y
        end do

        file_name = "Theoretical_values.txt"
        call export_matrix(theo_vals, file_name)
        
    end subroutine cantilever_1D_point_load
end module reference_solutions