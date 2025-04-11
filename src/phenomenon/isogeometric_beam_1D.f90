module isogeometric_beam_1D
    implicit none
contains
    subroutine gauss_legendre_nodes_weights(num_gauss_pts, gauss_nodes, gauss_weights)
        ! Info got from https://pomax.github.io/bezierinfo/legendre-gauss.html
        integer, intent(in) :: num_gauss_pts
        real, dimension(:), intent(out), allocatable :: gauss_nodes, gauss_weights

        allocate(gauss_nodes(num_gauss_pts))
        allocate(gauss_weights(num_gauss_pts))

        if (num_gauss_pts == 2) then
            gauss_nodes = (/-0.5774, 0.5774/)
            gauss_weights = (/1.0, 1.0/)
        else if (num_gauss_pts == 3) then
            gauss_nodes = (/-0.7746, 0.0, 0.7746/)
            gauss_weights = (/0.5556, 0.8889, 0.5556/)
        else if (num_gauss_pts == 4) then
            gauss_nodes = (/-0.8611, -0.3400, 0.3400, 0.8611/)
            gauss_weights = (/0.3479, 0.6521, 0.6521, 0.3479/)
        else if (num_gauss_pts == 5) then
            gauss_nodes = (/-0.9062, -0.5385, 0.0, 0.5385, 0.9062/)
            gauss_weights = (/0.2369, 0.4786, 0.5689, 0.4786, 0.2369/)
        else
            print *, "Value error"
        end if
    end subroutine gauss_legendre_nodes_weights
    
    subroutine calculate_integration_points(gauss_pts, ua, ub, intg_pts)
        real, intent(in) :: ua, ub
        real, dimension(:), intent(in) :: gauss_pts
        real, dimension(:), intent(out) :: intg_pts

        integer :: i, num_gauss_pts

        num_gauss_pts = size(gauss_pts)

        do i = 1, num_gauss_pts
            intg_pts(i) = 0.5*(ub - ua)*gauss_pts(i) + 0.5*(ub + ua)
        end do

    end subroutine calculate_integration_points
    
    subroutine element_matrix_stiffness()
    end subroutine element_matrix_stiffness

    subroutine element_force_vector()
    end subroutine element_force_vector

    subroutine assemble_weak_form(U_reduced, num_gauss_pts)
        real, dimension(:), intent(in) :: U_reduced
        integer, intent(in) :: num_gauss_pts

        integer :: num_unique_knots, i
        real, dimension(:), allocatable :: intg_pts, gauss_nodes, gauss_weights

        num_unique_knots = size(U_reduced)

        call gauss_legendre_nodes_weights(num_gauss_pts, gauss_nodes, gauss_weights)

        allocate(intg_pts(num_gauss_pts))

        do i = 1, num_unique_knots - 1
            call calculate_integration_points(gauss_nodes, U_reduced(i), U_reduced(i+1), intg_pts)
            print *, intg_pts
        end do

        ! call element_matrix_stiffness(Ke)
        ! call element_force_vector(Fe)
    end subroutine assemble_weak_form
end module isogeometric_beam_1D