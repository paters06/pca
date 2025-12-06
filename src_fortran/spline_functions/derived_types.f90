module derived_types
    implicit none

    type :: boundary_condition
        ! id_patch: number of the patch where the boundary condition is applied (NOT INCLUDED DUE TO GFORTRAN BUG)
        ! dir: parameter where the boundary condition will be enforced
        ! UV_param: parameter value
        ! prescribed_val: value to be enforced
        character :: dir
        real :: UV_param
        real :: prescribed_val
    end type boundary_condition

    type :: nurbs_curve
        integer :: p
        real, dimension(:), allocatable :: U_knot
        real, dimension(:,:), allocatable :: control_points, weight_points
    end type nurbs_curve

    type :: nurbs_surface
        ! refn_input has the information of the type of refinement
        ! and the parameter to be applied with
        integer :: p, q
        real, dimension(:), allocatable :: U_knot, V_knot
        real, dimension(:,:), allocatable :: control_points, weight_points
        character(len=1), dimension(:,:), allocatable :: refn_input
        integer :: refn_flag
    end type nurbs_surface

    type :: interface_line
        character(len=1) :: dir
        real :: UV_param
    end type interface_line
end module derived_types