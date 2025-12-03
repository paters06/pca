module derived_types
    implicit none

    type :: boundary_condition
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
        integer :: p, q
        real, dimension(:), allocatable :: U_knot, V_knot
        real, dimension(:,:), allocatable :: control_points, weight_points
    end type nurbs_surface

    type :: interface_line
        character(len=1) :: dir
        real :: UV_param
    end type interface_line
end module derived_types