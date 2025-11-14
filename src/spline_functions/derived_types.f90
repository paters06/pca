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
end module derived_types