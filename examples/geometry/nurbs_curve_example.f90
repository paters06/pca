program nurbs_curve_example
    use nurbs_curve_module
    use utils
    implicit none

    integer :: p
    real :: Rmax
    real, dimension(5,2) :: P_array
    real, dimension(5,1) :: w_array
    real, dimension(8) :: U_array
    real, dimension(:,:), allocatable :: cpts, dcpts
    character(:), allocatable :: file_name

    Rmax = 1.0
    
    P_array = reshape((/Rmax, Rmax, 0., -Rmax, -Rmax, & 
                0., Rmax, Rmax, Rmax, 0./), shape(P_array))
    w_array = reshape((/1.,0.5*sqrt(2.),1.,0.5*sqrt(2.),1./), shape(w_array))
    U_array = (/0.,0.,0.,0.5,0.5,1.,1.,1./)
    p = 2

    call create_curve(P_array, w_array, U_array, p, cpts)
    call create_tangent_curve(P_array, w_array, U_array, p, dcpts)

    file_name = "cpts.txt"
    ! call export_matrix(cpts, file_name)

    file_name = "dcpts.txt"
    ! call export_matrix(dcpts, file_name)

end program nurbs_curve_example