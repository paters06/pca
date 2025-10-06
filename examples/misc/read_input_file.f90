program input_file
    use utils
    implicit none

    character(:), allocatable :: file_name
    character(len=50), dimension(:), allocatable :: line_array
    integer :: p
    real, dimension(:), allocatable :: U_knot
    real, dimension(:,:), allocatable :: ctrl_pts

    file_name = "input_file_curve.txt"
   
    call import_data(file_name, line_array)
    call convert_data(line_array, p, U_knot, ctrl_pts)

    print *, p
    call print_row_vector(U_knot)
    call print_matrix(ctrl_pts)

end program input_file