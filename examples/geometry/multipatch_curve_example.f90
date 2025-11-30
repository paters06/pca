program nurbs_curve_example
    use nurbs_curve_module
    use utils
    use derived_types
    use input_output
    implicit none

    character(:), allocatable :: file_name
    character(len=50) :: file_name_inp
    character(len=50), dimension(:), allocatable :: line_array
    character(len=1), dimension(:), allocatable :: ref_list
    type(nurbs_curve) :: input_nurbs_curve
    real, dimension(:,:), allocatable :: cpts

    call get_command_argument(1, file_name_inp)
    file_name = file_name_inp

    call import_data(file_name, line_array)
    call convert_data_to_curve(line_array, input_nurbs_curve, ref_list)

    call create_curve(input_nurbs_curve, cpts)
end program nurbs_curve_example