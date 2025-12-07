program refinement_example_9
    use utils
    use nurbs_surface_module
    use surface_refinement
    use input_output
    use derived_types
    implicit none

    character(:), allocatable :: file_name
    character(len=50), dimension(:), allocatable :: line_array
    character(len=1), dimension(:,:), allocatable :: ref_list
    integer :: refn_flag, interf_flag
    ! type(nurbs_surface) :: input_nurbs_surface, refined_surface
    type(nurbs_surface), dimension(:), allocatable :: input_patches
    type(interface_line), dimension(:), allocatable :: interf_var

    file_name = "input_file_surface.txt"

    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, input_patches, refn_flag, ref_list, interf_flag, interf_var)

    call print_nurbs_surface_info(input_patches(1))
    call print_string_matrix(input_patches(1)%refn_input)
    ! call print_string_matrix(ref_list)

    call surface_spline_refinement(input_patches(1))
    ! call print_nurbs_surface_info(refined_surface)
    ! call assess_surface_refinement(input_nurbs_surface, refined_surface)
end program refinement_example_9