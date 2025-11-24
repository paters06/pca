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
    type(nurbs_surface) :: input_nurbs_surface, refined_surface

    file_name = "input_file_surface.txt"

    call import_data(file_name, line_array)
    call convert_data_to_surface(line_array, input_nurbs_surface, ref_list)

    call print_nurbs_surface_info(input_nurbs_surface)
    call print_string_matrix(ref_list)

    call surface_spline_refinement(input_nurbs_surface, ref_list, refined_surface)
    call print_nurbs_surface_info(refined_surface)
    call assess_surface_refinement(input_nurbs_surface, refined_surface)
end program refinement_example_9