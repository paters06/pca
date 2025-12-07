program refinement_example_4
    use utils
    use curve_refinement
    use nurbs_curve
    implicit none

    character(:), allocatable :: file_name
    character(len=50), dimension(:), allocatable :: line_array

    real, dimension(:,:), allocatable :: ctrl_pts, P_pts, w_pts, Ph_pts, wh_pts, cpts_1, cpts_2
    real, dimension(:), allocatable :: U_knot, Ubar
    integer :: pinp, ph

    integer :: size_1, size_2, size_vec

    character(len=1), dimension(:), allocatable :: ref_list

    file_name = "input_file_curve.txt"
   
    call import_data(file_name, line_array)
    call convert_data(line_array, pinp, U_knot, ctrl_pts, ref_list)

    size_1 = size(ctrl_pts, 1)
    size_2 = size(ctrl_pts, 2) - 1
    size_vec = size(U_knot)

    allocate(P_pts(0:size_1-1,0:size_2-1))
    allocate(w_pts(0:size_1-1,0))

    P_pts = ctrl_pts(1:size_1,1:size_2)
    w_pts = reshape(ctrl_pts(1:size_1,size_2+1), (/size_1, 1/))

    call print_matrix(P_pts)
    call print_matrix(w_pts)
    print *, ref_list

    call spline_refinement(pinp, P_pts, w_pts, U_knot, ref_list, ph, Ph_pts, wh_pts, Ubar)

    call print_row_vector(Ubar)
    call print_matrix(Ph_pts)
    call print_matrix(wh_pts)

    call create_curve(P_pts, w_pts, U_knot, pinp, cpts_1)
    call create_curve(Ph_pts, wh_pts, Ubar, ph, cpts_2)

    call assess_refinement(cpts_1, cpts_2)

end program refinement_example_4