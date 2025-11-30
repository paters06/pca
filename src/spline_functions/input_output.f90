module input_output
    use utils
    implicit none
contains
    subroutine import_data(file_name, line_array)
        ! Example taken and adapted from
        ! Figure 5-10
        ! Fortran 95/2003 for Scientists and Engineers
        ! Chapman, 2008
        character(:), intent(in), allocatable :: file_name
        character(len=50) :: line
        character(:), allocatable :: format_line
        integer :: MAX_NUM_VALS = 9999
        character(len=50), dimension(:), allocatable :: line_array_temp
        character(len=*), dimension(:), allocatable, intent(out) :: line_array
        integer :: io, ierror, num_values
        logical :: exists

        num_values = 0

        allocate(line_array_temp(0:MAX_NUM_VALS-1))
        format_line = "(A)"

        inquire(file=file_name, exist=exists)

        if (exists) then
            open(newunit=io, file=file_name, status="old", action="read", iostat=ierror)
            
            if (ierror == 0) then
                do
                    read (io,fmt=format_line,iostat=ierror) line
                    line_array_temp(num_values) = line
                    if (ierror /= 0) exit
                    num_values = num_values + 1
                end do

                if (ierror > 0) then
                    print "(A, I6)", "An error ocurred reading line ", num_values + 1
                else
                    print "(A, I6, A)", "End of file reached. There were", num_values + 1, " values in the file"
                end if
            else
                print "(A, I6)", "Error opening file: IOSTAT ", ierror
            end if

            close(io)
        end if

        allocate(line_array(0:num_values-1))

        line_array = line_array_temp(0:num_values-1)
    end subroutine import_data

    subroutine convert_data_to_curve(line_array, input_curve, refn_array)
        use derived_types
        character(len=*), dimension(0:), intent(in) :: line_array
        integer :: p
        real, dimension(:), allocatable :: U_knot
        real, dimension(:,:), allocatable :: ctrl_pts
        type(nurbs_curve), intent(out) :: input_curve
        character(len=1), dimension(:), allocatable, intent(out), optional :: refn_array

        character(:), allocatable :: format_line, second_format_line
        ! character(len=50), dimension(:), allocatable :: label_array
        character(:), allocatable :: control_points_label, p_label, uknot_label, end_label, refinement_label
        character(:), allocatable :: patch_label
        integer :: i, label_counter
        integer :: control_points_flag, p_flag, uknot_flag, refinement_flag, patch_flag
        integer, dimension(:), allocatable :: label_position

        integer :: size_1, size_2
        real, dimension(:,:), allocatable :: P_pts, w_pts

        format_line = "(A)"
        second_format_line = "(A, I3)"

        control_points_label = "*CONTROL_POINTS"
        p_label = "*p_ORDER"
        uknot_label = "*U_KNOT"
        refinement_label = "*REFN"
        end_label = "**END_GEOMETRY**"

        label_counter = 0
        allocate(label_position(5))
        label_position = 0

        control_points_flag = 0
        p_flag = 0
        uknot_flag = 0
        refinement_flag = 0

        do i = 0, size(line_array)-1
            if (control_points_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                control_points_flag = 1
            else if (p_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                p_flag = 1
            else if (uknot_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                uknot_flag = 1
            else if (refinement_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                refinement_flag = 1
            end if
        end do

        p = 0
        
        ! Reading control points
        if (control_points_flag == 1) then
            call read_real_matrix(line_array(label_position(1)+1:label_position(2)-1), ",", ctrl_pts)
        end if
        
        ! Reading spline degree
        if (p_flag == 1) then
            read (line_array(label_position(2)+1),*) p
        end if

        ! Reading knot vector U
        if (uknot_flag == 1) then
            call read_real_row_array(line_array(label_position(3)+1), ",", U_knot)
        end if

        ! Reading refining types
        if (refinement_flag == 1) then
            call read_row_string(line_array(label_position(4)+1), ",", refn_array)
        end if

        size_1 = size(ctrl_pts, 1)
        size_2 = size(ctrl_pts, 2) - 1

        allocate(P_pts(0:size_1-1,0:size_2-1))
        allocate(w_pts(0:size_1-1,0))

        P_pts = ctrl_pts(1:size_1,1:size_2)
        w_pts = reshape(ctrl_pts(1:size_1,size_2+1), (/size_1, 1/))

        input_curve%p = p
        input_curve%U_knot = U_knot
        input_curve%control_points = P_pts
        input_curve%weight_points = w_pts
    end subroutine

    subroutine convert_data_to_surface(line_array, input_surf, refn_array)
        use derived_types
        character(len=*), dimension(0:), intent(in) :: line_array
        integer :: p, q
        real, dimension(:), allocatable :: U_knot, V_knot
        real, dimension(:,:), allocatable :: ctrl_pts
        type(nurbs_surface), intent(out) :: input_surf
        ! character(len=1), dimension(:), allocatable, intent(out), optional :: refn_array
        character(len=1), dimension(:,:), allocatable, intent(out), optional :: refn_array

        character(:), allocatable :: format_line, second_format_line
        ! character(len=50), dimension(:), allocatable :: label_array
        character(:), allocatable :: start_geom_label, control_points_label
        character(:), allocatable :: p_label, q_label, uknot_label, vknot_label
        character(:), allocatable :: end_geom_label, refinement_label
        
        integer :: i, label_counter, i_start, i_end
        integer :: start_geom_flag, control_points_flag, p_flag, q_flag
        integer :: uknot_flag, vknot_flag, refinement_flag, end_geom_flag
        integer, dimension(:), allocatable :: label_position

        integer :: size_1, size_2, size_vec_1, size_vec_2
        real, dimension(:,:), allocatable :: P_pts, w_pts

        format_line = "(A)"
        second_format_line = "(A, I3)"

        start_geom_label = "**START_GEOMETRY**"
        control_points_label = "*CONTROL_POINTS"
        p_label = "*p_ORDER"
        q_label = "*q_ORDER"
        uknot_label = "*U_KNOT"
        vknot_label = "*V_KNOT"
        refinement_label = "*REFN"
        end_geom_label = "**END_GEOMETRY**"

        label_counter = 0
        allocate(label_position(6))
        label_position = 0

        start_geom_flag = 0
        control_points_flag = 0
        p_flag = 0
        q_flag = 0
        uknot_flag = 0
        vknot_flag = 0
        refinement_flag = 0
        end_geom_flag = 0

        do i = 0, size(line_array)-1
            if (start_geom_label == line_array(i)) then
                i_start = i
                start_geom_flag = 1
            else if (end_geom_label == line_array(i)) then
                i_end = i
                end_geom_flag = 1
            end if
        end do

        do i = i_start, i_end
            if (control_points_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                control_points_flag = 1
            else if (p_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                p_flag = 1
            else if (q_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                q_flag = 1
            else if (uknot_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                uknot_flag = 1
            else if (vknot_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                vknot_flag = 1
            else if (refinement_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                refinement_flag = 1
            end if
        end do

        p = 0
        q = 0
        
        ! Reading control points
        if (control_points_flag == 1) then
            call read_real_matrix(line_array(label_position(1)+1:label_position(2)-1), ",", ctrl_pts)
        end if
        
        ! Reading spline degree
        if (p_flag == 1) then
            read (line_array(label_position(2)+1),*) p
        end if

        ! Reading spline degree
        if (q_flag == 1) then
            read (line_array(label_position(3)+1),*) q
        end if

        ! Reading knot vector U
        if (uknot_flag == 1) then
            call read_real_row_array(line_array(label_position(4)+1), ",", U_knot)
        end if

        ! Reading knot vector V
        if (vknot_flag == 1) then
            call read_real_row_array(line_array(label_position(5)+1), ",", V_knot)
        end if

        ! Reading refining types
        if (refinement_flag == 1) then
            ! call read_row_string(line_array(label_position(6)+1), ",", refn_array)
            call read_string_matrix(line_array(label_position(6)+1:i_end-1), ",", refn_array)
        end if

        size_1 = size(ctrl_pts, 1)
        size_2 = size(ctrl_pts, 2) - 1
        size_vec_1 = size(U_knot)
        size_vec_2 = size(V_knot)

        allocate(P_pts(0:size_1-1,0:size_2-1))
        allocate(w_pts(0:size_1-1,0))

        P_pts = ctrl_pts(1:size_1,1:size_2)
        w_pts = reshape(ctrl_pts(1:size_1,size_2+1), (/size_1, 1/))

        input_surf%p = p
        input_surf%q = q
        input_surf%U_knot = U_knot
        input_surf%V_knot = V_knot
        input_surf%control_points = P_pts
        input_surf%weight_points = w_pts
    end subroutine

    subroutine convert_data_to_solver(line_array, num_gauss_pts, kappa, bc_array)
        ! bc_array: array with information for the enforcement of boundary conditions
        use derived_types
        character(len=*), dimension(0:), intent(in) :: line_array
        integer, intent(out) :: num_gauss_pts
        real, intent(out) :: kappa

        character(:), allocatable :: format_line, second_format_line
        character(:), allocatable :: start_solver_label, intg_label
        character(:), allocatable :: diffusion_label, bc_label
        character(:), allocatable :: end_solver_label
        integer :: i, label_counter, i_start, i_end
        integer :: start_solver_flag, intg_flag, diffusion_flag, bc_flag
        integer :: end_solver_flag
        integer, dimension(:), allocatable :: label_position
        
        type(boundary_condition), dimension(:), allocatable, intent(out) :: bc_array

        format_line = "(A)"
        second_format_line = "(A, I3)"

        start_solver_label = "**START_SOLVER**"
        intg_label = "*GAUSS_NUM"
        diffusion_label = "*DIFFUSION"
        bc_label = "*ESSENTIAL_COND_DOFS"
        end_solver_label = "**END_SOLVER**"

        label_counter = 0
        allocate(label_position(4))
        label_position = 0

        start_solver_flag = 0
        intg_flag = 0
        diffusion_flag = 0
        bc_flag = 0
        end_solver_flag = 0

        do i = 0, size(line_array)-1
            if (start_solver_label == line_array(i)) then
                i_start = i
            else if (end_solver_label == line_array(i)) then
                i_end = i
            end if
        end do

        do i = i_start, i_end
            if (intg_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                intg_flag = 1
            else if (diffusion_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                diffusion_flag = 1
            else if (bc_label == line_array(i)) then
                label_counter = label_counter + 1
                label_position(label_counter) = i
                bc_flag = 1
            end if
        end do

        ! Reading number of gauss points for numerical integration
        if (intg_flag == 1) then
            read (line_array(label_position(1)+1),*) num_gauss_pts
        end if
        
        ! Reading diffusion conductivity value
        if (diffusion_flag == 1) then
            read (line_array(label_position(2)+1),*) kappa
        end if

        ! Reading essential boundary conditions
        if (bc_flag == 1) then
            call read_derived_type_matrix(line_array(label_position(3)+1:i_end-1), bc_array)
        end if
    end subroutine convert_data_to_solver

    subroutine read_integer_row_array(string, delim, array)
        ! Adapted from Beliavksy
        ! https://fortran-lang.discourse.group/t/allocation-of-array-for-reading-strings-from-a-text-file/5986/2
        
        character(len=*), intent(in) :: string
        character(len=*), intent(in) :: delim
        integer, dimension(:), allocatable, intent(out) :: array
        character(:), allocatable :: substring, remainder_text
        integer, parameter :: MAX_SIZE = 50
        integer, dimension(:), allocatable :: temp
        
        integer :: ipos, string_length, num_values

        num_values = 0
        allocate(temp(MAX_SIZE))
        temp = 0.
        
        remainder_text = string
        ! string_length = len(remainder_text)
        do 
            ipos = index(remainder_text, delim)
            if (ipos /= 0) then
                substring = remainder_text(1:ipos-1)
                ! print "('Number: ', a)", substring
                num_values =  num_values + 1
                read (substring,*) temp(num_values)
                remainder_text = remainder_text(ipos+1:)
                string_length = len(remainder_text)
                ! print "('Remainder length: ', I3)", string_length 
            else
                substring = remainder_text
                num_values =  num_values + 1
                read (substring,*) temp(num_values)
                ! print "('Last number: ', a)", substring
                exit
            end if
        end do

        allocate(array(num_values))
        array = temp(1:num_values)
    end subroutine read_integer_row_array

    subroutine read_real_row_array(string, delim, array)
        ! Adapted from Beliavksy
        ! https://fortran-lang.discourse.group/t/allocation-of-array-for-reading-strings-from-a-text-file/5986/2
        
        character(len=*), intent(in) :: string
        character(len=*), intent(in) :: delim
        real, dimension(:), allocatable, intent(out) :: array
        character(:), allocatable :: substring, remainder_text
        integer, parameter :: MAX_SIZE = 50
        real, dimension(:), allocatable :: temp
        
        integer :: ipos, string_length, num_values

        num_values = 0
        allocate(temp(MAX_SIZE))
        temp = 0.
        
        remainder_text = string
        ! string_length = len(remainder_text)
        do 
            ipos = index(remainder_text, delim)
            if (ipos /= 0) then
                substring = remainder_text(1:ipos-1)
                ! print "('Number: ', a)", substring
                num_values =  num_values + 1
                read (substring,*) temp(num_values)
                remainder_text = remainder_text(ipos+1:)
                string_length = len(remainder_text)
                ! print "('Remainder length: ', I3)", string_length 
            else
                substring = remainder_text
                num_values =  num_values + 1
                read (substring,*) temp(num_values)
                ! print "('Last number: ', a)", substring
                exit
            end if
        end do

        allocate(array(num_values))
        array = temp(1:num_values)
    end subroutine read_real_row_array

    subroutine read_real_matrix(strings, delim, matrix)
        character(len=*), dimension(:), intent(in) :: strings
        character(len=*), intent(in) :: delim
        real, dimension(:,:), allocatable, intent(out) :: matrix
        
        real, dimension(:), allocatable :: array
        real, dimension(:,:), allocatable :: temp_matrix

        integer :: i, num_rows

        integer, parameter :: MAX_SIZE = 50
        
        num_rows = size(strings)
        
        allocate(temp_matrix(num_rows, MAX_SIZE))
        temp_matrix = 0.

        do i = 1, size(strings)
            call read_real_row_array(strings(i), delim, array)
            temp_matrix(i,1:size(array)) = array
        end do

        matrix = temp_matrix(:,1:size(array))
    end subroutine read_real_matrix

    subroutine read_row_string(string, delim, array)
        ! Adapted from Beliavksy
        ! https://fortran-lang.discourse.group/t/allocation-of-array-for-reading-strings-from-a-text-file/5986/2
        
        character(len=*), intent(in) :: string
        character(len=*), intent(in) :: delim
        character(len=1), dimension(:), allocatable, intent(out) :: array
        character(:), allocatable :: substring, remainder_text
        integer, parameter :: MAX_SIZE = 50
        character(len=1), dimension(:), allocatable :: temp
        
        integer :: ipos, string_length, num_values

        num_values = 0
        allocate(temp(MAX_SIZE))
        
        remainder_text = string
        ! string_length = len(remainder_text)
        do 
            ipos = index(remainder_text, delim)
            if (ipos /= 0) then
                substring = remainder_text(1:ipos-1)
                ! print "('Number: ', a)", substring
                num_values =  num_values + 1
                read (substring,*) temp(num_values)
                remainder_text = remainder_text(ipos+1:)
                string_length = len(remainder_text)
                ! print "('Remainder length: ', I3)", string_length 
            else
                substring = remainder_text
                num_values =  num_values + 1
                read (substring,*) temp(num_values)
                ! print "('Last number: ', a)", substring
                exit
            end if
        end do

        allocate(array(num_values))
        array = temp(1:num_values)
    end subroutine read_row_string

    subroutine read_string_matrix(strings, delim, matrix)
        character(len=*), dimension(:), intent(in) :: strings
        character(len=*), intent(in) :: delim
        character(len=*), dimension(:,:), allocatable, intent(out) :: matrix

        character(len=1), dimension(:), allocatable :: array
        character(len=1), dimension(:,:), allocatable :: temp_matrix

        integer :: i, num_rows

        integer, parameter :: MAX_SIZE = 50
        
        num_rows = size(strings)
        
        allocate(temp_matrix(num_rows, MAX_SIZE))

        do i = 1, size(strings)
            call read_row_string(strings(i), delim, array)
            temp_matrix(i,1:size(array)) = array
        end do

        matrix = temp_matrix(:,1:size(array))
    end subroutine read_string_matrix

    subroutine read_derived_type_matrix(strings, derived_array)
        use derived_types
        character(len=*), dimension(:), intent(in) :: strings
        type(boundary_condition), dimension(:), allocatable, intent(out) :: derived_array
        integer, parameter :: MAX_SIZE = 50
        type(boundary_condition) :: temp

        integer :: num_values, i

        num_values = size(strings)
        allocate(derived_array(num_values))

        do i = 1, num_values
            read (strings(i),*) temp
            derived_array(i) = temp
            ! print "('Direction: ', A)", temp%dir
            ! print "('Parameter: ', F4.2)", temp%UV_param
            ! print "('Prescribed value: ', F6.2)", temp%prescribed_val
        end do
    end subroutine read_derived_type_matrix

    subroutine export_matrix(mat, file_name)
        ! This function is to be generalized for more than two columns
        real, intent(in) :: mat(:,:)
        integer :: num_rows, num_cols, i, j, io, num_characters
        logical :: exists
        character(:), intent(in), allocatable :: file_name

        num_rows = size(mat,1)
        num_cols = size(mat,2)
        num_characters = len(file_name)

        inquire(file=file_name, exist=exists)
        if (exists) then
            open(newunit=io, file=file_name, status="old", action="write")
            do i = 1, num_rows
                write(io,"(*(3X, F15.4))") (mat(i,j), j=1,num_cols)
            end do
            close(io)
        else
            open(newunit=io, file=file_name, status="new", action="write")
            do i = 1, num_rows
                write(io,"(*(3X, F15.4))") (mat(i,j), j=1,num_cols)
            end do
            close(io)
        end if
    end subroutine export_matrix
end module input_output