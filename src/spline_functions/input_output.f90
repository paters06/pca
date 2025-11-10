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

        ! do i = 0, size(line_array)-1
        !     print format_line, line_array(i)
        ! end do
    end subroutine import_data

    subroutine convert_data_to_curve(line_array, p, U_knot, ctrl_pts, refn_array)
        character(len=*), dimension(0:), intent(in) :: line_array
        integer, intent(out) :: p
        real, dimension(:), allocatable, intent(out) :: U_knot
        real, dimension(:,:), allocatable, intent(out) :: ctrl_pts
        character(len=1), dimension(:), allocatable, intent(out), optional :: refn_array

        character(:), allocatable :: format_line, second_format_line
        ! character(len=50), dimension(:), allocatable :: label_array
        character(:), allocatable :: control_points_label, p_label, uknot_label, end_label, refinement_label
        integer :: i, label_counter
        integer :: control_points_flag, p_flag, uknot_flag, refinement_flag
        integer, dimension(:), allocatable :: label_position

        format_line = "(A)"
        second_format_line = "(A, I3)"

        control_points_label = "*CONTROL_POINTS"
        p_label = "*p_ORDER"
        uknot_label = "*U_KNOT"
        refinement_label = "*REFN"
        end_label = "**END**"

        ! label_array = (/control_points_label, p_label, uknot_label, end_label/)

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

        ! do i = 1, 5
        !     print *, label_position(i)
        ! end do

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

        ! print *, p
        ! call print_row_vector(U_knot)
        ! call print_matrix(ctrl_pts)
    end subroutine

    subroutine convert_data_to_surface(line_array, p, q, U_knot, V_knot, ctrl_pts, refn_array)
        character(len=*), dimension(0:), intent(in) :: line_array
        integer, intent(out) :: p, q
        real, dimension(:), allocatable, intent(out) :: U_knot, V_knot
        real, dimension(:,:), allocatable, intent(out) :: ctrl_pts
        ! character(len=1), dimension(:), allocatable, intent(out), optional :: refn_array
        character(len=1), dimension(:,:), allocatable, intent(out), optional :: refn_array

        character(:), allocatable :: format_line, second_format_line
        ! character(len=50), dimension(:), allocatable :: label_array
        character(:), allocatable :: control_points_label, p_label, q_label, uknot_label, vknot_label, end_label, refinement_label
        integer :: i, label_counter
        integer :: control_points_flag, p_flag, q_flag, uknot_flag, vknot_flag, refinement_flag, end_flag
        integer, dimension(:), allocatable :: label_position

        format_line = "(A)"
        second_format_line = "(A, I3)"

        control_points_label = "*CONTROL_POINTS"
        p_label = "*p_ORDER"
        q_label = "*q_ORDER"
        uknot_label = "*U_KNOT"
        vknot_label = "*V_KNOT"
        refinement_label = "*REFN"
        end_label = "**END**"

        ! label_array = (/control_points_label, p_label, uknot_label, end_label/)

        label_counter = 0
        allocate(label_position(7))
        label_position = 0

        control_points_flag = 0
        p_flag = 0
        q_flag = 0
        uknot_flag = 0
        vknot_flag = 0
        refinement_flag = 0
        end_flag = 0

        do i = 0, size(line_array)-1
            if (control_points_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                control_points_flag = 1
            else if (p_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                p_flag = 1
            else if (q_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                q_flag = 1
            else if (uknot_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                uknot_flag = 1
            else if (vknot_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                vknot_flag = 1
            else if (refinement_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                refinement_flag = 1
            else if (end_label == line_array(i)) then
                label_counter =  label_counter + 1
                label_position(label_counter) = i
                end_flag = 1
            end if
        end do

        ! do i = 1, size(label_position)
        !     print *, label_position(i)
        ! end do

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
            call read_string_matrix(line_array(label_position(6)+1:label_position(7)-1), ",", refn_array)
        end if

        ! print *, p
        ! call print_row_vector(U_knot)
        ! call print_matrix(ctrl_pts)
    end subroutine

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

        ! call print_row_vector(array)
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
            ! print "(A)", strings(i)
            call read_real_row_array(strings(i), delim, array)
            temp_matrix(i,1:size(array)) = array
            ! call print_row_vector(array)
        end do

        matrix = temp_matrix(:,1:size(array))
        
        ! call print_matrix(matrix)
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
            ! print "(A)", strings(i)
            call read_row_string(strings(i), delim, array)
            temp_matrix(i,1:size(array)) = array
            ! call print_row_vector(array)
        end do

        matrix = temp_matrix(:,1:size(array))

        ! do i = 1, size(matrix,1)
        !     print *, matrix(i,:)
        ! end do
        
        ! call print_matrix(matrix)
    end subroutine read_string_matrix

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
                ! print '(*(3X, f15.4))', (mat(i,j), j=1,num_cols)
            end do
            close(io)
        else
            open(newunit=io, file=file_name, status="new", action="write")
            do i = 1, num_rows
                ! write(io,"(F6.3, A, F8.3)") mat(i,1), ", ", mat(i,2)
                write(io,"(*(3X, F15.4))") (mat(i,j), j=1,num_cols)
            end do
            close(io)
        end if
    end subroutine export_matrix
end module input_output