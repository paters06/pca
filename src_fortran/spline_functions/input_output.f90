module input_output
    use utils
    use derived_types
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
        integer :: io, ierror, id_values
        logical :: exists

        id_values = 0

        allocate(line_array_temp(MAX_NUM_VALS))
        format_line = "(A)"

        inquire(file=file_name, exist=exists)

        if (exists) then
            open(newunit=io, file=file_name, status="old", action="read", iostat=ierror)
            
            if (ierror == 0) then
                do
                    read (io,fmt=format_line,iostat=ierror) line
                    id_values = id_values + 1
                    line_array_temp(id_values) = line
                    if (ierror /= 0) exit
                end do

                if (ierror > 0) then
                    print "(A, I6)", "An error ocurred reading line ", id_values
                else
                    print "(A, I6, A)", "End of file reached. There were", id_values-1, " values in the file"
                end if
            else
                print "(A, I6)", "Error opening file: IOSTAT ", ierror
            end if

            close(io)
        end if

        allocate(line_array(1:id_values-1))

        line_array = line_array_temp(1:id_values-1)
    end subroutine import_data

    subroutine convert_data_to_curve(line_array, input_curve, refn_array)
        ! use derived_types
        character(len=*), dimension(0:), intent(in) :: line_array
        integer :: p
        real, dimension(:), allocatable :: U_knot
        real, dimension(:,:), allocatable :: ctrl_pts
        type(nurbs_curve), intent(out) :: input_curve
        character(len=1), dimension(:), allocatable, intent(out), optional :: refn_array

        character(:), allocatable :: format_line, second_format_line
        character(:), allocatable :: control_points_label, p_label, uknot_label, end_label, refinement_label
        integer :: i, label_counter
        integer :: control_points_flag, p_flag, uknot_flag, refinement_flag
        integer, dimension(:), allocatable :: label_position_temp, label_position

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
        allocate(label_position_temp(5))
        label_position_temp = 0

        control_points_flag = 0
        p_flag = 0
        uknot_flag = 0
        refinement_flag = 0

        do i = 1, size(line_array)
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

        allocate(label_position(label_counter))
        label_position = label_position_temp(label_counter)
        
        ! Reading control points
        if (control_points_flag == 1) then
            call read_real_matrix(line_array(label_position(1)+1:label_position(2)-1), ",", ctrl_pts)
        end if
        
        p = 0
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

    function find_geom_input_limits(string_a, string_b, string_array) result(string_limits)
        use utils
        character(len=*), intent(in) :: string_a, string_b
        character(len=*), dimension(:), intent(in) :: string_array
        integer, dimension(:,:), allocatable :: string_limits_temp, string_limits
        integer :: i, MAX_SIZE, counter_a, counter_b
        
        MAX_SIZE = size(string_array)
        allocate(string_limits_temp(MAX_SIZE,2))
        string_limits_temp = -1

        counter_a = 0
        counter_b = 0
        do i = 1, MAX_SIZE
            if (string_array(i) == string_a) then
                counter_a = counter_a + 1
                string_limits_temp(counter_a, 1) = i
            end if
            if (string_array(i) == string_b) then
                counter_b = counter_b + 1
                string_limits_temp(counter_b, 2) = i
            end if
        end do

        if (counter_a == counter_b) then
            allocate(string_limits(counter_a,2))
            string_limits = string_limits_temp(1:counter_a,:)
        end if
    end function find_geom_input_limits

    function find_string_in_array(string, string_array, i_start, i_end) result (idx)
        ! If string is not found, return -1
        character(len=*), intent(in) :: string
        character(len=*), dimension(:), intent(in) :: string_array
        integer, intent(in) :: i_start, i_end
        integer :: num_elements, i, idx

        idx = -1
        num_elements = size(string_array)

        do i = i_start, i_end
            if (string_array(i) == string) then
                idx = i
            end if
        end do
    end function find_string_in_array

    subroutine convert_data_to_surface(line_array, input_patches, refn_flag, refn_array, interf_flag, interf_var)
        ! use derived_types
        character(len=*), dimension(:), intent(in) :: line_array
        character(:), dimension(:), allocatable :: patch_input_array
        integer :: p, q
        real, dimension(:), allocatable :: U_knot, V_knot
        real, dimension(:,:), allocatable :: ctrl_pts
        type(nurbs_surface), dimension(:), allocatable, intent(out) :: input_patches
        character(len=1), dimension(:,:), allocatable, intent(out), optional :: refn_array
        type(interface_line), dimension(:), allocatable, intent(out), optional :: interf_var

        character(:), allocatable :: format_line, second_format_line
        character(:), allocatable :: start_geom_label, control_points_label
        character(:), allocatable :: p_label, q_label, uknot_label, vknot_label
        character(:), allocatable :: end_geom_label, refinement_label, interface_label
        integer, dimension(:,:), allocatable :: geom_input_limits
        
        integer :: i, i_start, i_end, num_patches
        integer :: start_geom_flag, end_geom_flag
        integer, intent(out) :: refn_flag, interf_flag

        integer :: idx_control_points, idx_p, idx_q, idx_UP, idx_VP, idx_refn
        integer :: idx_interf

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
        interface_label = "*INTERFACE"
        end_geom_label = "**END_GEOMETRY**"

        start_geom_flag = 0
        refn_flag = 0
        end_geom_flag = 0

        geom_input_limits = find_geom_input_limits(start_geom_label, end_geom_label, line_array)
        num_patches = size(geom_input_limits,1)

        allocate(input_patches(num_patches))

        do i = 1, num_patches
            i_start = geom_input_limits(i,1)
            i_end = geom_input_limits(i,2)
            patch_input_array = line_array(i_start:i_end)

            idx_control_points = find_string_in_array(control_points_label,line_array, i_start, i_end)
            idx_p = find_string_in_array(p_label,line_array, i_start, i_end)
            idx_q = find_string_in_array(q_label,line_array, i_start, i_end)
            idx_UP = find_string_in_array(uknot_label,line_array, i_start, i_end)
            idx_VP = find_string_in_array(vknot_label,line_array, i_start, i_end)
            idx_refn = find_string_in_array(refinement_label,line_array, i_start, i_end)
            idx_interf = find_string_in_array(interface_label,line_array, i_start, i_end)

            ! Reading control points
            if (idx_control_points /= -1) then
                call read_real_matrix(line_array(idx_control_points+1:idx_p-1), ",", ctrl_pts)
            end if

            print "(I3, I3)", lbound(ctrl_pts,1), ubound(ctrl_pts,1)
            
            ! Reading spline degree
            p = 0
            if (idx_p /= -1) then
                read (line_array(idx_p+1),*) p
            end if

            ! Reading spline degree
            q = 0
            if (idx_q /= -1) then
                read (line_array(idx_q+1),*) q
            end if

            ! Reading knot vector U
            if (idx_UP /= -1) then
                call read_real_row_array(line_array(idx_UP+1), ",", U_knot)
            end if

            ! Reading knot vector V
            if (idx_VP /= -1) then
                call read_real_row_array(line_array(idx_VP+1), ",", V_knot)
            end if

            ! Reading refining types
            if (idx_refn /= -1) then
                refn_flag = 1
                if (idx_interf == -1) then
                    call read_string_matrix(line_array(idx_refn+1:i_end-1), ",", refn_array)
                else
                    call read_string_matrix(line_array(idx_refn+1:idx_interf-1), ",", refn_array)
                end if
            end if

            ! Reading interface information
            if (idx_interf /= -1) then
                interf_flag = 1
                call read_interface_info(line_array(idx_interf+1:i_end-1), interf_var)
            end if

            size_1 = size(ctrl_pts, 1)
            size_2 = size(ctrl_pts, 2) - 1
            size_vec_1 = size(U_knot)
            size_vec_2 = size(V_knot)

            allocate(P_pts(0:size_1-1,0:size_2-1))
            allocate(w_pts(0:size_1-1,0))

            P_pts = ctrl_pts(1:size_1,1:size_2)
            w_pts = reshape(ctrl_pts(1:size_1,size_2+1), (/size_1, 1/))

            input_patches(i)%p = p
            input_patches(i)%q = q
            input_patches(i)%U_knot = U_knot
            input_patches(i)%V_knot = V_knot
            input_patches(i)%control_points = P_pts
            input_patches(i)%weight_points = w_pts
            input_patches(i)%refn_flag = refn_flag
            
            if (refn_flag == 1) then
                input_patches(i)%refn_input = refn_array
            end if

            deallocate(P_pts)
            deallocate(w_pts)
        end do
    end subroutine

    subroutine convert_data_to_solver(line_array, num_gauss_pts, kappa, bc_array, id_patch)
        ! bc_array: array with information for the enforcement of boundary conditions
        ! use derived_types
        character(len=*), dimension(:), intent(in) :: line_array
        integer, intent(out) :: num_gauss_pts
        real, intent(out) :: kappa

        character(:), allocatable :: format_line, second_format_line
        character(:), allocatable :: start_solver_label, intg_label
        character(:), allocatable :: kappa_label, bc_label
        character(:), allocatable :: end_solver_label
        integer :: i_start, i_end
        integer :: idx_intg, idx_kappa, idx_bc
        integer, dimension(:), allocatable, intent(out) :: id_patch
        
        type(boundary_condition), dimension(:), allocatable, intent(out) :: bc_array

        integer, dimension(:,:), allocatable :: solver_input_limits

        format_line = "(A)"
        second_format_line = "(A, I3)"

        start_solver_label = "**START_SOLVER**"
        intg_label = "*GAUSS_NUM"
        kappa_label = "*DIFFUSION"
        bc_label = "*ESSENTIAL_COND_DOFS"
        end_solver_label = "**END_SOLVER**"

        solver_input_limits = find_geom_input_limits(start_solver_label, end_solver_label, line_array)

        i_start = solver_input_limits(1,1)
        i_end = solver_input_limits(1,2)
        idx_intg = find_string_in_array(intg_label, line_array, i_start, i_end)
        idx_kappa = find_string_in_array(kappa_label, line_array, i_start, i_end)
        idx_bc = find_string_in_array(bc_label,line_array, i_start, i_end)

        ! Reading number of gauss points for numerical integration
        if (idx_intg /= -1) then
            read (line_array(idx_intg+1),*) num_gauss_pts
        end if
        
        ! Reading diffusion conductivity value
        if (idx_kappa /= -1) then
            read (line_array(idx_kappa+1),*) kappa
        end if

        ! Reading essential boundary conditions
        if (idx_bc /= -1) then
            call read_derived_type_matrix(line_array(idx_bc+1:i_end-1), bc_array, id_patch)
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

    subroutine read_derived_type_matrix(strings, derived_array, id_patches)
        character(len=*), dimension(:), intent(in) :: strings
        type(boundary_condition), dimension(:), allocatable, intent(out) :: derived_array
        integer, dimension(:), allocatable, intent(out) :: id_patches
        type(boundary_condition) :: temp

        integer :: num_patches, i

        num_patches = size(strings)/2
        allocate(derived_array(num_patches))
        allocate(id_patches(num_patches))
        id_patches = 0

        do i = 1, num_patches
            read (strings(2*i-1),*) id_patches(i)
            read (strings(2*i),*) temp
            derived_array(i) = temp
            print "('Patch ID: ', I3)", id_patches(i)
            print "('Direction: ', A)", temp%dir
            print "('Parameter: ', F4.2)", temp%UV_param
            print "('Prescribed value: ', F6.2)", temp%prescribed_val
        end do
    end subroutine read_derived_type_matrix

    subroutine read_interface_info(strings, interf_var)
        ! use derived_types
        character(len=*), dimension(:), intent(in) :: strings
        type(interface_line), dimension(:), allocatable, intent(out) :: interf_var
        type(interface_line) :: temp
        ! integer, parameter :: MAX_SIZE = 50
        integer :: num_values, i

        num_values = size(strings)
        allocate(interf_var(num_values))

        do i = 1, num_values
            read (strings(i),*) temp
            interf_var(i) = temp
            print "('Direction: ', A)", temp%dir
            print "('Parameter: ', F4.2)", temp%UV_param
        end do
    end subroutine

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