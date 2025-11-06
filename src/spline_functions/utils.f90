module utils
    implicit none
contains
    subroutine linspace_real(min_val, max_val, num_points, ndim, vec)
        real, intent(in) :: min_val, max_val
        integer, intent(in) :: num_points, ndim
        real, dimension(:,:), intent(out), allocatable :: vec
        integer :: i

        allocate(vec(num_points, ndim))

        vec = 0

        do i = 1, num_points
            vec(i,1) = ((max_val - min_val)/(num_points-1))*(i-1) + min_val
        end do
    end subroutine linspace_real

    subroutine linspace_intg(min_val, max_val, step, vec)
        integer, intent(in) :: min_val, max_val
        integer, intent(in) :: step
        integer, dimension(:), intent(out), allocatable :: vec
        integer :: i, num_points

        num_points = (max_val-min_val)/step + 1
        allocate(vec(num_points))

        do i = 1, num_points
            vec(i) = min_val + step*(i-1)
        end do
    end subroutine linspace_intg

    subroutine print_integer(num)
        integer, intent(in) :: num
        print "(I3)", num
    end subroutine print_integer

    subroutine print_real(num)
        real, intent(in) :: num
        print "(F6.3)", num
    end subroutine print_real
    
    subroutine print_matrix(mat)
        real, intent(in) :: mat(:,:)
        integer :: num_rows, num_cols, i, j

        num_rows = size(mat,1)
        num_cols = size(mat,2)

        print *, "===================="
        do i = 1, num_rows
            ! print '(*(3X, f8.5))', (mat(i,j), j=1,num_cols)
            print '(*(3X, f15.4))', (mat(i,j), j=1,num_cols)
        end do
        print *, "===================="
    end subroutine print_matrix

    subroutine print_row_vector(row_vec)
        real, intent(in) :: row_vec(:)
        integer :: i, num_cols

        num_cols = size(row_vec,1)

        print *, "===================="
        print '(*(3X, F8.5))', (row_vec(i), i=1,num_cols)
        print *, "===================="
    end subroutine print_row_vector

    subroutine print_row_vector_intg(row_vec)
        integer, intent(in) :: row_vec(:)
        integer :: i, num_cols

        num_cols = size(row_vec,1)

        print *, "===================="
        print '(*(3X, I5))', (row_vec(i), i=1,num_cols)
        print *, "===================="
    end subroutine print_row_vector_intg

    subroutine print_row_vector_bool(row_vec)
        logical, intent(in) :: row_vec(:)
        integer :: i, num_cols

        num_cols = size(row_vec,1)

        print *, "===================="
        print '(*(3X, L5))', (row_vec(i), i=1,num_cols)
        print *, "===================="
    end subroutine print_row_vector_bool

    subroutine print_column_vector(col_vec)
        real, intent(in) :: col_vec(:,:)
        integer :: num_rows, i

        num_rows = size(col_vec,1)

        print *, "===================="
        do i = 1, num_rows
            print '((3X, f15.3))', col_vec(i, 1)
        end do
        print *, "===================="
    end subroutine print_column_vector

    subroutine print_string_matrix(str_mat)
        character(len=*), dimension(:,:), intent(in) :: str_mat
        integer :: num_rows, num_cols, i, j

        num_rows = size(str_mat,1)
        num_cols = size(str_mat,2)

        print *, "===================="
        do i = 1, num_rows
            ! print '(*(3X, f8.5))', (mat(i,j), j=1,num_cols)
            print '(*(1X, A))', (str_mat(i,j), j=1,num_cols)
        end do
        print *, "===================="
    end subroutine print_string_matrix

    subroutine export_matrix(mat, file_name)
        ! This function is to be generalized for more than two columns
        real, intent(in) :: mat(:,:)
        integer :: num_rows, num_cols, i, io, num_characters
        logical :: exists
        character(:), intent(in), allocatable :: file_name

        num_rows = size(mat,1)
        num_cols = size(mat,2)
        num_characters = len(file_name)

        inquire(file=file_name, exist=exists)
        if (exists) then
            open(newunit=io, file=file_name, status="old", action="write")
            do i = 1, num_rows
                write(io,"(F6.3, A, F8.3)") mat(i,1), ", ", mat(i,2)
            end do
            close(io)
        else
            open(newunit=io, file=file_name, status="new", action="write")
            do i = 1, num_rows
                write(io,"(F6.3, A, F8.3)") mat(i,1), ", ", mat(i,2)
            end do
            close(io)
        end if

    end subroutine export_matrix

    subroutine compute_element_midvalues(array, p, midvalues)
        real, dimension(0:), intent(in) :: array
        integer, intent(in) :: p
        real, dimension(:), allocatable, intent(out) :: midvalues
        integer :: r

        r = size(array) - 1

        allocate(midvalues(0:r-1))

        midvalues = 0.5*(array(p:r-p-1) + array(p+1:r-p))
    end subroutine compute_element_midvalues

    subroutine calculate_surface_elements(Uarray, Varray, surface_elem_bounds)
        real, dimension(0:), intent(in) :: Uarray
        real, dimension(0:), intent(in) :: Varray
        integer :: i, j, U_size, V_size, num_elems_U, num_elems_V
        real, dimension(:), allocatable :: unique_Uarray, unique_Varray, temp_Uarray, temp_Varray
        integer :: i_unique, i_elems, num_bounds
        real, dimension(:,:), allocatable, intent(out) :: surface_elem_bounds

        U_size = size(Uarray) - 1
        V_size = size(Varray) - 1
        num_elems_U = 0
        num_elems_V = 0

        allocate(temp_Uarray(0:U_size))
        allocate(temp_Varray(0:V_size))
        temp_Uarray = 0.
        temp_Varray = 0.

        i_unique = 0

        do i = 0, U_size-1
            if (abs(Uarray(i+1) - Uarray(i)) > 1e-5) then
                num_elems_U = num_elems_U + 1
                temp_Uarray(i_unique) = Uarray(i)
                temp_Uarray(i_unique+1) = Uarray(i+1)
                i_unique = i_unique + 1
            end if
        end do

        i_unique = 0

        do i = 0, V_size-1
            if (abs(Varray(i+1) - Varray(i)) > 1e-5) then
                num_elems_V = num_elems_V + 1
                temp_Varray(i_unique) = Varray(i)
                temp_Varray(i_unique+1) = Varray(i+1)
                i_unique = i_unique + 1
            end if
        end do

        print "(A, I3)", "Number of elements in U direction: ", num_elems_U
        print "(A, I3)", "Number of elements in V direction: ", num_elems_V

        allocate(unique_Uarray(0:num_elems_U))
        allocate(unique_Varray(0:num_elems_V))
        
        unique_Uarray = temp_Uarray(0:num_elems_U)
        unique_Varray = temp_Varray(0:num_elems_V)

        num_bounds = num_elems_U*num_elems_V

        allocate(surface_elem_bounds(0:num_bounds-1,0:3))

        i_elems = 0
        
        do j = 0, num_elems_V-1
            do i = 0, num_elems_U-1
                surface_elem_bounds(i_elems, 0) = unique_Uarray(i)
                surface_elem_bounds(i_elems, 1) = unique_Varray(j)
                surface_elem_bounds(i_elems, 2) = unique_Uarray(i+1)
                surface_elem_bounds(i_elems, 3) = unique_Varray(j+1)
                i_elems = i_elems + 1
            end do
        end do
    end subroutine calculate_surface_elements

    function inverse_matrix_2(mat) result(inv_mat)
        real, dimension(0:1,0:1), intent(in) :: mat
        real, dimension(0:1,0:1) :: inv_mat
        real :: det_mat

        inv_mat = 0.

        det_mat = mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0)
        inv_mat(0,0) = (1./det_mat)*mat(1,1)
        inv_mat(0,1) = -(1./det_mat)*mat(0,1)
        inv_mat(1,0) = -(1./det_mat)*mat(1,0)
        inv_mat(1,1) = (1./det_mat)*mat(0,0)
    end function inverse_matrix_2

    function inverse_matrix_3(mat) result(inv_mat)
        real, dimension(0:2,0:2), intent(in) :: mat
        real, dimension(0:2,0:2) :: inv_mat
        real :: A, B, C, D, E, F, G, H, I, det_mat

        A = (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1))
        B = -(mat(1,0)*mat(2,2) - mat(1,2)*mat(2,0))
        C = (mat(1,0)*mat(2,1) - mat(1,1)*mat(2,0))
        D = -(mat(0,1)*mat(2,2) - mat(0,2)*mat(2,1))
        E = (mat(0,0)*mat(2,2) - mat(0,2)*mat(2,0))
        F = -(mat(0,0)*mat(2,1) - mat(0,1)*mat(2,0))
        G = (mat(0,1)*mat(1,2) - mat(0,2)*mat(1,1))
        H = -(mat(0,0)*mat(1,2) - mat(0,2)*mat(1,0))
        I = (mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0))

        det_mat = mat(0,0)*A + mat(0,1)*B + mat(0,2)*C

        inv_mat(0,0) = (1./det_mat)*A
        inv_mat(1,0) = (1./det_mat)*B
        inv_mat(2,0) = (1./det_mat)*C
        inv_mat(0,1) = (1./det_mat)*D
        inv_mat(1,1) = (1./det_mat)*E
        inv_mat(2,1) = (1./det_mat)*F
        inv_mat(0,2) = (1./det_mat)*G
        inv_mat(1,2) = (1./det_mat)*H
        inv_mat(2,2) = (1./det_mat)*I
    end function inverse_matrix_3
end module utils