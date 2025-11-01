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
end module utils