module utils
    implicit none
contains
    subroutine linspace(min_val, max_val, num_points, ndim, vec)
        real, intent(in) :: min_val, max_val
        integer, intent(in) :: num_points, ndim
        real, dimension(:,:), intent(out), allocatable :: vec
        integer :: i

        allocate(vec(num_points, ndim))

        vec = 0

        do i = 1, num_points
            vec(i,1) = ((max_val - min_val)/(num_points-1))*(i-1) + min_val
        end do
    end subroutine linspace

    subroutine print_matrix(mat)
        real, intent(in) :: mat(:,:)
        integer :: num_rows, num_cols, i, j

        num_rows = size(mat,1)
        num_cols = size(mat,2)

        print *, "===================="
        do i = 1, num_rows
            print '(*(3X, f8.5))', (mat(i,j), j=1,num_cols)
        end do
        print *, "===================="
    end subroutine print_matrix

    subroutine print_row_vector(row_vec)
        real, intent(in) :: row_vec(:)
        integer :: i, num_cols

        num_cols = size(row_vec,1)

        print *, "===================="
        print '(*(3X, f8.5))', (row_vec(i), i=1,num_cols)
        print *, "===================="
    end subroutine print_row_vector

    subroutine print_column_vector(col_vec)
        real, intent(in) :: col_vec(:,:)
        integer :: num_rows, i

        num_rows = size(col_vec,1)

        print *, "===================="
        do i = 1, num_rows
            print '((3X, f9.5))', col_vec(i, 1)
        end do
        print *, "===================="
    end subroutine print_column_vector

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
end module utils