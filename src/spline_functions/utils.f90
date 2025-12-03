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

    subroutine print_integer_matrix(mat)
        integer, intent(in) :: mat(:,:)
        integer :: num_rows, num_cols, i, j

        num_rows = size(mat,1)
        num_cols = size(mat,2)

        print *, "===================="
        do i = 1, num_rows
            ! print '(*(3X, f8.5))', (mat(i,j), j=1,num_cols)
            print '(*(3X, I5))', (mat(i,j), j=1,num_cols)
        end do
        print *, "===================="
    end subroutine print_integer_matrix

    subroutine print_row_vector(row_vec)
        real, intent(in) :: row_vec(:)
        integer :: i, num_cols

        num_cols = size(row_vec,1)

        print *, "===================="
        print '(*(3X, F10.4))', (row_vec(i), i=1,num_cols)
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

    function get_unique_values(full_arr) result(unique_arr)
        real, dimension(0:), intent(in) :: full_arr
        real, dimension(:), allocatable :: temp_arr, unique_arr
        integer :: i, arr_size, i_unique

        arr_size = size(full_arr)
        i_unique = 0

        allocate(temp_arr(0:arr_size-1))

        do i = 0, arr_size-2
            if (abs(full_arr(i+1) - full_arr(i)) > 1e-5) then
                temp_arr(i_unique) = full_arr(i)
                temp_arr(i_unique+1) = full_arr(i+1)
                i_unique = i_unique + 1
            end if
        end do

        allocate(unique_arr(0:i_unique))
        unique_arr = temp_arr(0:i_unique)
    end function get_unique_values

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
        real, dimension(:), allocatable :: unique_Uarray, unique_Varray
        integer :: i_elems, num_bounds
        real, dimension(:,:), allocatable, intent(out) :: surface_elem_bounds

        U_size = size(Uarray) - 1
        V_size = size(Varray) - 1
        num_elems_U = 0
        num_elems_V = 0

        allocate(unique_Uarray(0:num_elems_U))
        allocate(unique_Varray(0:num_elems_V))
        
        unique_Uarray = get_unique_values(Uarray)
        unique_Varray = get_unique_values(Varray)

        num_elems_U = size(unique_Uarray) - 1
        num_elems_V = size(unique_Varray) - 1

        print "(A, I3)", "Number of elements in U direction: ", num_elems_U
        print "(A, I3)", "Number of elements in V direction: ", num_elems_V

        num_bounds = num_elems_U*num_elems_V

        allocate(surface_elem_bounds(0:num_bounds-1,0:3))

        i_elems = 0
        
        do j = 1, num_elems_V
            do i = 1, num_elems_U
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

    subroutine compute_connectivity_matrices(p, q, n, m, IEN, INC)
        ! Adapted from Algorithm 7. Appendix A Connectivity arrays from
        ! Isogeometric Analysis, Cottrell, Hughes and Bazilevs (2009)
        ! --------------------------------------------------------------------
        ! Variables
        ! p: polynomial order in U direction
        ! q: polynomial order in V direction
        ! n: number of functions in U direction
        ! m: number of functions in V direction
        integer, intent(in) :: p, q, n, m
        integer, dimension(:,:), allocatable, intent(out) :: IEN, INC

        integer :: num_elements, num_global_functions, num_local_functions
        integer :: e, A, BB, b, i, j, iloc, jloc

        num_elements = (n-p)*(m-q)
        num_global_functions = n*m
        num_local_functions = (p+1)*(q+1)

        allocate(INC(num_global_functions,2))
        allocate(IEN(num_elements,num_elements))
        INC = 0
        IEN = 0

        e = 0
        A = 0
        BB = 0
        b = 0

        do j = 1, m
            do i = 1, n
                A = A+1
                INC(A,1) = i
                INC(A,2) = j
                if (i >= p+1) then
                    if (j >= q+1) then
                        e = e+1
                        do jloc = 0, q
                            do iloc = 0, p
                                BB = A - jloc*n - iloc
                                b = jloc*(p+1) + iloc + 1
                                IEN(b,e) = BB
                            end do
                        end do
                    end if
                end if
            end do
        end do
    end subroutine compute_connectivity_matrices

    function compute_global_dof(i,j,nu) result(dof)
        integer, intent(in) :: i,j,nu
        integer :: dof
        dof = j*(nu+1) + i
    end function compute_global_dof

    function compute_parametric_boundary_points(arr1, arr2, arr3, arr4) result(bound_param_pts)
        ! arr_i structure: Umin, Vmin, Umax, Vmax
        !
        !
        real, dimension(:), intent(in) :: arr1, arr2, arr3, arr4
        real, dimension(:,:), allocatable :: bound_param_pts
        integer :: num_points, i_pt, i
        real :: u, v, u_max, u_min, v_max, v_min

        num_points = 41

        allocate(bound_param_pts(0:(4*num_points)-1,0:1))
        bound_param_pts = 0.0
        i_pt = 0

        ! u = ((u_max - u_min)/num_points)*(i-1) + u_min
        ! v = ((v_max - v_min)/num_points)*(i-1) + v_min

        u_min = arr1(1)
        v_min = arr1(2)
        u_max = arr1(3)
        v_max = arr1(4)

        first_boundary: do i = 1, num_points
            u = ((u_max - u_min)/num_points)*(i-1) + u_min
            v = v_min
            bound_param_pts(i_pt,0) = u
            bound_param_pts(i_pt,1) = v
            i_pt = i_pt + 1
        end do first_boundary

        u_min = arr2(1)
        v_min = arr2(2)
        u_max = arr2(3)
        v_max = arr2(4)

        second_boundary: do i = 1, num_points
            u = u_max
            v = ((v_max - v_min)/num_points)*(i-1) + v_min
            bound_param_pts(i_pt,0) = u
            bound_param_pts(i_pt,1) = v
            i_pt = i_pt + 1
        end do second_boundary

        u_min = arr3(1)
        v_min = arr3(2)
        u_max = arr3(3)
        v_max = arr3(4)

        third_boundary: do i = 1, num_points
            u = ((u_max - u_min)/num_points)*(i-1) + u_min
            v = v_max
            bound_param_pts(i_pt,0) = u
            bound_param_pts(i_pt,1) = v
            i_pt = i_pt + 1
        end do third_boundary

        u_min = arr4(1)
        v_min = arr4(2)
        u_max = arr4(3)
        v_max = arr4(4)

        fourth_boundary: do i = 1, num_points
            u = u_min
            v = ((v_max - v_min)/num_points)*(i-1) + v_min
            bound_param_pts(i_pt,0) = u
            bound_param_pts(i_pt,1) = v
            i_pt = i_pt + 1
        end do fourth_boundary
    end function compute_parametric_boundary_points

    subroutine get_boundary_conditions_dof(surf, bc_array, id_disp, u_pres)
        use derived_types
        integer :: p, q
        real, dimension(:), allocatable :: UP, VP
        type(nurbs_surface), intent(in) :: surf
        type(boundary_condition), dimension(:), allocatable, intent(in) :: bc_array
        integer, dimension(:), allocatable, intent(out) :: id_disp
        real, dimension(:), allocatable, intent(out) :: u_pres

        integer :: i, j, r, s, iu, jv, nu, nv, id, num_bc, i_temp_1
        integer, dimension(:), allocatable :: i_arr, temp_1
        real, dimension(:), allocatable :: temp_2

        integer, parameter :: MAX_SIZE = 5000

        p = surf%p
        q = surf%q
        UP = surf%U_knot
        VP = surf%V_knot

        r = size(UP) - 1
        s = size(VP) - 1
        nu = r - p - 1
        nv = s - q - 1

        num_bc = size(bc_array)

        allocate(temp_1(0:MAX_SIZE-1))
        allocate(temp_2(0:MAX_SIZE-1))
        temp_1 = 0
        temp_2 = 0.0
        i_temp_1 = 0

        do j = 1, num_bc
            if (bc_array(j)%dir == "U") then
                allocate(i_arr(0:nu))
                i_arr = (/(i, i = 0, nu)/)
                if (abs(bc_array(j)%UV_param) < 1e-5) then
                    jv = 0
                end if
                if (abs(bc_array(j)%UV_param-1.0) < 1e-5) then
                    jv = nv
                end if
    
                do i = 0, nu
                    id = compute_global_dof(i_arr(i),jv,nu)
                    ! print "(A, I3, A, F6.1)", "Global dof:", id, " |BC value:", bc_arr(j)
                    temp_1(i_temp_1) = id
                    temp_2(i_temp_1) = bc_array(j)%prescribed_val
                    i_temp_1 = i_temp_1 + 1
                end do
                deallocate(i_arr)
            else if (bc_array(j)%dir == "V") then
                allocate(i_arr(0:nv))
                i_arr = (/(i, i = 0, nv)/)
                if (abs(bc_array(j)%UV_param) < 1e-5) then
                    iu = 0
                end if
                if (abs(bc_array(j)%UV_param-1.0) < 1e-5) then
                    iu = nu
                end if
                do i = 0, nv
                    id = compute_global_dof(iu,i_arr(i),nu)
                    ! print "(A, I3, A, F6.1)", "Global dof:", id, " |BC value:", bc_arr(j)
                    temp_1(i_temp_1) = id
                    temp_2(i_temp_1) = bc_array(j)%prescribed_val
                    i_temp_1 = i_temp_1 + 1
                end do
                deallocate(i_arr)
            end if
        end do

        allocate(id_disp(i_temp_1-2))
        allocate(u_pres(i_temp_1-2))
        id_disp = temp_1(0:i_temp_1-2)
        u_pres = temp_2(0:i_temp_1-2)
    end subroutine get_boundary_conditions_dof

    subroutine print_nurbs_surface_info(surf)
        use derived_types
        type(nurbs_surface), intent(in) :: surf

        print "(A)", "Surface control points"
        call print_matrix(surf%control_points)
        print "(A)", "Surface weight points"
        call print_matrix(surf%weight_points)
        print "(A)", "U knot vector"
        call print_row_vector(surf%U_knot)
        print "(A)", "V knot vector"
        call print_row_vector(surf%V_knot)
        print "(A, I3, A, I3)", "p degree:", surf%p, " |  q degree:", surf%q
    end subroutine print_nurbs_surface_info

    function find_string_in_array(string, string_array) result (idx)
        ! If string is not found, return -1
        character(len=*), intent(in) :: string
        character(len=*), dimension(0:), intent(in) :: string_array
        integer :: num_elements, i, idx

        idx = -1
        num_elements = size(string_array)

        do i = 0, num_elements-1
            if (string_array(i) == string) then
                idx = i
            end if
        end do
    end function find_string_in_array
end module utils