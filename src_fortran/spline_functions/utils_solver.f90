module utils_solver
    use utils
    use derived_types
    implicit none
contains
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

    subroutine compute_connectivity_matrices_2(nx, ny, p, q, i_patch, ni_dof, num_interfaces, elem_mat)
        integer, intent(in) :: nx, ny, p, q, i_patch, ni_dof
        integer, dimension(2), intent(in) :: num_interfaces
        integer :: i_elem, i, j, ii, jj, num_elem
        integer :: c1, num_nonzero, i_nonzero
        integer, dimension(:,:), allocatable, intent(out) :: elem_mat

        if (i_patch == 1) then
            num_elem = nx*ny
            num_nonzero = (p+1)*(q+1)

            allocate(elem_mat(0:num_nonzero,0:num_elem-1))
            elem_mat = 0

            i_elem = 0
            do j = 0, ny - 1
                do i = 0, nx - 1
                    elem_mat(0, i_elem) = i_elem
                    i_nonzero = 1
                    do jj = 0, q
                        do ii = 0, p
                            c1 = compute_global_dof(i + ii, j + jj, nx)
                            elem_mat(i_nonzero, i_elem) = c1
                            i_nonzero = i_nonzero + 1
                        end do
                    end do
                    i_elem = i_elem + 1
                end do
            end do
        else
            num_elem = nx*ny
            num_nonzero = (p+1)*(q+1)

            allocate(elem_mat(0:num_nonzero,0:num_elem-1))
            elem_mat = 0

            i_elem = 0
            do j = 0, ny - 1
                do i = 0, nx - 1
                    elem_mat(0, i_elem) = i_elem
                    i_nonzero = 1
                    if (i > 0) then
                        do jj = 0, q
                            do ii = 0, p
                                c1 = compute_multipatch_global_dof(i + ii, j + jj, nx, num_interfaces, ni_dof, i_patch)
                                elem_mat(i_nonzero, i_elem) = c1
                                i_nonzero = i_nonzero + 1
                            end do
                        end do
                    else
                        do jj = 0, q
                            ! Interface dofs
                            c1 = compute_global_dof(nx, j + jj, nx)
                            elem_mat(i_nonzero, i_elem) = c1
                            i_nonzero = i_nonzero + 1

                            ! Other patch dofs
                            c1 = compute_multipatch_global_dof(1, j + jj, nx, num_interfaces, ni_dof, i_patch)
                            elem_mat(i_nonzero, i_elem) = c1
                            i_nonzero = i_nonzero + 1
                        end do
                    end if
                    i_elem = i_elem + 1
                end do
            end do
        end if
        call print_integer_matrix(elem_mat)
        ! call print_integer_matrix(patch_nodes)
    end subroutine compute_connectivity_matrices_2

    subroutine compute_patch_nodes(num_patches, input_surf, interf_var, num_dofs_full, patch_nodes)
        use derived_types, only: nurbs_surface
        integer :: num_patches
        type(nurbs_surface), dimension(num_patches), intent(in) :: input_surf
        type(interface_boundary), dimension(:), allocatable, intent(in) :: interf_var
        integer :: i_patch, i_nodes, i, j, c1, c11, i_nodes_full, k_temp
        integer :: nu, nv, ni_dof, p, q, r, s, num_numbering
        integer, dimension(2) :: num_interfaces
        real, dimension(:), allocatable :: UP, VP
        integer, intent(out) :: num_dofs_full
        integer, dimension(:,:), allocatable, intent(out) :: patch_nodes

        num_dofs_full = 0
        num_numbering = 0

        do i_patch = 1, num_patches
            p = input_surf(i_patch)%p
            q = input_surf(i_patch)%q
            UP = input_surf(i_patch)%U_knot
            VP = input_surf(i_patch)%V_knot
            num_interfaces = input_surf(i_patch)%num_interfaces

            r = size(UP) - 1
            s = size(VP) - 1
            nu = r - p - 1
            nv = s - q - 1

            if (i_patch == 1) then
                num_dofs_full = num_dofs_full + (nu+1)*nv + nu
            else
                k_temp = (nu+1-num_interfaces(1))*(nv-num_interfaces(2)) + nu
                num_dofs_full = num_dofs_full + k_temp
            end if

            num_numbering = num_numbering + (nu+1)*nv + nu
        end do

        ! The first dof starts at 0, that is why it is added 1
        num_dofs_full = num_dofs_full + 1
        num_numbering = num_numbering + num_patches

        allocate(patch_nodes(0:num_numbering-1,0:2))
        i_nodes_full = 0

        patch_loop: do i_patch = 1, num_patches
            p = input_surf(i_patch)%p
            q = input_surf(i_patch)%q
            UP = input_surf(i_patch)%U_knot
            VP = input_surf(i_patch)%V_knot
        
            r = size(UP) - 1
            s = size(VP) - 1
            nu = r - p - 1
            nv = s - q - 1

            if (i_patch == 1) then
                i_nodes = 0
                do j = 0, nv
                    do i = 0, nu
                        c1 = compute_global_dof(i, j, nu)
                        patch_nodes(i_nodes_full, 0) = i_patch
                        patch_nodes(i_nodes_full, 1) = i_nodes
                        patch_nodes(i_nodes_full, 2) = c1
                        ! print "(I3, I3)", i_patch, i_nodes
                        i_nodes = i_nodes + 1
                        i_nodes_full = i_nodes_full + 1
                    end do
                end do
                ni_dof = ni_dof + (nu+1)*nv + nu
            else
                i_nodes = 0
                do j = 0, nv
                    do i = 0, nu
                        if (i > 0) then
                            c11 = compute_multipatch_global_dof(i, j, nu, num_interfaces, ni_dof, i_patch)
                        else
                            c11 = compute_global_dof(nu, j, nu)
                        end if
                        patch_nodes(i_nodes_full, 0) = i_patch
                        patch_nodes(i_nodes_full, 1) = i_nodes
                        patch_nodes(i_nodes_full, 2) = c11
                        ! print "(I3, I3)", i_patch, i_nodes
                        i_nodes = i_nodes + 1
                        i_nodes_full = i_nodes_full + 1
                    end do
                end do
                ni_dof = ni_dof + (nu)*nv + nu
            end if
        end do patch_loop

        call print_integer_matrix(patch_nodes)
    end subroutine compute_patch_nodes

    subroutine generate_node_numbering(input_surf, node_numbering)
        type(nurbs_surface), dimension(:), allocatable, intent(in) :: input_surf
        integer, dimension(:,:), allocatable, intent(out) :: node_numbering
        integer :: num_total_nodes, num_patches, i_node, patch_node_i
        integer :: num_patch_nodes, k_temp, i_nodes_full
        integer :: i_patch, p, q, r, s, nu, nv, i, j, i_start, i_end
        real, dimension(:), allocatable :: UP, VP
        integer, dimension(2) :: num_interfaces

        num_total_nodes = get_multipatch_number_nodes(input_surf)
        num_patches = size(input_surf)
        num_patch_nodes = 0
        i_nodes_full = 0

        allocate(node_numbering(0:num_total_nodes-1,0:2))
        node_numbering = 0

        patch_loop: do i_patch = 1, num_patches
            p = input_surf(i_patch)%p
            q = input_surf(i_patch)%q
            UP = input_surf(i_patch)%U_knot
            VP = input_surf(i_patch)%V_knot
            num_interfaces = input_surf(i_patch)%num_interfaces
        
            r = size(UP) - 1
            s = size(VP) - 1
            nu = r - p - 1
            nv = s - q - 1

            i_node = 0

            ! ITERATE THROUGH INTERFACES AND SET THE INDEX BOUNDS
            ! USE IMPLEMENTATION FROM get_boundary_conditions_dof
            ! AS REFERENCE

            i_start = 0 + num_interfaces(1)
            i_end = nu - num_interfaces(1)

            do j = 0, nv
                do i = 0, nu
                    ! patch_node_i = compute_global_dof(i, j, nu)
                    if (i_patch > 1) then
                        if (i /= i_start) then
                            patch_node_i = compute_multipatch_global_dof(i, j, nu, num_interfaces, num_patch_nodes, i_patch)
                        else
                            patch_node_i = compute_multipatch_global_dof(0, j, nu, num_interfaces, num_patch_nodes, 1)
                        end if
                    else
                        patch_node_i = compute_multipatch_global_dof(i, j, nu, num_interfaces, num_patch_nodes, i_patch)
                    end if
                    node_numbering(i_nodes_full, 0) = i_patch
                    node_numbering(i_nodes_full, 1) = i_node
                    node_numbering(i_nodes_full, 2) = patch_node_i
                    ! print "(I3, I3)", i_patch, i_node
                    i_node = i_node + 1
                    i_nodes_full = i_nodes_full + 1
                end do
            end do
            k_temp = (nu+1-num_interfaces(1))*(nv-num_interfaces(2)) + nu
            num_patch_nodes = num_patch_nodes + k_temp
        end do patch_loop
        call print_integer_matrix(node_numbering)
    end subroutine generate_node_numbering

    function compute_global_dof(i, j, nu) result(dof)
        integer, intent(in) :: i, j, nu
        integer :: dof
        dof = j*(nu+1) + i
    end function compute_global_dof

    function compute_multipatch_global_dof(i, j, nu, num_interfaces, ni_dof, i_patch) result(dof)
        integer, intent(in) :: i, j, nu, ni_dof, i_patch
        integer, dimension(2), intent(in) :: num_interfaces
        integer :: dof
        if (i_patch > 1) then
            dof = j*(nu+1-num_interfaces(1)) + i + ni_dof
        else
            dof = j*(nu+1-num_interfaces(1)) + i
        end if
    end function compute_multipatch_global_dof

    function get_multipatch_number_nodes(input_surf) result(num_total_nodes)
        integer :: num_patches
        type(nurbs_surface), dimension(:), allocatable, intent(in) :: input_surf
        integer :: num_total_nodes

        integer :: i_patch, p, q, r, s, nu, nv, k_temp
        integer, dimension(2) :: num_interfaces
        real, dimension(:), allocatable :: UP, VP
        
        num_patches = size(input_surf)
        num_total_nodes = 0

        do i_patch = 1, num_patches
            p = input_surf(i_patch)%p
            q = input_surf(i_patch)%q
            UP = input_surf(i_patch)%U_knot
            VP = input_surf(i_patch)%V_knot
            num_interfaces = input_surf(i_patch)%num_interfaces

            r = size(UP) - 1
            s = size(VP) - 1
            nu = r - p - 1
            nv = s - q - 1

            if (i_patch == 1) then
                num_total_nodes = num_total_nodes + (nu+1)*nv + nu
            else
                ! k_temp = (nu+1-num_interfaces(1))*(nv-num_interfaces(2)) + nu
                k_temp = (nu+1)*(nv+1)
                num_total_nodes = num_total_nodes + k_temp
            end if
        end do

        ! The first dof starts at 0, that is why it is added 1
        num_total_nodes = num_total_nodes + 1
    end function get_multipatch_number_nodes

    subroutine get_boundary_conditions_dof(input_surf, bc_array, id_patches, id_disp, u_pres, ctrl_pts_pres)
        ! ctrl_pts_pres: Control points with prescribed essential boundary conditions
        use derived_types
        integer :: p, q
        real, dimension(:), allocatable :: UP, VP
        integer, dimension(:), allocatable, intent(in) :: id_patches
        type(nurbs_surface), dimension(:), allocatable, intent(in) :: input_surf
        type(boundary_condition), dimension(:), allocatable, intent(in) :: bc_array
        
        integer, dimension(:), allocatable, intent(out) :: id_disp
        real, dimension(:), allocatable, intent(out) :: u_pres
        real, dimension(:,:), allocatable, intent(out) :: ctrl_pts_pres
        real :: u_val

        integer :: i, r, s, iu, jv, nu, nv, id, i_temp_1, id_i, i_bc
        integer :: i_patch, num_patches, ni_dof, num_bc
        integer, dimension(:), allocatable :: i_arr, temp_1
        integer, dimension(2) :: num_interfaces
        real, dimension(:), allocatable :: temp_2
        real, dimension(:,:), allocatable :: temp_3, ctrl_pts_i
        real, dimension(:,:), allocatable :: P_pts

        integer, parameter :: MAX_SIZE = 5000
        allocate(temp_1(0:MAX_SIZE-1))
        allocate(temp_2(0:MAX_SIZE-1))
        allocate(temp_3(0:MAX_SIZE-1,0:2))
        temp_1 = 0
        temp_2 = 0.0
        temp_3 = 0.0
        i_temp_1 = 0

        ! num_patches = size(input_surf)
        num_bc = size(bc_array)

        ni_dof = 0

        patch_loop: do i_bc = 1, num_bc
            print "(A, I3)", "BC #", i_bc
            i_patch = id_patches(i_bc)
            print "(A, I3)", "Patch #", i_patch
            p = input_surf(i_patch)%p
            q = input_surf(i_patch)%q
            UP = input_surf(i_patch)%U_knot
            VP = input_surf(i_patch)%V_knot
            P_pts = input_surf(i_patch)%control_points
            num_interfaces = input_surf(i_patch)%num_interfaces
            
            allocate(ctrl_pts_i(0:size(P_pts,1)-1,0:size(P_pts,2)-1))
            ctrl_pts_i = input_surf(i_patch)%control_points

            r = size(UP) - 1
            s = size(VP) - 1
            nu = r - p - 1
            nv = s - q - 1

            ! print "(I3, I3)", lbound(ctrl_pts_i,1), ubound(ctrl_pts_i,1)

            if (bc_array(i_bc)%dir == "U") then
                allocate(i_arr(0:nu))
                u_val = bc_array(i_bc)%prescribed_val
                i_arr = (/(i, i = 0, nu)/)
                if (abs(bc_array(i_bc)%UV_param) < 1e-5) then
                    jv = 0
                end if
                if (abs(bc_array(i_bc)%UV_param-1.0) < 1e-5) then
                    jv = nv
                end if

                do i = 0, nu
                    id = compute_global_dof(i_arr(i),jv,nu)
                    id_i = compute_multipatch_global_dof(i_arr(i), jv, nu, num_interfaces, ni_dof, i_patch)
                    print "(A, I3, A, I3, A, F6.1)", "Patch ID:", i_patch, "|Global dof:", id_i, " |BC value:", u_val
                    temp_1(i_temp_1) = id_i
                    temp_2(i_temp_1) = u_val
                    temp_3(i_temp_1,:) = ctrl_pts_i(id,:)
                    i_temp_1 = i_temp_1 + 1
                end do
                deallocate(i_arr)
            else if (bc_array(i_bc)%dir == "V") then
                allocate(i_arr(0:nv))
                u_val = bc_array(i_bc)%prescribed_val
                i_arr = (/(i, i = 0, nv)/)
                if (abs(bc_array(i_bc)%UV_param) < 1e-5) then
                    iu = 0
                end if
                if (abs(bc_array(i_bc)%UV_param-1.0) < 1e-5) then
                    iu = nu
                end if

                ! print "(I3, I3)", ni_dof, i_patch

                do i = 0, nv
                    id = compute_global_dof(iu,i_arr(i),nu)
                    id_i = compute_multipatch_global_dof(iu, i_arr(i), nu, num_interfaces, ni_dof, i_patch)
                    print "(A, I3, A, I3, A, F6.1)", "Patch ID:", i_patch, "|Global dof:", id_i, " |BC value:", u_val
                    ! print "(F6.3)", ctrl_pts_i(id,:)
                    temp_1(i_temp_1) = id_i
                    temp_2(i_temp_1) = u_val
                    temp_3(i_temp_1,:) = ctrl_pts_i(id,:)
                    i_temp_1 = i_temp_1 + 1
                end do
                deallocate(i_arr)
            end if
            ni_dof = ni_dof + size(ctrl_pts_i,1) - 1
            deallocate(ctrl_pts_i)
        end do patch_loop

        allocate(id_disp(i_temp_1))
        allocate(u_pres(i_temp_1))
        allocate(ctrl_pts_pres(i_temp_1,0:2))
        id_disp = temp_1(0:i_temp_1-1)
        u_pres = temp_2(0:i_temp_1-1)
        ctrl_pts_pres = temp_3(0:i_temp_1-1,:)
        ! call print_matrix(ctrl_pts_pres)

        ! DOFS FOR MULTIPATCH SCHEMES
        call print_row_vector_intg(id_disp)
    end subroutine get_boundary_conditions_dof
end module utils_solver