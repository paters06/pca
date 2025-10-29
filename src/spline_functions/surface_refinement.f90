module surface_refinement
    use utils
    implicit none
contains
    subroutine surface_knot_insertion(p, UP, q, VP, Pw_net, uv, dir, nq, UQ, mq, VQ, Qw_net)
        ! Algorithm A5.3 form the NURBS Book
        ! Surface knot insertion
        ! ---------------------------------------------------------------------
        ! Input: np, p, UP, mp, q, VP, Pw, dir, uv, k, s, r
        ! Output: nq, UQ, mq, VQ, Qw
        ! ---------------------------------------------------------------------
        ! Variables:
        ! uv is the u-v knot to insert
        ! np+1 are the rows of the Pw_net control net
        ! mp+1 are the columns of the Pw_net control net
        ! nq+1 are the rows of the Qw_net control net
        ! mq+1 are the columns of the Qw_net control net
        ! UP and VP are the U and V knot vectors before the insertion
        ! UQ and VQ are the U and V knot vectors after the insertion
        ! Pw_net and Qw_net are the control nets before and after insertion
        ! k: knot span
        ! s: knot multiplicity
        ! r: number of insertions
        use bspline_basis_functions, only: find_span

        integer, intent(in) :: p, q
        real, intent(in) :: uv
        real, dimension(0:), intent(in) :: UP, VP
        real, dimension(0:,0:,0:), intent(in) :: Pw_net
        character(len=*), intent(in) :: dir

        integer, intent(out) :: nq, mq
        real, dimension(:), allocatable, intent(out) :: UQ, VQ
        real, dimension(:,:,:), allocatable, intent(out) :: Qw_net

        integer :: np, mp, r, s, k, rp, sp
        integer :: j, L, i, row
        real, dimension(:,:), allocatable :: Rw
        real, dimension(:,:), allocatable :: alpha

        r = 1
        s = 0
        
        U_dir_refn: if (dir == "U") then
            np = size(Pw_net,1) - 1
            mp = size(Pw_net,2) - 1

            rp = np + p + 1

            nq = np + r
            mq = mp

            allocate(alpha(0:p-s-1,r))
            allocate(UQ(0:nq+p+1))
            allocate(Rw(0:p,0:3))
            allocate(Qw_net(0:np+1,0:mp,0:3))
            
            call find_span(rp, p, uv, UP, k)
            
            !! Load u-vector as in A5.1
            do i = 0, k
                UQ(i) = UP(i)
            end do
    
            do i = 1, r
                UQ(k+i) = uv
            end do
    
            do i = k+1, rp
                UQ(i+r) = UP(i)
            end do

            VQ = VP
            ! Save the alphas
            do j = 1, r
                L = k-p+j
                do i = 0, p-j+s
                    alpha(i,j) = (uv-UP(L+i))/(UP(i+k+1)-UP(L+i))
                end do
            end do
            ! For each row do
            do row = 0, mp
                ! Save unaltered control points
                do i = 0, k-p
                    Qw_net(i,row,:) = Pw_net(i,row,:)
                end do

                do i = k-s, np
                    print *, i
                    Qw_net(i+r,row,:) = Pw_net(i,row,:)
                end do

                ! Load auxiliary control points
                do i = 0, p-s
                    Rw(i,:) =  Pw_net(k-p+i,row,:)
                end do

                ! Insert the knor r times
                do j = 1, r
                    L = k-p+j
                    do i = 0, p-j-s
                        Rw(i,:) = alpha(i,j)*Rw(i+1,:) + (1.0-alpha(i,j))*Rw(i,:)
                    end do
                    Qw_net(L,row,:) = Rw(0,:)
                    Qw_net(k+r-j-s,row,:) = Rw(p-j-s,:)
                end do

                ! Load the remaining control points
                do i=L+1, k-s-1
                    Qw_net(i,row,:) = Rw(i-L,:)
                end do
            end do
        else if (dir == "V") then
            np = size(Pw_net,1) - 1
            mp = size(Pw_net,2) - 1

            nq = np
            mq = mp + r

            sp = np + q + 1

            allocate(alpha(0:p-s-1,r))
            allocate(VQ(0:mq+q+1))
            allocate(Rw(0:q,0:3))
            allocate(Qw_net(0:np,0:mp+1,0:3))
            
            call find_span(sp, q, uv, VP, k)
            
            !! Load u-vector as in A5.1
            do i = 0, k
                VQ(i) = VP(i)
            end do
    
            do i = 1, r
                VQ(k+i) = uv
            end do
    
            do i = k+1, sp
                VQ(i+r) = VP(i)
            end do

            UQ = UP

            ! Save the alphas
            do j = 1, r
                L = k-p+j
                do i = 0, p-j+s
                    alpha(i,j) = (uv-VP(L+i))/(VP(i+k+1)-VP(L+i))
                end do
            end do
            ! For each row do
            do row = 0, np
                ! Save unaltered control points
                do i = 0, k-p
                    Qw_net(row,i,:) = Pw_net(row,i,:)
                end do

                do i = k-s, mp
                    Qw_net(row,i+r,:) = Pw_net(row,i,:)
                end do

                ! Load auxiliary control points
                do i = 0, p-s
                    Rw(i,:) =  Pw_net(row, k-p+i,:)
                end do

                ! Insert the knor r times
                do j = 1, r
                    L = k-p+j
                    do i = 0, p-j-s
                        Rw(i,:) = alpha(i,j)*Rw(i+1,:) + (1.0-alpha(i,j))*Rw(i,:)
                    end do
                    Qw_net(row,L,:) = Rw(0,:)
                    Qw_net(row,k+r-j-s,:) = Rw(p-j-s,:)
                end do

                ! Load the remaining control points
                do i=L+1, k-s-1
                    Qw_net(row,i,:) = Rw(i-L,:)
                end do
            end do
        end if U_dir_refn
    end subroutine surface_knot_insertion

    subroutine surface_knot_refinement(p, UP, q, VP, Pw_net, X_array, dir, Ubar, Vbar, Qw_net)
        ! Algorithm A5.5 form the NURBS Book
        ! Refine surface vector
        ! ---------------------------------------------------------------------
        ! Input: n, p, U, m, q, V, Pw, X, r, dir
        ! Output: Ubar, Vbar, Qw
        ! ---------------------------------------------------------------------
        ! Variables:
        ! np+1 are the rows of the Pw_net control net
        ! mp+1 are the columns of the Pw_net control net
        ! nq+1 are the rows of the Qw_net control net
        ! mq+1 are the columns of the Qw_net control net
        ! UP and VP are the U and V knot vectors before the refinement
        ! Ubar and Vbar are the U and V knot vectors after the refinement
        ! Pw_net and Qw_net are the control nets before and after refinement
        ! s: knot multiplicity
        ! r: number of insertions
        use bspline_basis_functions, only: find_span

        integer, intent(in) :: p, q
        real, dimension(0:), intent(in) :: UP, VP, X_array
        real, dimension(0:,0:,0:), intent(in) :: Pw_net
        character(len=*), intent(in) :: dir
        
        real, dimension(:), allocatable, intent(out) :: Ubar, Vbar
        real, dimension(:,:,:), allocatable, intent(out) :: Qw_net
        
        integer :: r, s, a, b, num_ref_pts
        integer :: row, j, k, i, np, mp, rp, l, ind, sp
        ! real, dimension(:,:), allocatable :: alpha
        real :: alfa

        ! s = 0

        num_ref_pts = size(X_array)
        r = size(X_array) - 1

        dir_refn: if (dir == "U") then
            np = size(Pw_net,1) - 1
            mp = size(Pw_net,2) - 1

            rp = np + p + 1

            call find_span(rp, p, X_array(0), UP, a)
            call find_span(rp, p, X_array(num_ref_pts-1), UP, b)

            b = b+1

            allocate(Qw_net(0:np+1+r,0:mp,0:3))
            allocate(Ubar(0:size(UP) + r))

            !! Load u-vector as in A5.1
            do j = 0, a
                Ubar(j) = UP(j)
            end do
    
            do j = b+p, mp
                Ubar(j+r+1) = UP(j)
            end do

            Vbar = VP

            ! Save unaltered control points
            do row = 0, mp
                do k = 0, a-p
                    Qw_net(k,row,:) = Pw_net(k,row,:)
                end do

                do k = b-1, np
                    Qw_net(k+r+1,row,:) = Pw_net(k,row,:)
                end do
            end do

            i = b + p - 1
            k = b + p + r

            do j = r, 0, -1

                do 
                    if (X_array(j) <= UP(i)) then
                        if (i > a) then
                            Ubar(k) = UP(i)
                            do row = 0, mp
                                Qw_net(k-p-1,row,:) = Pw_net(i-p-1,row,:)
                            end do
                            k = k-1
                            i = i-1
                        else
                            exit
                        end if
                    else
                        exit
                    end if
                end do

                do row = 0, mp
                    Qw_net(k-p-1,row,:) = Qw_net(k-p,row,:)
                end do

                do l = 1, p
                    ind = k-p+l
                    alfa = Ubar(k+l) - X_array(j)
                    if (abs(alfa) < 1e-5) then
                        do row = 0, mp
                            Qw_net(ind-1,row,:) = Qw_net(ind,row,:)
                        end do
                    else
                        alfa = alfa/(Ubar(k+l) - UP(i-p+l))
                        do row = 0, mp
                            Qw_net(ind-1,row,:) = alfa*Qw_net(ind-1,row,:) + (1.0-alfa)*Qw_net(ind,row,:)
                        end do
                    end if
                end do
                Ubar(k) = X_array(j)
                k = k - 1
            end do
        else if (dir == "V") then
            np = size(Pw_net,1) - 1
            mp = size(Pw_net,2) - 1

            sp = np + q + 1

            call find_span(sp, q, X_array(0), VP, a)
            call find_span(sp, q, X_array(num_ref_pts-1), VP, b)

            b = b+1

            allocate(Qw_net(0:np,0:mp+1+r,0:3))
            allocate(Vbar(0:size(VP) + r))

            !! Load u-vector as in A5.1
            do j = 0, a
                Vbar(j) = VP(j)
            end do
    
            do j = b+p, np
                Vbar(j+r+1) = VP(j)
            end do

            Ubar = UP

            ! Save unaltered control points
            do row = 0, np
                do k = 0, a-q
                    Qw_net(row,k,:) = Pw_net(row,k,:)
                end do

                do k = b-1, np
                    Qw_net(row,k+r+1,:) = Pw_net(row,k,:)
                end do
            end do

            i = b + q - 1
            k = b + q + r

            do j = r, 0, -1

                do 
                    if (X_array(j) <= VP(i)) then
                        if (i > a) then
                            Vbar(k) = VP(i)
                            do row = 0, np
                                Qw_net(row,k-q-1,:) = Pw_net(row,i-q-1,:)
                            end do
                            k = k-1
                            i = i-1
                        else
                            exit
                        end if
                    else
                        exit
                    end if
                end do

                do row = 0, np
                    Qw_net(row,k-p-1,:) = Qw_net(row,k-p,:)
                end do

                do l = 1, q
                    ind = k-q+l
                    alfa = Vbar(k+l) - X_array(j)
                    if (abs(alfa) < 1e-5) then
                        do row = 0, np
                            Qw_net(row,ind-1,:) = Qw_net(row,ind,:)
                        end do
                    else
                        alfa = alfa/(Vbar(k+l) - VP(i-q+l))
                        do row = 0, mp
                            Qw_net(row,ind-1,:) = alfa*Qw_net(row,ind-1,:) + (1.0-alfa)*Qw_net(row,ind,:)
                        end do
                    end if
                end do
                Vbar(k) = X_array(j)
                k = k - 1
            end do
        end if dir_refn

        ! do i = 0, 3
        !     call print_matrix(Qw_net(:,:,i))
        ! end do
    end subroutine surface_knot_refinement

    subroutine surface_degree_elevation
    end subroutine surface_degree_elevation

    subroutine surface_h_refinement(p, q, P_pts, w_pts, UP, VP, dir, ph, Ph_pts, wh_pts, Ubar, Vbar)
        ! Input: 
        !   p: previous spline degree
        !   P_pts: previous control points
        !   w_pts: previous weights
        !   U_knot: previous knot vector
        ! Output: 
        !   ph: spline degree after h-refinement
        !   Ph_pts: control points after h-refinement
        !   wh_pts: weights after h-refinement
        !   Ubar: knot vector after h-refinement
        use nurbs_curve
        use nurbs_surface
        use utils
        integer, intent(in) :: p, q
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(0:), intent(in) :: UP, VP
        character(len=*), intent(in) :: dir
        integer, intent(out) :: ph
        real, dimension(:,:), allocatable, intent(out) :: Ph_pts
        real, dimension(:,:), allocatable, intent(out) :: wh_pts
        real, dimension(:), allocatable, intent(out) :: Ubar, Vbar

        real, dimension(:), allocatable :: Ured, Vred, X_array
        real, dimension(:,:), allocatable :: Pw_pts, Qw_pts
        integer :: np, mp, mred, nq, mq, nred, r, s, nu, nv
        real, dimension(:,:,:), allocatable :: Pw_net, Qw_net

        r = size(UP) - 1
        s = size(VP) - 1
        nu = r - p - 1
        nv = s - q - 1
        
        call weighted_control_points(P_pts, w_pts, Pw_pts)
        call create_control_net(nu, nv, Pw_pts, Pw_net)

        ! mp = size(UP) - 1
        ! mq = size(VP) - 1
        ! np = rp - p - 1
        ! nq = sq - q - 1

        mred = r-(2*p)
        nred = s-(2*q)

        allocate(Ured(0:mred))
        allocate(Vred(0:nred))
        allocate(X_array(0:mred-1))
        
        ! X_array = 0.5*(Ured(0:mred-1) + Ured(1:mred))
        
        if (dir == "U") then
            Ured = UP(p:r-p)
            X_array = 0.5*(Ured(0:mred-1) + Ured(1:mred))
            call surface_knot_refinement(p, UP, q, VP, Pw_net, X_array, dir, Ubar, Vbar, Qw_net)
        end if

        call create_control_list(Qw_net, Qw_pts)
        call geometric_control_points(Qw_pts, Ph_pts, wh_pts)
        ph = p
    end subroutine surface_h_refinement

    subroutine surface_p_refinement
    end subroutine surface_p_refinement

    subroutine surface_k_refinement
    end subroutine surface_k_refinement

    subroutine surface_spline_refinement
    end subroutine surface_spline_refinement

    subroutine assess_surface_refinement(spts_1, spts_2)
        use utils
        real, dimension(:,:), intent(in) :: spts_1, spts_2
        real, dimension(:,:), allocatable :: diff
        real, dimension(:), allocatable :: norm_diff
        logical, dimension(:), allocatable :: norm_diff_bool
        real :: tol = 1e-5
        real :: average_error
        
        diff = spts_1-spts_2

        norm_diff = norm2(diff,2)

        norm_diff_bool = norm_diff < tol

        if (all(norm_diff_bool)) then
            print *, "Successful refinement"
        else
            print *, "Failure on refinement"
        end if

        average_error = sum(norm_diff)/size(norm_diff)

        print *, "Average L2 error: "
        print *, average_error
    end subroutine assess_surface_refinement
end module surface_refinement