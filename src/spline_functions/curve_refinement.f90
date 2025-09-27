module curve_refinement
    implicit none
contains
    subroutine knot_insertion(p, U_knot, Pw_pts, u, UQ, Qw_pts)
        ! Algorithm 5.1 from the NURBS Book
        ! Input: np, p, UP, Pw, u, k, s, r
        ! Output: nq, UQ, Qw_pts
        ! s: knot multiplicity
        ! r: number of insertions
        use bspline_basis_functions
        use utils
        
        integer, intent(in) :: p
        real, intent(in) :: u
        real, dimension(:), intent(in) :: U_knot
        real, dimension(:,:), intent(in) :: Pw_pts

        real, dimension(:), allocatable, intent(out) :: UQ
        real, dimension(:,:), allocatable, intent(out) :: Qw_pts
        
        integer :: np, mp, nq, i, L, j, r, s, k
        real :: alpha
        real, dimension(:,:), allocatable :: Rw_pts

        r = 1
        s = 0

        np = (size(U_knot) - 1) - p - 1
        mp = np + p + 1
        nq = np + r

        call find_span(mp, p, u, U_knot, k)

        allocate(UQ(size(U_knot) + 1))
        allocate(Rw_pts(p+1, size(Pw_pts,2)))
        allocate(Qw_pts(size(Pw_pts,1)+1, size(Pw_pts,2)))

        Rw_pts = 0.0
        Qw_pts = 0.0

        do i = 1, k+1
            UQ(i) = U_knot(i)
        end do

        do i = 2, r+1
            UQ(k+i) = u
        end do

        do i = k+2, mp+1
            UQ(i+r) = U_knot(i)
        end do

        ! Save unaltered control points
        do i = 1, k-p+1
            Qw_pts(i,:) = Pw_pts(i,:)
        end do

        do i = k-s+1, np+1
            Qw_pts(i+r,:) = Pw_pts(i,:)
        end do

        do i = 1, p-s+1
            Rw_pts(i,:) = Pw_pts(k-p+i,:)
        end do

        ! Insert the knot r times
        do j = 1, r
            L = k-p+j

            ! indices plus one:
            do i = 1, p-j-s+1
                alpha = (u - U_knot(L+i))/(U_knot(i+k+1)- U_knot(L+i))
                Rw_pts(i,:) = alpha*Rw_pts(i+1,:) + (1.0 - alpha)*Rw_pts(i,:)
            end do

            ! indices plus one
            Qw_pts(L+1,:) = Rw_pts(1,:)
            Qw_pts(k+r-j-s+1,:) = Rw_pts(p-j-s+1,:)
        end do

        do i = L+1, k-s-1
            Qw_pts(i,:) = Rw_pts(i-L,:)
        end do
    end subroutine knot_insertion

    subroutine knot_refinement(p, U_knot, X_array, Pw_pts, Qw_pts, Ubar)
        ! Algorithm 5.4 from the NURBS Book
        ! Refine curve knot vector
        ! Input: n, p, U_knot, Pw_pts, X_array, r
        ! Output: Ubar, Qw_pts
        use bspline_basis_functions
        use utils
        
        real, dimension(:), intent(in) :: U_knot, X_array
        real, dimension(:,:), intent(in) :: Pw_pts
        integer, intent(in) :: p

        real, dimension(:), allocatable, intent(out) :: Ubar
        real, dimension(:,:), allocatable, intent(out) :: Qw_pts

        integer :: a, b, num_ref_pts, mi
        integer :: n, r, m, j, i, k, l, ind
        real :: alpha

        num_ref_pts = size(X_array)
        mi = size(U_knot,1)

        call find_span(mi, p, X_array(1), U_knot, a)
        call find_span(mi, p, X_array(num_ref_pts), U_knot, b)

        print *, a
        print *, b

        b = b+1

        n = size(Pw_pts,1) - 1
        r = size(X_array) - 1
        m = n + p + 1
        allocate(Qw_pts(n + 1 + r + 1,size(Pw_pts,2)))
        allocate(Ubar(size(U_knot) + r + 1))
        Qw_pts = 0.0
        Ubar = 0.0

        do j = 1, a-p+1
            Qw_pts(j,:) = Pw_pts(j,:)
        end do

        do j = b, n+1
            Qw_pts(j+r+1,:) = Pw_pts(j,:)
        end do

        do j = 1, a+1
            Ubar(j) = U_knot(j)
        end do

        do j = b+p+1, m+1
            Ubar(j+r+1) = U_knot(j)
        end do

        i = b + p - 1
        k = b + p + r
        do j = r+1,1,-1
            do while (X_array(j) <= U_knot(i) .and. i > a)
                Qw_pts(k-p-1,:) = Pw_pts(i-p-1,:)
                Ubar(k) = U_knot(i)
                k = k - 1
                i = i - 1
            end do

            Qw_pts(k-p-1,:) = Qw_pts(k-p,:)
            do l = 1, p
                ind = k - p + l
                alpha = Ubar(k+l) - X_array(j)
                if (abs(alpha) < 1e-5) then
                    Qw_pts(ind-1,:) = Qw_pts(ind,:)
                else
                    alpha = alpha/(Ubar(k+l) - U_knot(i-p+l))
                    ! print *, k+l
                    ! print *, i-p+l
                    Qw_pts(ind-1,:) = alpha*Qw_pts(ind-1,:) + (1.0 - alpha)*Qw_pts(ind,:)
                end if
            end do

            Ubar(k) = X_array(j)
            k = k - 1
        end do

        call print_row_vector(U_knot)
        call print_row_vector(Ubar)
    end subroutine knot_refinement

    subroutine pre_spline_decomposition
    end subroutine pre_spline_decomposition

    subroutine spline_splitting
    end subroutine spline_splitting

    subroutine spline_splitting_v2
    end subroutine spline_splitting_v2

    subroutine degree_elevation
    end subroutine degree_elevation

    subroutine h_refinement
    end subroutine h_refinement

    subroutine p_refinement
    end subroutine p_refinement

    subroutine k_refinement
    end subroutine k_refinement

    subroutine spline_refinement
    end subroutine spline_refinement

    subroutine assess_refinement(cpts_1, cpts_2)
        use utils
        real, dimension(:,:), intent(in) :: cpts_1, cpts_2
        real, dimension(:,:), allocatable :: diff
        real, dimension(:), allocatable :: norm_diff
        logical, dimension(:), allocatable :: norm_diff_bool
        real :: tol = 1e-5
        
        diff = cpts_1-cpts_2

        norm_diff = norm2(diff,2)

        norm_diff_bool = norm_diff < tol

        if (all(norm_diff_bool)) then
            print *, "Successful refinement"
        else
            print *, "Failure on refinement"
        end if
    end subroutine assess_refinement

end module curve_refinement