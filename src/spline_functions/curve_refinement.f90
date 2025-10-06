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
        real, dimension(0:), intent(in) :: U_knot
        real, dimension(0:,0:), intent(in) :: Pw_pts

        real, dimension(:), allocatable, intent(out) :: UQ
        real, dimension(:,:), allocatable, intent(out) :: Qw_pts
        
        integer :: np, mp, nq, i, L, j, r, s, k
        real :: alpha
        real, dimension(:,:), allocatable :: Rw_pts

        r = 1
        s = 0

        mp = size(U_knot) - 1
        np = mp - p - 1
        ! mp = np + p + 1
        nq = np + r

        call find_span(mp, p, u, U_knot, k)

        allocate(UQ(0:nq+p+1))
        allocate(Rw_pts(0:p, size(Pw_pts,2)))
        allocate(Qw_pts(0:size(Pw_pts,1), size(Pw_pts,2)))

        Rw_pts = 0.0
        Qw_pts = 0.0

        do i = 0, k
            UQ(i) = U_knot(i)
        end do

        do i = 1, r
            UQ(k+i) = u
        end do

        do i = k+1, mp
            UQ(i+r) = U_knot(i)
        end do

        ! Save unaltered control points
        do i = 0, k-p
            Qw_pts(i,:) = Pw_pts(i,:)
        end do

        do i = k-s, np
            Qw_pts(i+r,:) = Pw_pts(i,:)
        end do

        do i = 0, p-s
            Rw_pts(i,:) = Pw_pts(k-p+i,:)
        end do

        ! Insert the knot r times
        do j = 1, r
            L = k-p+j

            do i = 0, p-j-s
                alpha = (u - U_knot(L+i))/(U_knot(i+k+1)- U_knot(L+i))
                Rw_pts(i,:) = alpha*Rw_pts(i+1,:) + (1.0 - alpha)*Rw_pts(i,:)
            end do

            Qw_pts(L,:) = Rw_pts(0,:)
            Qw_pts(k+r-j-s,:) = Rw_pts(p-j-s,:)
        end do

        ! Load remaining control points
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
        
        real, dimension(0:), intent(in) :: U_knot, X_array
        real, dimension(0:,0:), intent(in) :: Pw_pts
        integer, intent(in) :: p

        real, dimension(:), allocatable, intent(out) :: Ubar
        real, dimension(:,:), allocatable, intent(out) :: Qw_pts

        integer :: a, b, num_ref_pts, m
        integer :: n, r, j, i, k, l, ind
        real :: alpha

        num_ref_pts = size(X_array)
        m = size(U_knot,1) - 1

        call find_span(m, p, X_array(0), U_knot, a)
        call find_span(m, p, X_array(num_ref_pts-1), U_knot, b)

        b = b+1

        n = size(Pw_pts,1) - 1
        r = size(X_array) - 1
        m = n + p + 1
        allocate(Qw_pts(0:n + 1 + r,size(Pw_pts,2)))
        allocate(Ubar(0:size(U_knot) + r))
        Qw_pts = 0.0
        Ubar = 0.0

        do j = 0, a-p+0
            Qw_pts(j,:) = Pw_pts(j,:)
        end do

        do j = b-1, n
            Qw_pts(j+r+1,:) = Pw_pts(j,:)
        end do

        do j = 0, a
            Ubar(j) = U_knot(j)
        end do

        do j = b+p, m
            Ubar(j+r+1) = U_knot(j)
        end do

        i = b + p - 1
        k = b + p + r
        do j = r,0,-1

            do 
                if (X_array(j) <= U_knot(i)) then
                    if (i > a) then
                        Qw_pts(k-p-1,:) = Pw_pts(i-p-1,:)
                        Ubar(k) = U_knot(i)
                        k = k - 1
                        i = i - 1
                    else
                        exit
                    end if
                else
                    exit
                end if
            end do

            Qw_pts(k-p-1,:) = Qw_pts(k-p,:)
            do l = 1, p
                ind = k - p + l
                alpha = Ubar(k+l) - X_array(j)
                if (abs(alpha) < 1e-5) then
                    Qw_pts(ind-1,:) = Qw_pts(ind,:)
                else
                    alpha = alpha/(Ubar(k+l) - U_knot(i-p+l))
                    Qw_pts(ind-1,:) = alpha*Qw_pts(ind-1,:) + (1.0 - alpha)*Qw_pts(ind,:)
                end if
            end do

            Ubar(k) = X_array(j)
            k = k - 1
        end do
    end subroutine knot_refinement

    subroutine pre_spline_decomposition
    end subroutine pre_spline_decomposition

    subroutine spline_splitting
    end subroutine spline_splitting

    subroutine spline_splitting_v2
    end subroutine spline_splitting_v2

    subroutine degree_elevation(p, U_knot, Pw, t, Qw, Uh)
        ! Algorithm 5.9 from the NURBS book
        ! Degree eleate a curve t times
        ! Input: n, p, U, Pw, t
        ! Output: nh, Uh, Qw
        use utils
        use nurbs_curve

        integer, intent(in) :: p
        real, dimension(0:), intent(in) :: U_knot
        real, dimension(0:,0:), intent(in) :: Pw
        integer, intent(in) :: t
        real, dimension(:,:), allocatable, intent(out) :: Qw
        real, dimension(:), allocatable, intent(out) :: Uh
        
        integer :: n, m, ph, ph2, mpi, i, j, k, s
        integer :: mh, kind, r, a, b, cind, oldr, lbz, rbz, save_var
        integer :: first, last, kj, tr, nh, mul
        real :: bc, bc1, bc2
        real :: inv, ua, ub, numer, den, bet, alf, gam
        real, dimension(:), allocatable :: alfs, Uh_prev
        real, dimension(:,:), allocatable :: bezalfs, bpts, ebpts, Nextbpts, Qw_prev

        n = size(Pw,1) - 1
        m = n + p + 1
        ph = p + t
        ph2 = ph/2

        allocate(bezalfs(0:p+t,0:p))
        allocate(bpts(0:p,size(Pw,2)))
        allocate(ebpts(0:p+t,size(Pw,2)))
        allocate(Nextbpts(0:p-2,size(Pw,2)))
        allocate(alfs(0:p-2))
        allocate(Uh_prev(0:2*size(U_knot)-1))
        allocate(Qw_prev(0:size(Uh_prev)-ph-1,0:size(Pw,2)-1))

        bpts = 0.0
        ebpts = 0.0
        Nextbpts = 0.0
        alfs = 0.0
        Uh_prev = 0.0
        Qw_prev = 0.0

        ! Compute bezier degree elevation coefficients
        bezalfs(0,0) = 1.0
        bezalfs(ph,p) = 1.0

        do i = 1, ph2
            call binomial(ph, i, bc)
            inv = 1.0/bc
            mpi = min(p,i)
            
            do j = max(0,i-t), mpi
                call binomial(p, j, bc1)
                call binomial(t, i-j, bc2)
                bezalfs(i,j) = inv*bc1*bc2
            end do
        end do

        do i = ph2+1, ph-1
            mpi = min(p,i)
            do j = max(0,i-t), mpi
                bezalfs(i,j) = bezalfs(ph-i,p-j)
            end do
        end do

        mh = ph
        kind = ph + 1
        r = -1
        a = p
        b = p + 1
        cind = 1
        ua = U_knot(0)
        Qw_prev(0,:) = Pw(0,:)

        do i = 0, ph
            Uh_prev(i) = ua
        end do

        ! Initialize first bezier segment
        do i = 0, p
            bpts(i,:) = Pw(i,:)
        end do

        ! Big loop through knot vector
        do while (b < m)
            i = b
            do 
                if ((b < m)) then
                    if (abs(U_knot(b+1) - U_knot(b)) < 1e-5) then
                        b = b + 1
                    else
                        exit
                    end if
                else
                    exit
                end if
            end do

            mul = b - i + 1
            mh = mh + mul + t
            ub = U_knot(b)
            oldr = r
            r = p - mul
            ! Insert knot u(b) r times
            if (oldr > 0) then
                lbz = (oldr + 2)/2
            else
                lbz = 1
            end if

            if (r > 0) then
                rbz = ph - (r + 1)/2
            else
                rbz = ph
            end if

            ! Insert knot to get bezier segment
            if (r > 0) then
                numer = ub - ua
                ! for k in range(p,mul,-1):
                do k = p, mul+1, -1
                    alfs(k-mul-1) = numer/(U_knot(a+k) - ua)
                end do

                do j = 1, r
                    save_var = r - j
                    s = mul + j
                    do k = p, s,-1
                        bpts(k,:) = alfs(k-s)*bpts(k,:) + (1.0 - alfs(k-s))*bpts(k-1,:)
                    end do

                    Nextbpts(save_var,:) = bpts(p,:)
                end do
            end if

            ! End of insert knot
            ! Degree elevate bezier
            ! Only points lbz,...,ph are used below
            do i = lbz, ph
                ebpts(i,:) = 0.0
                mpi = min(p,i)
                do j = max(0,i-t), mpi
                    ebpts(i,:) = ebpts(i,:) + bezalfs(i,j)*bpts(j,:)
                end do
            end do

            ! print('Bezier + t control points of current segment')
            ! print(ebpts)
            ! End of degree elevating bezier
            ! Must remove knot u = U(a) oldr times
            if (oldr > 1) then
                first = kind - 2
                last = kind
                den = ub - ua
                bet = (ub - Uh_prev(kind-1))/den
                ! Knot removal loop
                do tr = 1, oldr-1
                    i = first
                    j = last
                    kj = j - kind + 1
                    ! Loop and compute the new
                    ! Control points for one removal step
                    do while ((j - i) > tr)
                        if (i < cind) then
                            alf = (ub - Uh_prev(i))/(ua - Uh_prev(i))
                            Qw_prev(i,:) = alf*Qw_prev(i,:) + (1.0 - alf)*Qw_prev(i-1,:)
                        end if

                        if (j >= lbz) then
                            if ((j - tr) <= (kind - ph + oldr)) then
                                gam = (ub - Uh_prev(j - tr))/den
                                ebpts(kj,:) = gam*ebpts(kj,:) + (1.0 - gam)*ebpts(kj+1,:)
                            else
                                ebpts(kj,:) = bet*ebpts(kj,:) + (1.0 - bet)*ebpts(kj+1,:)
                            end if
                        end if
                        i = i + 1
                        j = j - 1
                        kj = kj - 1
                    end do

                    first = first - 1
                    last = last + 1
                end do
            end if

            ! End of removing knot u = U(a)
            ! Load the knot ua
            if (a /= p) then
                do i = 0, ph-oldr-1
                    Uh_prev(kind) = ua
                    kind = kind + 1
                end do
            end if

            ! Load control points into Qw
            do j = lbz, rbz
                Qw_prev(cind,:) = ebpts(j,:)
                cind = cind + 1
            end do

            if (b < m) then
                ! Set up for the next loop through
                do j = 0, r-1
                    bpts(j,:) = Nextbpts(j,:)
                end do

                do j = r, p
                    bpts(j,:) = Pw(b-p+j,:)
                end do

                a = b
                b = b + 1
                ua = ub
            else
                ! End knot
                do i = 0, ph
                    Uh_prev(kind+i) = ub
                end do
            end if

        end do
        ! End of while loop (b < m)
        nh = mh - ph - 1

        allocate(Uh(0:mh))
        allocate(Qw(0:nh,0:size(Pw,2)-1))
        Uh = 0.0
        
        Uh = Uh_prev(0:mh)
        Qw = Qw_prev(0:nh,:)
    end subroutine degree_elevation

    subroutine h_refinement(p, P_pts, w_pts, U_knot, ph, Ph_pts, wh_pts, Ubar)
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
        use utils
        integer, intent(in) :: p
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(0:), intent(in) :: U_knot
        integer, intent(out) :: ph
        real, dimension(:,:), allocatable, intent(out) :: Ph_pts
        real, dimension(:,:), allocatable, intent(out) :: wh_pts
        real, dimension(:), allocatable, intent(out) :: Ubar

        real, dimension(:), allocatable :: Ured, X_array
        real, dimension(:,:), allocatable :: Pw_pts, Qw_pts
        integer :: n, m, mred

        call weighted_control_points(P_pts, w_pts, Pw_pts)

        m = size(U_knot) - 1
        n = m - p - 1
        mred = m-(2*p)

        allocate(Ured(0:mred))
        allocate(X_array(0:mred-1))
        
        Ured = U_knot(p:m-p)

        X_array = 0.5*(Ured(0:mred-1) + Ured(1:mred))

        call knot_refinement(p, U_knot, X_array, Pw_pts, Qw_pts, Ubar)
        call geometric_control_points(Qw_pts, Ph_pts, wh_pts)

        ph = p
    end subroutine h_refinement

    subroutine p_refinement(p, P_pts, w_pts, U_knot, pp, Pp_pts, wp_pts, Up)
        ! Input: 
        !   p: previous spline degree
        !   P_pts: previous control points
        !   w_pts: previous weights
        !   U_knot: previous knot vector
        ! Output: 
        !   pp: spline degree after p-refinement
        !   Pp_pts: control points after p-refinement
        !   wp_pts: weights after p-refinement
        !   Up: knot vector after p-refinement

        ! Pw = rbs.weightedControlPoints(P,w)
        ! Pw = np.hstack((P,np.ones((P.shape[0],1))))
        ! Pw *= w
        use nurbs_curve
        use utils
        integer, intent(in) :: p
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(:), intent(in) :: U_knot
        integer, intent(out) :: pp
        real, dimension(:,:), allocatable, intent(out) :: Pp_pts
        real, dimension(:,:), allocatable, intent(out) :: wp_pts
        real, dimension(:), allocatable, intent(out) :: Up

        real, dimension(:,:), allocatable :: Pw_pts, Qw_pts
        integer :: t

        call weighted_control_points(P_pts, w_pts, Pw_pts)

        t = 1
        pp = p+t
        call degree_elevation(p, U_knot, Pw_pts, t, Qw_pts, Up)

        call geometric_control_points(Qw_pts, Pp_pts, wp_pts)
    end subroutine p_refinement

    subroutine k_refinement(p, P_pts, w_pts, U_knot, pk, Pk_pts, wk_pts, Uk)
        ! Input: 
        !   pp: previous spline degree
        !   Pp_pts: previous control points
        !   wp_pts: previous weights
        !   Up: previous knot vector
        ! Output: 
        !   pk: spline degree after refinement
        !   Pk_pts: control points after refinement
        !   wk_pts: weights after refinement
        !   Uk: knot vector after refinement
        integer, intent(in) :: p
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(:), intent(in) :: U_knot
        integer, intent(out) :: pk
        real, dimension(:,:), allocatable, intent(out) :: Pk_pts
        real, dimension(:,:), allocatable, intent(out) :: wk_pts
        real, dimension(:), allocatable, intent(out) :: Uk

        integer :: pp
        real, dimension(:,:), allocatable :: Pp_pts, wp_pts
        real, dimension(:), allocatable :: Up
        
        call p_refinement(p, P_pts, w_pts, U_knot, pp, Pp_pts, wp_pts, Up)
        call h_refinement(pp, Pp_pts, wp_pts, Up, pk, Pk_pts, wk_pts, Uk)
    end subroutine k_refinement

    subroutine spline_refinement(p, P_pts, w_pts, U_knot, ref_list, pref, Pref_pts, wref_pts, Uref)
        ! Input: 
        !   p: previous spline degree
        !   P_pts: previous control points
        !   w_pts: previous weights
        !   U: previous knot vector
        ! Output: 
        !   pref: spline degree after refinement
        !   Pref_pts: control points after refinement
        !   wref_pts: weights after refinement
        !   Uref: knot vector after refinement
        integer, intent(in) :: p
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(:), intent(in) :: U_knot
        character(len=1), dimension(:), intent(in) :: ref_list
        integer, intent(out) :: pref
        real, dimension(:,:), allocatable, intent(out) :: Pref_pts
        real, dimension(:,:), allocatable, intent(out) :: wref_pts
        real, dimension(:), allocatable, intent(out) :: Uref

        integer :: h_ref, p_ref, k_ref, i
        
        integer :: pprev
        real, dimension(:,:), allocatable :: Pprev_pts
        real, dimension(:,:), allocatable :: wprev_pts
        real, dimension(:), allocatable :: Uprev_knot

        h_ref = 0
        p_ref = 0
        k_ref = 0

        pprev = p
        Pprev_pts = P_pts
        wprev_pts = w_pts
        Uprev_knot = U_knot

        do i = 1, size(ref_list)
            ! print *, ref_list(i)
            if (ref_list(i) == 'h') then
                ! Uout,pout,Pout,wout = hRefinement(Uin,pin,Pin,win)
                call h_refinement(pprev, Pprev_pts, wprev_pts, Uprev_knot, pref, Pref_pts, wref_pts, Uref)
                h_ref = h_ref + 1
            else if (ref_list(i) == 'p') then
                ! Uout,pout,Pout,wout = pRefinement(Uin,pin,Pin,win)
                call p_refinement(pprev, Pprev_pts, wprev_pts, Uprev_knot, pref, Pref_pts, wref_pts, Uref)
                p_ref = p_ref + 1
            else if (ref_list(i) == 'k') then
                ! Uout,pout,Pout,wout = kRefinement(Uin,pin,Pin,win)
                call k_refinement(pprev, Pprev_pts, wprev_pts, Uprev_knot, pref, Pref_pts, wref_pts, Uref)
                k_ref = k_ref + 1
            else
                print *, "Invalid option"
            end if

            pprev = pref
            Pprev_pts = Pref_pts
            wprev_pts = wref_pts
            Uprev_knot = Uref
        end do

        print "(A, I3)", 'Number of h-refinements:', h_ref
        print "(A, I3)", 'Number of p-refinements:', p_ref
        print "(A, I3)", 'Number of k-refinements:', k_ref
    end subroutine spline_refinement

    subroutine assess_refinement(cpts_1, cpts_2)
        use utils
        real, dimension(:,:), intent(in) :: cpts_1, cpts_2
        real, dimension(:,:), allocatable :: diff
        real, dimension(:), allocatable :: norm_diff
        logical, dimension(:), allocatable :: norm_diff_bool
        real :: tol = 1e-5
        real :: average_error
        
        diff = cpts_1-cpts_2

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
    end subroutine assess_refinement

end module curve_refinement