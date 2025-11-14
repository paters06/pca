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

            sp = nq + q + 1

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
        
        integer :: r, a, b, num_ref_pts
        integer :: row, j, k, i, np, mp, rp, l, ind, sp
        ! real, dimension(:,:), allocatable :: alpha
        real :: alfa

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
    
            do j = b+p, rp
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

            sp = mp + q + 1

            call find_span(sp, q, X_array(0), VP, a)
            call find_span(sp, q, X_array(num_ref_pts-1), VP, b)

            b = b+1

            allocate(Qw_net(0:np,0:mp+1+r,0:3))
            allocate(Vbar(0:size(VP) + r))

            !! Load u-vector as in A5.1
            do j = 0, a
                Vbar(j) = VP(j)
            end do
    
            do j = b+q, sp
                Vbar(j+r+1) = VP(j)
            end do

            Ubar = UP

            ! Save unaltered control points
            do row = 0, np
                do k = 0, a-q
                    Qw_net(row,k,:) = Pw_net(row,k,:)
                end do

                do k = b-1, mp
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
                    Qw_net(row,k-q-1,:) = Qw_net(row,k-q,:)
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
                        do row = 0, np
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

    subroutine surface_degree_elevation(p, UP, q, VP, Pw_net, dir, t, ph, Uh, qh, Vh, Qw_net)
        ! Algorithm A5.10 from the NURBS book
        ! Degree elevate a surface t times
        ! ---------------------------------------------------------------------
        ! Input: n, p, U, m, q, V, Pw, dir, t
        ! Output: nh, Uh, mh, Vh, Qw
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
        use nurbs_curve_module
        integer, intent(in) :: p, q, t
        real, dimension(0:), intent(in) :: UP, VP
        real, dimension(0:,0:,0:), intent(in) :: Pw_net
        character(len=*), intent(in) :: dir
        integer, intent(out) :: ph, qh
        real, dimension(:), allocatable, intent(out) :: Uh, Vh
        real, dimension(:,:,:), allocatable, intent(out) :: Qw_net

        integer :: np, mp, rp, sq, ph2, qh2, i, mpi, j, kind, r, a, b, cind
        integer :: mul, oldr, lbz, rbz, k, save_var, s, first, last, tr, kj
        integer :: rh, sh, nh, mqi
        real :: bc, inv, bc1, bc2, ua, ub, numer, den, bet, alf, gam, va, vb
        real, dimension(:), allocatable :: alfs, Uh_prev, Vh_prev
        real, dimension(:,:), allocatable :: bezalfs
        real, dimension(:,:,:), allocatable :: bpts, ebpts, Nextbpts, Qw_prev

        dir_refn: if (dir == "U") then
            np = size(Pw_net,1) - 1
            mp = size(Pw_net,2) - 1
            rp = np + p + 1
            ph = p + t
            ph2 = ph/2

            allocate(bezalfs(0:p+t,0:p))
            allocate(bpts(0:p,0:mp,0:3))
            allocate(ebpts(0:p+t,0:mp,0:3))
            allocate(Nextbpts(0:p-2,size(Pw_net,2),0:3))
            allocate(alfs(0:p-2))
            allocate(Uh_prev(0:2*size(UP)-1))
            allocate(Qw_prev(0:size(Uh_prev)-ph-1,0:mp,0:3))

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

            rh = ph
            kind = ph + 1
            r = -1
            a = p
            b = p + 1
            cind = 1
            ua = UP(0)
            Qw_prev(0,:,:) = Pw_net(0,:,:)

            do i = 0, ph
                Uh_prev(i) = ua
            end do

            ! Initialize first bezier segment
            do i = 0, p
                bpts(i,:,:) = Pw_net(i,:,:)
            end do

            ! Big loop through knot vector
            do while (b < rp)
                i = b
                do 
                    if ((b < rp)) then
                        if (abs(UP(b+1) - UP(b)) < 1e-5) then
                            b = b + 1
                        else
                            exit
                        end if
                    else
                        exit
                    end if
                end do

                mul = b - i + 1
                rh = rh + mul + t
                ub = UP(b)
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
                        alfs(k-mul-1) = numer/(UP(a+k) - ua)
                    end do

                    do j = 1, r
                        save_var = r - j
                        s = mul + j
                        do k = p, s,-1
                            bpts(k,:,:) = alfs(k-s)*bpts(k,:,:) + (1.0 - alfs(k-s))*bpts(k-1,:,:)
                        end do

                        Nextbpts(save_var,:,:) = bpts(p,:,:)
                    end do
                end if

                ! End of insert knot
                ! Degree elevate bezier
                ! Only points lbz,...,ph are used below
                do i = lbz, ph
                    ebpts(i,:,:) = 0.0
                    mpi = min(p,i)
                    do j = max(0,i-t), mpi
                        ebpts(i,:,:) = ebpts(i,:,:) + bezalfs(i,j)*bpts(j,:,:)
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
                                Qw_prev(i,:,:) = alf*Qw_prev(i,:,:) + (1.0 - alf)*Qw_prev(i-1,:,:)
                            end if

                            if (j >= lbz) then
                                if ((j - tr) <= (kind - ph + oldr)) then
                                    gam = (ub - Uh_prev(j - tr))/den
                                    ebpts(kj,:,:) = gam*ebpts(kj,:,:) + (1.0 - gam)*ebpts(kj+1,:,:)
                                else
                                    ebpts(kj,:,:) = bet*ebpts(kj,:,:) + (1.0 - bet)*ebpts(kj+1,:,:)
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
                    Qw_prev(cind,:,:) = ebpts(j,:,:)
                    cind = cind + 1
                end do

                if (b < rp) then
                    ! Set up for the next loop through
                    do j = 0, r-1
                        bpts(j,:,:) = Nextbpts(j,:,:)
                    end do

                    do j = r, p
                        bpts(j,:,:) = Pw_net(b-p+j,:,:)
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
            nh = rh - ph - 1

            allocate(Uh(0:rh))
            allocate(Qw_net(0:nh,0:mp,0:3))
            Uh = 0.0
            
            Uh = Uh_prev(0:rh)
            Qw_net = Qw_prev(0:nh,:,:)

            Vh = VP
            qh = q
        else if (dir == "V") then
            np = size(Pw_net,1) - 1
            mp = size(Pw_net,2) - 1
            sq = mp + q + 1
            qh = q + t
            qh2 = qh/2

            ! Fix allocations
            allocate(bezalfs(0:q+t,0:q))
            allocate(bpts(0:np,0:q,0:3))
            allocate(ebpts(0:np,0:q+t,0:3))
            allocate(Nextbpts(0:np,0:q-2,0:3))
            allocate(alfs(0:q-2))
            allocate(Vh_prev(0:2*size(VP)-1))
            allocate(Qw_prev(0:np,0:size(Vh_prev)-qh-1,0:3))

            bpts = 0.0
            ebpts = 0.0
            Nextbpts = 0.0
            alfs = 0.0
            Vh_prev = 0.0
            Qw_prev = 0.0

            ! Compute bezier degree elevation coefficients
            bezalfs(0,0) = 1.0
            bezalfs(qh,q) = 1.0

            do i = 1, qh2
                call binomial(qh, i, bc)
                inv = 1.0/bc
                mqi = min(q,i)
                
                do j = max(0,i-t), mqi
                    call binomial(q, j, bc1)
                    call binomial(t, i-j, bc2)
                    bezalfs(i,j) = inv*bc1*bc2
                end do
            end do

            do i = qh2+1, qh-1
                mqi = min(q,i)
                do j = max(0,i-t), mqi
                    bezalfs(i,j) = bezalfs(qh-i,q-j)
                end do
            end do

            sh = qh
            kind = qh + 1
            r = -1
            a = q
            b = q + 1
            cind = 1
            va = VP(0)
            Qw_prev(:,0,:) = Pw_net(:,0,:)

            do i = 0, qh
                Vh_prev(i) = va
            end do

            ! Initialize first bezier segment
            do i = 0, q
                bpts(:,i,:) = Pw_net(:,i,:)
            end do

            ! Big loop through knot vector
            do while (b < sq)
                i = b
                do 
                    if ((b < sq)) then
                        if (abs(VP(b+1) - VP(b)) < 1e-5) then
                            b = b + 1
                        else
                            exit
                        end if
                    else
                        exit
                    end if
                end do

                mul = b - i + 1
                sh = sh + mul + t
                vb = VP(b)
                oldr = r
                r = q - mul
                
                ! Insert knot u(b) r times
                if (oldr > 0) then
                    lbz = (oldr + 2)/2
                else
                    lbz = 1
                end if

                if (r > 0) then
                    rbz = qh - (r + 1)/2
                else
                    rbz = qh
                end if

                ! Insert knot to get bezier segment
                if (r > 0) then
                    numer = vb - va
                    ! for k in range(p,mul,-1):
                    do k = q, mul+1, -1
                        alfs(k-mul-1) = numer/(VP(a+k) - va)
                    end do

                    do j = 1, r
                        save_var = r - j
                        s = mul + j
                        do k = q, s, -1
                            bpts(:,k,:) = alfs(k-s)*bpts(:,k,:) + (1.0 - alfs(k-s))*bpts(:,k-1,:)
                        end do

                        Nextbpts(:,save_var,:) = bpts(:,q,:)
                    end do
                end if

                ! End of insert knot
                ! Degree elevate bezier
                ! Only points lbz,...,qh are used below
                do i = lbz, qh
                    ebpts(:,i,:) = 0.0
                    mqi = min(q,i)
                    do j = max(0,i-t), mqi
                        ebpts(:,i,:) = ebpts(:,i,:) + bezalfs(i,j)*bpts(:,j,:)
                    end do
                end do

                ! print('Bezier + t control points of current segment')
                ! print(ebpts)
                ! End of degree elevating bezier
                ! Must remove knot v = V(a) oldr times
                if (oldr > 1) then
                    first = kind - 2
                    last = kind
                    den = vb - va
                    bet = (vb - Vh_prev(kind-1))/den
                    ! Knot removal loop
                    do tr = 1, oldr-1
                        i = first
                        j = last
                        kj = j - kind + 1
                        ! Loop and compute the new
                        ! Control points for one removal step
                        do while ((j - i) > tr)
                            if (i < cind) then
                                alf = (vb - Vh_prev(i))/(va - Vh_prev(i))
                                Qw_prev(:,i,:) = alf*Qw_prev(:,i,:) + (1.0 - alf)*Qw_prev(:,i-1,:)
                            end if

                            if (j >= lbz) then
                                if ((j - tr) <= (kind - qh + oldr)) then
                                    gam = (vb - Vh_prev(j - tr))/den
                                    ebpts(:,kj,:) = gam*ebpts(:,kj,:) + (1.0 - gam)*ebpts(:,kj+1,:)
                                else
                                    ebpts(:,kj,:) = bet*ebpts(:,kj,:) + (1.0 - bet)*ebpts(:,kj+1,:)
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

                ! End of removing knot v = V(a)
                ! Load the knot ua
                if (a /= q) then
                    do i = 0, qh-oldr-1
                        Vh_prev(kind) = va
                        kind = kind + 1
                    end do
                end if

                ! Load control points into Qw
                do j = lbz, rbz
                    Qw_prev(:,cind,:) = ebpts(:,j,:)
                    cind = cind + 1
                end do

                if (b < sq) then
                    ! Set up for the next loop through
                    do j = 0, r-1
                        bpts(:,j,:) = Nextbpts(:,j,:)
                    end do

                    do j = r, q
                        bpts(:,j,:) = Pw_net(:,b-q+j,:)
                    end do

                    a = b
                    b = b + 1
                    va = vb
                else
                    ! End knot
                    do i = 0, qh
                        Vh_prev(kind+i) = vb
                    end do
                end if

            end do
            ! End of while loop (b < m)
            nh = sh - qh - 1

            allocate(Vh(0:sh))
            allocate(Qw_net(0:np,0:nh,0:3))
            Vh = 0.0
            
            Vh = Vh_prev(0:sh)
            Qw_net = Qw_prev(:,0:nh,:)

            Uh = UP
            ph = p
        end if dir_refn
    end subroutine surface_degree_elevation

    subroutine surface_h_refinement(p, q, P_pts, w_pts, UP, VP, ref_list, ph, qh, Ph_pts, wh_pts, Ubar, Vbar)
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
        use nurbs_curve_module
        use nurbs_surface_module
        use utils
        integer, intent(in) :: p, q
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(0:), intent(in) :: UP, VP
        ! character(len=*), intent(in) :: dir
        character(len=*), dimension(:,:), allocatable, intent(in) :: ref_list
        integer, intent(out) :: ph, qh
        real, dimension(:,:), allocatable, intent(out) :: Ph_pts
        real, dimension(:,:), allocatable, intent(out) :: wh_pts
        real, dimension(:), allocatable, intent(out) :: Ubar, Vbar

        real, dimension(:), allocatable :: X_array, Utemp, Vtemp
        real, dimension(:,:), allocatable :: Pw_pts, Qw_pts
        integer :: r, s, nu, nv, i
        real, dimension(:,:,:), allocatable :: Pw_net, Qw_net, Pw_temp
        character(:), allocatable :: dir

        r = size(UP) - 1
        s = size(VP) - 1
        nu = r - p - 1
        nv = s - q - 1
        
        call weighted_control_points(P_pts, w_pts, Pw_pts)
        call create_control_net(nu, nv, Pw_pts, Pw_net)

        Utemp = UP
        Vtemp = VP
        Pw_temp = Pw_net

        do i = 1, size(ref_list,1)
            dir = ref_list(i,2)
            if (dir == "U") then
                call compute_element_midvalues(Utemp, p, X_array)
            else if (dir == "V") then
                call compute_element_midvalues(Vtemp, q, X_array)
            end if
            
            call surface_knot_refinement(p, Utemp, q, Vtemp, Pw_temp, X_array, dir, Ubar, Vbar, Qw_net)

            Utemp = Ubar
            Vtemp = Vbar
            Pw_temp = Qw_net
        end do

        call create_control_list(Qw_net, Qw_pts)
        call geometric_control_points(Qw_pts, Ph_pts, wh_pts)
        ph = p
        qh = q
    end subroutine surface_h_refinement

    subroutine surface_p_refinement(p, q, P_pts, w_pts, UP, VP, ref_list, ph, qh, Ph_pts, wh_pts, Ubar, Vbar)
        ! Input: 
        !   p: previous spline degree
        !   q: previous spline degree
        !   P_pts: previous control points
        !   w_pts: previous weights
        !   U_knot: previous knot vector
        !   V_knot: previous knot vector
        ! Output: 
        !   ph: spline degree after p-refinement
        !   qh: spline degree after p-refinement
        !   Ph_pts: control points after p-refinement
        !   wh_pts: weights after p-refinement
        !   Ubar: knot vector after p-refinement
        !   Vbar: knot vector after p-refinement
        use nurbs_curve_module
        use nurbs_surface_module
        use utils
        integer, intent(in) :: p, q
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(0:), intent(in) :: UP, VP
        ! character(len=*), intent(in) :: dir
        character(len=*), dimension(:,:), allocatable, intent(in) :: ref_list
        integer, intent(out) :: ph, qh
        real, dimension(:,:), allocatable, intent(out) :: Ph_pts
        real, dimension(:,:), allocatable, intent(out) :: wh_pts
        real, dimension(:), allocatable, intent(out) :: Ubar, Vbar

        real, dimension(:), allocatable :: Utemp, Vtemp
        real, dimension(:,:), allocatable :: Pw_pts, Qw_pts
        integer :: r, s, nu, nv, i, t, ptemp, qtemp
        real, dimension(:,:,:), allocatable :: Pw_net, Qw_net, Pw_temp
        character(:), allocatable :: dir

        t = 1
        
        r = size(UP) - 1
        s = size(VP) - 1
        nu = r - p - 1
        nv = s - q - 1
        
        call weighted_control_points(P_pts, w_pts, Pw_pts)
        call create_control_net(nu, nv, Pw_pts, Pw_net)

        ptemp = p
        qtemp = q
        Utemp = UP
        Vtemp = VP
        Pw_temp = Pw_net

        do i = 1, size(ref_list,1)
            dir = ref_list(i,2)
            
            call surface_degree_elevation(ptemp, Utemp, qtemp, Vtemp, Pw_temp, dir, t, ph, Ubar, qh, Vbar, Qw_net)

            ptemp = ph
            qtemp = qh
            Utemp = Ubar
            Vtemp = Vbar
            Pw_temp = Qw_net
        end do

        call create_control_list(Qw_net, Qw_pts)
        call geometric_control_points(Qw_pts, Ph_pts, wh_pts)
    end subroutine surface_p_refinement

    subroutine surface_k_refinement(p, q, P_pts, w_pts, Uknot, Vknot, ref_list, pk, qk, Pk_pts, wk_pts, Uk, Vk)
        ! Input: 
        !   p: previous spline degree
        !   q: previous spline degree
        !   P_pts: previous control points
        !   w_pts: previous weights
        !   U_knot: previous knot vector
        !   V_knot: previous knot vector
        ! Output: 
        !   pp: spline degree after p-refinement
        !   qp: spline degree after p-refinement
        !   Pp_pts: control points after p-refinement
        !   wp_pts: weights after p-refinement
        !   Up: knot vector after p-refinement
        !   Vp: knot vector after p-refinement
        !   pk: spline degree after k-refinement
        !   qk: spline degree after k-refinement
        !   Pk_pts: control points after k-refinement
        !   wk_pts: weights after k-refinement
        !   Uk: knot vector after k-refinement
        !   Vk: knot vector after k-refinement
        use nurbs_curve_module
        use nurbs_surface_module
        use utils
        integer, intent(in) :: p, q
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(0:), intent(in) :: Uknot, Vknot
        character(len=*), dimension(:,:), allocatable, intent(in) :: ref_list
        integer, intent(out) :: pk, qk
        real, dimension(:,:), allocatable, intent(out) :: Pk_pts
        real, dimension(:,:), allocatable, intent(out) :: wk_pts
        real, dimension(:), allocatable, intent(out) :: Uk, Vk

        integer :: i, ptemp, qtemp, t, r, s, nu, nv, ph, qh
        character(:), allocatable :: dir

        real, dimension(:), allocatable :: Utemp, Vtemp, X_array, Uh, Vh
        real, dimension(:,:), allocatable :: Pw_pts, Qw_pts
        real, dimension(:,:,:), allocatable :: Pw_net, Qkw_net, Pw_temp, Qhw_net

        t = 1
        
        r = size(Uknot) - 1
        s = size(Vknot) - 1
        nu = r - p - 1
        nv = s - q - 1

        call weighted_control_points(P_pts, w_pts, Pw_pts)
        call create_control_net(nu, nv, Pw_pts, Pw_net)

        ptemp = p
        qtemp = q
        Utemp = Uknot
        Vtemp = Vknot
        Pw_temp = Pw_net

        do i = 1, size(ref_list,1)
            dir = ref_list(i,2)
            
            call surface_degree_elevation(ptemp, Utemp, qtemp, Vtemp, Pw_temp, dir, t, ph, Uh, qh, Vh, Qhw_net)
            
            if (dir == "U") then
                call compute_element_midvalues(Uh, ph, X_array)
            else if (dir == "V") then
                call compute_element_midvalues(Vh, qh, X_array)
            end if
            
            call surface_knot_refinement(ph, Uh, qh, Vh, Qhw_net, X_array, dir, Uk, Vk, Qkw_net)

            ptemp = ph
            qtemp = qh
            Utemp = Uk
            Vtemp = Vk
            Pw_temp = Qkw_net
        end do

        pk = ptemp
        qk = qtemp

        call create_control_list(Qkw_net, Qw_pts)
        call geometric_control_points(Qw_pts, Pk_pts, wk_pts)
    end subroutine surface_k_refinement

    subroutine surface_spline_refinement(p, q, P_pts, w_pts, Uknot, Vknot, ref_list, pk, qk, Pk_pts, wk_pts, Uk, Vk)
        use nurbs_curve_module
        use nurbs_surface_module
        use utils
        integer, intent(in) :: p, q
        real, dimension(:,:), intent(in) :: P_pts
        real, dimension(:,:), intent(in) :: w_pts
        real, dimension(0:), intent(in) :: Uknot, Vknot
        character(len=*), dimension(:,:), allocatable, intent(in) :: ref_list
        integer, intent(out) :: pk, qk
        real, dimension(:,:), allocatable, intent(out) :: Pk_pts
        real, dimension(:,:), allocatable, intent(out) :: wk_pts
        real, dimension(:), allocatable, intent(out) :: Uk, Vk

        integer :: i, ptemp, qtemp, t, r, s, nu, nv, ph, qh
        character(:), allocatable :: dir, refn_type

        real, dimension(:), allocatable :: Utemp, Vtemp, X_array, Uh, Vh
        real, dimension(:,:), allocatable :: Pw_pts, Qw_pts
        real, dimension(:,:,:), allocatable :: Pw_net, Qkw_net, Pw_temp, Qhw_net

        t = 1
        
        r = size(Uknot) - 1
        s = size(Vknot) - 1
        nu = r - p - 1
        nv = s - q - 1

        call weighted_control_points(P_pts, w_pts, Pw_pts)
        call create_control_net(nu, nv, Pw_pts, Pw_net)

        ptemp = p
        qtemp = q
        Utemp = Uknot
        Vtemp = Vknot
        Pw_temp = Pw_net

        do i = 1, size(ref_list,1)
            refn_type = ref_list(i,1)
            dir = ref_list(i,2)
            
            if (refn_type == "h") then
                
                if (dir == "U") then
                    call compute_element_midvalues(Utemp, ptemp, X_array)
                else if (dir == "V") then
                    call compute_element_midvalues(Vtemp, qtemp, X_array)
                end if
                
                call surface_knot_refinement(ptemp, Utemp, qtemp, Vtemp, Pw_temp, X_array, dir, Uk, Vk, Qkw_net)
                ph = ptemp
                qh = qtemp

            else if (refn_type == "p") then
                
                call surface_degree_elevation(ptemp, Utemp, qtemp, Vtemp, Pw_temp, dir, t, ph, Uk, qh, Vk, Qkw_net)

            else if (refn_type == "k") then
                
                call surface_degree_elevation(ptemp, Utemp, qtemp, Vtemp, Pw_temp, dir, t, ph, Uh, qh, Vh, Qhw_net)
            
                if (dir == "U") then
                    call compute_element_midvalues(Uh, ph, X_array)
                else if (dir == "V") then
                    call compute_element_midvalues(Vh, qh, X_array)
                end if
                
                call surface_knot_refinement(ph, Uh, qh, Vh, Qhw_net, X_array, dir, Uk, Vk, Qkw_net)

            else
                exit
            end if

            ptemp = ph
            qtemp = qh
            Utemp = Uk
            Vtemp = Vk
            Pw_temp = Qkw_net
        end do

        pk = ptemp
        qk = qtemp
        Uk = Utemp
        Vk = Vtemp

        call create_control_list(Qkw_net, Qw_pts)
        call geometric_control_points(Qw_pts, Pk_pts, wk_pts)
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