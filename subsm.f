c> \file subsm.f

      subroutine subsm ( n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy,
     +                   theta, xx, gg,
     +                   col, head, iword, wv, wn, iprint, info )
      implicit none
      integer          n, m, nsub, col, head, iword, iprint, info, 
     +                 ind(nsub), nbd(n)
      double precision theta, 
     +                 l(n), u(n), x(n), d(n), xp(n), xx(n), gg(n),
     +                 ws(n, m), wy(n, m), 
     +                 wv(2*m), wn(2*m, 2*m)

c     **********************************************************************
c
c     This routine contains the major changes in the updated version.
c     The changes are described in the accompanying paper
c
c      Jose Luis Morales, Jorge Nocedal
c      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large-Scale
c       Bound Constrained Optimization". Decemmber 27, 2010.
c
c             J.L. Morales  Departamento de Matematicas, 
c                           Instituto Tecnologico Autonomo de Mexico
c                           Mexico D.F.
c
c             J, Nocedal    Department of Electrical Engineering and
c                           Computer Science.
c                           Northwestern University. Evanston, IL. USA
c
c                           January 17, 2011
c
c      **********************************************************************
c                           
c
c     Subroutine subsm
c
c     Given xcp, l, u, r, an index set that specifies
c       the active set at xcp, and an l-BFGS matrix B 
c       (in terms of WY, WS, SY, WT, head, col, and theta), 
c       this subroutine computes an approximate solution
c       of the subspace problem
c
c       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
c
c             subject to l<=x<=u
c                       x_i=xcp_i for all i in A(xcp)
c                     
c       along the subspace unconstrained Newton direction 
c       
c          d = -(Z'BZ)^(-1) r.
c
c       The formula for the Newton direction, given the L-BFGS matrix
c       and the Sherman-Morrison formula, is
c
c          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
c 
c       where
c                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c
c     Note that this procedure for computing d differs 
c     from that described in [1]. One can show that the matrix K is
c     equal to the matrix M^[-1]N in that paper.
c
c     n is an integer variable.
c       On entry n is the dimension of the problem.
c       On exit n is unchanged.
c
c     m is an integer variable.
c       On entry m is the maximum number of variable metric corrections
c         used to define the limited memory matrix.
c       On exit m is unchanged.
c
c     nsub is an integer variable.
c       On entry nsub is the number of free variables.
c       On exit nsub is unchanged.
c
c     ind is an integer array of dimension nsub.
c       On entry ind specifies the coordinate indices of free variables.
c       On exit ind is unchanged.
c
c     l is a double precision array of dimension n.
c       On entry l is the lower bound of x.
c       On exit l is unchanged.
c
c     u is a double precision array of dimension n.
c       On entry u is the upper bound of x.
c       On exit u is unchanged.
c
c     nbd is a integer array of dimension n.
c       On entry nbd represents the type of bounds imposed on the
c         variables, and must be specified as follows:
c         nbd(i)=0 if x(i) is unbounded,
c                1 if x(i) has only a lower bound,
c                2 if x(i) has both lower and upper bounds, and
c                3 if x(i) has only an upper bound.
c       On exit nbd is unchanged.
c
c     x is a double precision array of dimension n.
c       On entry x specifies the Cauchy point xcp. 
c       On exit x(i) is the minimizer of Q over the subspace of
c                                                        free variables. 
c
c     d is a double precision array of dimension n.
c       On entry d is the reduced gradient of Q at xcp.
c       On exit d is the Newton direction of Q. 
c
c    xp is a double precision array of dimension n.
c       used to safeguard the projected Newton direction 
c
c    xx is a double precision array of dimension n
c       On entry it holds the current iterate
c       On output it is unchanged

c    gg is a double precision array of dimension n
c       On entry it holds the gradient at the current iterate
c       On output it is unchanged
c
c     ws and wy are double precision arrays;
c     theta is a double precision variable;
c     col is an integer variable;
c     head is an integer variable.
c       On entry they store the information defining the
c                                          limited memory BFGS matrix:
c         ws(n,m) stores S, a set of s-vectors;
c         wy(n,m) stores Y, a set of y-vectors;
c         theta is the scaling factor specifying B_0 = theta I;
c         col is the number of variable metric corrections stored;
c         head is the location of the 1st s- (or y-) vector in S (or Y).
c       On exit they are unchanged.
c
c     iword is an integer variable.
c       On entry iword is unspecified.
c       On exit iword specifies the status of the subspace solution.
c         iword = 0 if the solution is in the box,
c                 1 if some bound is encountered.
c
c     wv is a double precision working array of dimension 2m.
c
c     wn is a double precision array of dimension 2m x 2m.
c       On entry the upper triangle of wn stores the LEL^T factorization
c         of the indefinite matrix
c
c              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                  [L_a -R_z           theta*S'AA'S ]
c                                                    where E = [-I  0]
c                                                              [ 0  I]
c       On exit wn is unchanged.
c
c     iprint is an INTEGER variable that must be set by the user.
c       It controls the frequency and type of output generated:
c        iprint<0    no output is generated;
c        iprint=0    print only one line at the last iteration;
c        0<iprint<99 print also f and |proj g| every iprint iterations;
c        iprint=99   print details of every iteration except n-vectors;
c        iprint=100  print also the changes of active set and final x;
c        iprint>100  print details of every iteration including x and g;
c       When iprint > 0, the file iterate.dat will be created to
c                        summarize the iteration.
c
c     info is an integer variable.
c       On entry info is unspecified.
c       On exit info = 0       for normal return,
c                    = nonzero for abnormal return 
c                                  when the matrix K is ill-conditioned.
c
c     Subprograms called:
c
c       Linpack dtrsl.
c
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c
c
c                           *  *  *
c
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
c
c
c     ************

      integer          pointr,m2,col2,ibd,jy,js,i,j,k
      double precision alpha, xk, dk, temp1, temp2 
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
c
      double precision dd_p

      if (nsub .le. 0) return
      if (iprint .ge. 99) write (6,1001)

c     Compute wv = W'Zd.

      pointr = head 
      do 20 i = 1, col
         temp1 = zero
         temp2 = zero
         do 10 j = 1, nsub
            k = ind(j)
            temp1 = temp1 + wy(k,pointr)*d(j)
            temp2 = temp2 + ws(k,pointr)*d(j)
  10     continue
         wv(i) = temp1
         wv(col + i) = theta*temp2
         pointr = mod(pointr,m) + 1
  20  continue
 
c     Compute wv:=K^(-1)wv.

      m2 = 2*m
      col2 = 2*col
      call dtrsl(wn,m2,col2,wv,11,info)
      if (info .ne. 0) return
      do 25 i = 1, col
         wv(i) = -wv(i)
  25     continue
      call dtrsl(wn,m2,col2,wv,01,info)
      if (info .ne. 0) return
 
c     Compute d = (1/theta)d + (1/theta**2)Z'W wv.
 
      pointr = head
      do 40 jy = 1, col
         js = col + jy
         do 30 i = 1, nsub
            k = ind(i)
            d(i) = d(i) + wy(k,pointr)*wv(jy)/theta     
     +                  + ws(k,pointr)*wv(js)
  30     continue
         pointr = mod(pointr,m) + 1
  40  continue

      call dscal( nsub, one/theta, d, 1 )
c 
c-----------------------------------------------------------------
c     Let us try the projection, d is the Newton direction

      iword = 0

      call dcopy ( n, x, 1, xp, 1 )
c
      do 50 i=1, nsub
         k  = ind(i)
         dk = d(i)
         xk = x(k)
         if ( nbd(k) .ne. 0 ) then
c
            if ( nbd(k).eq.1 ) then          ! lower bounds only
               x(k) = max( l(k), xk + dk )
               if ( x(k).eq.l(k) ) iword = 1
            else 
c     
               if ( nbd(k).eq.2 ) then       ! upper and lower bounds
                  xk   = max( l(k), xk + dk ) 
                  x(k) = min( u(k), xk )
                  if ( x(k).eq.l(k) .or. x(k).eq.u(k) ) iword = 1
               else
c
                  if ( nbd(k).eq.3 ) then    ! upper bounds only
                     x(k) = min( u(k), xk + dk )
                     if ( x(k).eq.u(k) ) iword = 1
                  end if 
               end if
            end if
c            
         else                                ! free variables
            x(k) = xk + dk
         end if 
 50   continue
c
      if ( iword.eq.0 ) then
         go to 911
      end if
c
c     check sign of the directional derivative
c
      dd_p = zero
      do 55 i=1, n
         dd_p  = dd_p + (x(i) - xx(i))*gg(i)
 55   continue
      if ( dd_p .gt.zero ) then
         call dcopy( n, xp, 1, x, 1 )
         write(6,*) ' Positive dir derivative in projection '
         write(6,*) ' Using the backtracking step '
      else
         go to 911
      endif
c
c-----------------------------------------------------------------
c
      alpha = one
      temp1 = alpha
      ibd   = 0 
      do 60 i = 1, nsub
         k = ind(i)
         dk = d(i)
         if (nbd(k) .ne. 0) then
            if (dk .lt. zero .and. nbd(k) .le. 2) then
               temp2 = l(k) - x(k)
               if (temp2 .ge. zero) then
                  temp1 = zero
               else if (dk*alpha .lt. temp2) then
                  temp1 = temp2/dk
               endif
            else if (dk .gt. zero .and. nbd(k) .ge. 2) then
               temp2 = u(k) - x(k)
               if (temp2 .le. zero) then
                  temp1 = zero
               else if (dk*alpha .gt. temp2) then
                  temp1 = temp2/dk
               endif
            endif
            if (temp1 .lt. alpha) then
               alpha = temp1
               ibd = i
            endif
         endif
 60   continue
      
      if (alpha .lt. one) then
         dk = d(ibd)
         k = ind(ibd)
         if (dk .gt. zero) then
            x(k) = u(k)
            d(ibd) = zero
         else if (dk .lt. zero) then
            x(k) = l(k)
            d(ibd) = zero
         endif
      endif
      do 70 i = 1, nsub
         k    = ind(i)
         x(k) = x(k) + alpha*d(i)
 70   continue
cccccc
 911  continue

      if (iprint .ge. 99) write (6,1004)

 1001 format (/,'----------------SUBSM entered-----------------',/)
 1004 format (/,'----------------exit SUBSM --------------------',/)

      return

      end
