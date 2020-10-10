c> \file cauchy.f

c> \brief For given x, l, u, g (with sbgnrm > 0), and a limited memory
c>        BFGS matrix B defined in terms of matrices WY, WS, WT, and
c>        scalars head, col, and theta, this subroutine computes the
c>        generalized Cauchy point (GCP), defined as the first local
c>        minimizer of the quadratic
c>
c>                   Q(x + s) = g's + 1/2 s'Bs
c>
c>        along the projected gradient direction P(x-tg,l,u).
c>        The routine returns the GCP in xcp. 
c>
c> @param n On entry n is the dimension of the problem.<br/>
c>          On exit n is unchanged.
c>
c> @param x On entry x is the starting point for the GCP computation.<br/>
c>          On exit x is unchanged.
c>
c> @param l On entry l is the lower bound of x.<br/>
c>          On exit l is unchanged.
c>
c> @param u On entry u is the upper bound of x.<br/>
c>          On exit u is unchanged.
c>
c> @param nbd On entry nbd represents the type of bounds imposed on the
c>            variables, and must be specified as follows:
c>            nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
c>                       <li>1 if x(i) has only a lower bound,</li>
c>                       <li>2 if x(i) has both lower and upper bounds, and</li>
c>                       <li>3 if x(i) has only an upper bound.</li></ul>
c>            On exit nbd is unchanged.
c>
c> @param g On entry g is the gradient of f(x). g must be a nonzero vector.<br/>
c>          On exit g is unchanged.
c>
c> @param iorder iorder will be used to store the breakpoints in the piecewise
c>               linear path and free variables encountered.<br/>
c>               On exit,<ul><li>iorder(1),...,iorder(nleft)
c>                               are indices of breakpoints
c>                               which have not been encountered;</li>
c>                           <li>iorder(nleft+1),...,iorder(nbreak)
c>                               are indices of
c>                               encountered breakpoints; and</li>
c>                           <li>iorder(nfree),...,iorder(n)
c>                               are indices of variables which
c>                               have no bound constraits along the search direction.</li></ul>
c>
c> @param iwhere On entry iwhere indicates only the permanently fixed (iwhere=3)
c>                 or free (iwhere= -1) components of x.<br/>
c>               On exit iwhere records the status of the current x variables.
c>                 iwhere(i)=<ul><li>-3  if x(i) is free and has bounds, but is not moved</li>
c>                               <li> 0   if x(i) is free and has bounds, and is moved</li>
c>                               <li> 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)</li>
c>                               <li> 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)</li>
c>                               <li> 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)</li>
c>                               <li>-1  if x(i) is always free, i.e., it has no bounds.</li></ul>
c>
c> @param t working array; will be used to store the break points.
c>
c> @param d the Cauchy direction P(x-tg)-x
c>
c> @param xcp returns the GCP on exit
c>
c> @param m On entry m is the maximum number of variable metric corrections 
c>            used to define the limited memory matrix.<br/>
c>          On exit m is unchanged.
c>
c> @param ws On entry this stores S, a set of s-vectors, that defines the
c>           limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param wy On entry this stores Y, a set of y-vectors, that defines the
c>           limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param sy On entry this stores S'Y, that defines the
c>           limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param wt On entry this stores the
c>           Cholesky factorization of (theta*S'S+LD^(-1)L'), that defines the
c>           limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
c>              On exit theta is unchanged.
c>
c> @param col On entry col is the actual number of variable metric
c>              corrections stored so far.<br/>
c>            On exit col is unchanged.
c>
c> @param head On entry head is the location of the first s-vector (or y-vector)
c>               in S (or Y).<br/>
c>             On exit col is unchanged.
c>
c> @param p will be used to store the vector p = W^(T)d.
c>
c> @param c will be used to store the vector c = W^(T)(xcp-x).
c>
c> @param wbp will be used to store the row of W corresponding
c>              to a breakpoint.
c>
c> @param v working array
c>
c> @param nseg On exit nseg records the number of quadratic segments explored
c>               in searching for the GCP.
c> 
c> @param iprint variable that must be set by the user.<br/>
c>       It controls the frequency and type of output generated:
c>       <ul><li>iprint<0    no output is generated;</li>
c>           <li>iprint=0    print only one line at the last iteration;</li>
c>           <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
c>           <li>iprint=99   print details of every iteration except n-vectors;</li>
c>           <li>iprint=100  print also the changes of active set and final x;</li>
c>           <li>iprint>100  print details of every iteration including x and g;</li></ul>
c>       When iprint > 0, the file iterate.dat will be created to
c>                        summarize the iteration.
c>
c> @param sbgnrm On entry sbgnrm is the norm of the projected gradient at x.<br/>
c>               On exit sbgnrm is unchanged.
c>
c> @param info On entry info is 0.
c>             On exit info =<ul><li>0       for normal return,</li>
c>                               <li>= nonzero for abnormal return when the the system
c>                                     used in routine bmv is singular.</li></ul>
c> @param epsmch machine precision epsilon
      subroutine cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, 
     +                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, 
     +                  v, nseg, iprint, sbgnrm, info, epsmch)
      implicit none
      integer          n, m, head, col, nseg, iprint, info, 
     +                 nbd(n), iorder(n), iwhere(n)
      double precision theta, epsmch,
     +                 x(n), l(n), u(n), g(n), t(n), d(n), xcp(n),
     +                 wy(n, col), ws(n, col), sy(m, m),
     +                 wt(m, m), p(2*m), c(2*m), wbp(2*m), v(2*m)
c
c     References:
c
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c       Subroutines for Large Scale Bound Constrained Optimization''
c       Tech. Report, NAM-11, EECS Department, Northwestern University,
c       1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
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

      logical          xlower,xupper,bnded
      integer          i,j,col2,nfree,nbreak,pointr,
     +                 ibp,nleft,ibkmin,iter
      double precision f1,f2,dt,dtm,tsum,dibp,zibp,dibp2,bkmin,
     +                 tu,tl,wmc,wmp,wmw,ddot,tj,tj0,neggi,sbgnrm,
     +                 f2_org
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
 
c     Check the status of the variables, reset iwhere(i) if necessary;
c       compute the Cauchy direction d and the breakpoints t; initialize
c       the derivative f1 and the vector p = W'd (for theta = 1).
 
      if (sbgnrm .le. zero) then
         if (iprint .ge. 0) write (6,*) 'Subgnorm = 0.  GCP = X.'
         call dcopy(n,x,1,xcp,1)
         return
      endif 
      bnded = .true.
      nfree = n + 1
      nbreak = 0
      ibkmin = 0
      bkmin = zero
      col2 = 2*col
      f1 = zero
      if (iprint .ge. 99) write (6,3010)

c     We set p to zero and build it up as we determine d.

      do 20 i = 1, col2
         p(i) = zero
  20  continue 

c     In the following loop we determine for each variable its bound
c        status and its breakpoint, and update p accordingly.
c        Smallest breakpoint is identified.

      do 50 i = 1, n 
         neggi = -g(i)      
         if (iwhere(i) .ne. 3 .and. iwhere(i) .ne. -1) then
c             if x(i) is not a constant and has bounds,
c             compute the difference between x(i) and its bounds.
            if (nbd(i) .le. 2) tl = x(i) - l(i)
            if (nbd(i) .ge. 2) tu = u(i) - x(i)

c           If a variable is close enough to a bound
c             we treat it as at bound.
            xlower = nbd(i) .le. 2 .and. tl .le. zero
            xupper = nbd(i) .ge. 2 .and. tu .le. zero

c              reset iwhere(i).
            iwhere(i) = 0
            if (xlower) then
               if (neggi .le. zero) iwhere(i) = 1
            else if (xupper) then
               if (neggi .ge. zero) iwhere(i) = 2
            else
               if (abs(neggi) .le. zero) iwhere(i) = -3
            endif
         endif 
         pointr = head
         if (iwhere(i) .ne. 0 .and. iwhere(i) .ne. -1) then
            d(i) = zero
         else
            d(i) = neggi
            f1 = f1 - neggi*neggi
c             calculate p := p - W'e_i* (g_i).
            do 40 j = 1, col
               p(j) = p(j) +  wy(i,pointr)* neggi
               p(col + j) = p(col + j) + ws(i,pointr)*neggi
               pointr = mod(pointr,m) + 1
  40        continue 
            if (nbd(i) .le. 2 .and. nbd(i) .ne. 0
     +                        .and. neggi .lt. zero) then
c                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tl/(-neggi)
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else if (nbd(i) .ge. 2 .and. neggi .gt. zero) then
c                                 x(i) + d(i) is bounded; compute t(i).
               nbreak = nbreak + 1
               iorder(nbreak) = i
               t(nbreak) = tu/neggi
               if (nbreak .eq. 1 .or. t(nbreak) .lt. bkmin) then
                  bkmin = t(nbreak)
                  ibkmin = nbreak
               endif
            else
c                x(i) + d(i) is not bounded.
               nfree = nfree - 1
               iorder(nfree) = i
               if (abs(neggi) .gt. zero) bnded = .false.
            endif
         endif
  50  continue 
 
c     The indices of the nonzero components of d are now stored
c       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
c       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
 
      if (theta .ne. one) then
c                   complete the initialization of p for theta not= one.
         call dscal(col,theta,p(col+1),1)
      endif
 
c     Initialize GCP xcp = x.

      call dcopy(n,x,1,xcp,1)

      if (nbreak .eq. 0 .and. nfree .eq. n + 1) then
c                  is a zero vector, return with the initial xcp as GCP.
         if (iprint .gt. 100) write (6,1010) (xcp(i), i = 1, n)
         return
      endif    
 
c     Initialize c = W'(xcp - x) = 0.
  
      do 60 j = 1, col2
         c(j) = zero
  60  continue 
 
c     Initialize derivative f2.
 
      f2 =  -theta*f1 
      f2_org  =  f2
      if (col .gt. 0) then
         call bmv(m,sy,wt,col,p,v,info)
         if (info .ne. 0) return
         f2 = f2 - ddot(col2,v,1,p,1)
      endif
      dtm = -f1/f2
      tsum = zero
      nseg = 1
      if (iprint .ge. 99) 
     +   write (6,*) 'There are ',nbreak,'  breakpoints '
 
c     If there are no breakpoints, locate the GCP and return. 
 
      if (nbreak .eq. 0) goto 888
             
      nleft = nbreak
      iter = 1
 
 
      tj = zero
 
c------------------- the beginning of the loop -------------------------
 
 777  continue
 
c     Find the next smallest breakpoint;
c       compute dt = t(nleft) - t(nleft + 1).
 
      tj0 = tj
      if (iter .eq. 1) then
c         Since we already have the smallest breakpoint we need not do
c         heapsort yet. Often only one breakpoint is used and the
c         cost of heapsort is avoided.
         tj = bkmin
         ibp = iorder(ibkmin)
      else
         if (iter .eq. 2) then
c             Replace the already used smallest breakpoint with the
c             breakpoint numbered nbreak > nlast, before heapsort call.
            if (ibkmin .ne. nbreak) then
               t(ibkmin) = t(nbreak)
               iorder(ibkmin) = iorder(nbreak)
            endif 
c        Update heap structure of breakpoints
c           (if iter=2, initialize heap).
         endif
         call hpsolb(nleft,t,iorder,iter-2)
         tj = t(nleft)
         ibp = iorder(nleft)  
      endif 
         
      dt = tj - tj0
 
      if (dt .ne. zero .and. iprint .ge. 100) then
         write (6,4011) nseg,f1,f2
         write (6,5010) dt
         write (6,6010) dtm
      endif          
 
c     If a minimizer is within this interval, locate the GCP and return. 
 
      if (dtm .lt. dt) goto 888
 
c     Otherwise fix one variable and
c       reset the corresponding component of d to zero.
    
      tsum = tsum + dt
      nleft = nleft - 1
      iter = iter + 1
      dibp = d(ibp)
      d(ibp) = zero
      if (dibp .gt. zero) then
         zibp = u(ibp) - x(ibp)
         xcp(ibp) = u(ibp)
         iwhere(ibp) = 2
      else
         zibp = l(ibp) - x(ibp)
         xcp(ibp) = l(ibp)
         iwhere(ibp) = 1
      endif
      if (iprint .ge. 100) write (6,*) 'Variable  ',ibp,'  is fixed.'
      if (nleft .eq. 0 .and. nbreak .eq. n) then
c                                             all n variables are fixed,
c                                                return with xcp as GCP.
         dtm = dt
         goto 999
      endif
 
c     Update the derivative information.
 
      nseg = nseg + 1
      dibp2 = dibp**2
 
c     Update f1 and f2.
 
c        temporarily set f1 and f2 for col=0.
      f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp
      f2 = f2 - theta*dibp2

      if (col .gt. 0) then
c                          update c = c + dt*p.
         call daxpy(col2,dt,p,1,c,1)
 
c           choose wbp,
c           the row of W corresponding to the breakpoint encountered.
         pointr = head
         do 70 j = 1,col
            wbp(j) = wy(ibp,pointr)
            wbp(col + j) = theta*ws(ibp,pointr)
            pointr = mod(pointr,m) + 1
  70     continue 
 
c           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
         call bmv(m,sy,wt,col,wbp,v,info)
         if (info .ne. 0) return
         wmc = ddot(col2,c,1,v,1)
         wmp = ddot(col2,p,1,v,1) 
         wmw = ddot(col2,wbp,1,v,1)
 
c           update p = p - dibp*wbp. 
         call daxpy(col2,-dibp,wbp,1,p,1)
 
c           complete updating f1 and f2 while col > 0.
         f1 = f1 + dibp*wmc
         f2 = f2 + 2.0d0*dibp*wmp - dibp2*wmw
      endif

      f2 = max(epsmch*f2_org,f2)
      if (nleft .gt. 0) then
         dtm = -f1/f2
         goto 777
c                 to repeat the loop for unsearched intervals. 
      else if(bnded) then
         f1 = zero
         f2 = zero
         dtm = zero
      else
         dtm = -f1/f2
      endif 

c------------------- the end of the loop -------------------------------
 
 888  continue
      if (iprint .ge. 99) then
         write (6,*)
         write (6,*) 'GCP found in this segment'
         write (6,4010) nseg,f1,f2
         write (6,6010) dtm
      endif 
      if (dtm .le. zero) dtm = zero
      tsum = tsum + dtm
 
c     Move free variables (i.e., the ones w/o breakpoints) and 
c       the variables whose breakpoints haven't been reached.
 
      call daxpy(n,tsum,d,1,xcp,1)
 
 999  continue
 
c     Update c = c + dtm*p = W'(x^c - x) 
c       which will be used in computing r = Z'(B(x^c - x) + g).
 
      if (col .gt. 0) call daxpy(col2,dtm,p,1,c,1)
      if (iprint .gt. 100) write (6,1010) (xcp(i),i = 1,n)
      if (iprint .ge. 99) write (6,2010)

 1010 format ('Cauchy X =  ',/,(4x,1p,6(1x,d11.4)))
 2010 format (/,'---------------- exit CAUCHY----------------------',/)
 3010 format (/,'---------------- CAUCHY entered-------------------')
 4010 format ('Piece    ',i3,' --f1, f2 at start point ',1p,2(1x,d11.4))
 4011 format (/,'Piece    ',i3,' --f1, f2 at start point ',
     +        1p,2(1x,d11.4))
 5010 format ('Distance to the next break point =  ',1p,d11.4)
 6010 format ('Distance to the stationary point =  ',1p,d11.4) 
 
      return
 
      end
