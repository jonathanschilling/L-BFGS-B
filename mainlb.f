c> \file mainlb.f

c> \brief This subroutine solves bound constrained optimization problems by
c>        using the compact formula of the limited memory BFGS updates.
c>
c> This subroutine solves bound constrained optimization problems by
c> using the compact formula of the limited memory BFGS updates.
c>   
c> @param n On entry n is the number of variables.<br/>
c>          On exit n is unchanged.
c>
c> @param m On entry m is the maximum number of variable metric
c>             corrections allowed in the limited memory matrix.<br/>
c>          On exit m is unchanged.
c>
c> @param x On entry x is an approximation to the solution.<br/>
c>          On exit x is the current approximation.
c>
c> @param l On entry l is the lower bound of x.<br/>
c>          On exit l is unchanged.
c>
c> @param u On entry u is the upper bound of x.<br/>
c>          On exit u is unchanged.
c>
c> @param nbd On entry nbd represents the type of bounds imposed on the
c>               variables, and must be specified as follows:
c>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
c>                          <li>1 if x(i) has only a lower bound,</li>
c>                          <li>2 if x(i) has both lower and upper bounds,</li>
c>                          <li>3 if x(i) has only an upper bound.</li></ul>
c>            On exit nbd is unchanged.
c>
c> @param f On first entry f is unspecified.<br/>
c>          On final exit f is the value of the function at x.
c>
c> @param g On first entry g is unspecified.<br/>
c>          On final exit g is the value of the gradient at x.
c>
c> @param factr On entry factr >= 0 is specified by the user.  The iteration
c>                 will stop when<br/>
c>                        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch<br/>
c>                 where epsmch is the machine precision, which is automatically
c>                 generated by the code.<br/>
c>              On exit factr is unchanged.
c>
c> @param pgtol On entry pgtol >= 0 is specified by the user.  The iteration
c>                 will stop when<br/>
c>                        max{|proj g_i | i = 1, ..., n} <= pgtol<br/>
c>                 where pg_i is the ith component of the projected gradient.<br/>
c>              On exit pgtol is unchanged.
c> 
c> @param ws On entry this stores S, a set of s-vectors, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param wy On entry this stores Y, a set of y-vectors, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param sy On entry this stores S'Y, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param ss On entry this stores S'S, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param wt On entry this stores the
c>              Cholesky factorization of (theta*S'S+LD^(-1)L'), that defines the
c>              limited memory BFGS matrix. See eq. (2.26) in [3].<br/>
c>           On exit this array is unchanged.
c>
c> @param wn working array used to store the LEL^T factorization of the indefinite matrix
c>                     K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c>                         [L_a -R_z           theta*S'AA'S ]<br/>
c>           where     E = [-I  0]
c>                         [ 0  I]
c>
c> @param snd working array used to store the lower triangular part of
c>                      N = [Y' ZZ'Y   L_a'+R_z']
c>                          [L_a +R_z  S'AA'S   ]
c>                 
c> @param z working array used at different times to store the Cauchy point and
c>          the Newton point.
c> @param r working array
c> @param d working array
c> @param t working array
c> @param xp working array used to safeguard the projected Newton direction
c> @param wa working array
c>
c> @param index In subroutine freev, index is used to store the free and fixed
c>              variables at the Generalized Cauchy Point (GCP).
c>
c> @param iwhere working array used to record
c>               the status of the vector x for GCP computation.<br/>
c>               iwhere(i)=<ul><li>0 or -3 if x(i) is free and has bounds,</li>
c>                             <li> 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)</li>
c>                             <li> 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)</li>
c>                             <li> 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)</li>
c>                             <li>-1       if x(i) is always free, i.e., no bounds on it.</li></ul>
c>
c> @param indx2 working array<br/>
c>              Within subroutine cauchy, indx2 corresponds to the array iorder.<br/>
c>              In subroutine freev, a list of variables entering and leaving
c>              the free set is stored in indx2, and it is passed on to
c>              subroutine formk with this information.
c>
c> @param task working string indicating
c>             the current job when entering and leaving this subroutine.
c>
c> @param iprint It controls the frequency and type of output generated:<ul>
c>               <li>iprint<0    no output is generated;</li>
c>               <li>iprint=0    print only one line at the last iteration;</li>
c>               <li>0<iprint<99 print also f and |proj g| every iprint iterations;</li>
c>               <li>iprint=99   print details of every iteration except n-vectors;</li>
c>               <li>iprint=100  print also the changes of active set and final x;</li>
c>               <li>iprint>100  print details of every iteration including x and g;</li></ul>
c>               When iprint > 0, the file iterate.dat will be created to
c>                                summarize the iteration.
c>
c> @param csave working string
c>
c> @param lsave working array
c>
c> @param isave working array
c>
c> @param dsave working array
      subroutine mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy,
     +                  sy, ss, wt, wn, snd, z, r, d, t, xp, wa, 
     +                  index, iwhere, indx2, task,
     +                  iprint, csave, lsave, isave, dsave)
      implicit none
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, nbd(n), index(n),
     +                 iwhere(n), indx2(n), isave(23)
      double precision f, factr, pgtol,
     +                 x(n), l(n), u(n), g(n), z(n), r(n), d(n), t(n), 
c-jlm-jn
     +                 xp(n), 
     +                 wa(8*m), 
     +                 ws(n, m), wy(n, m), sy(m, m), ss(m, m), 
     +                 wt(m, m), wn(2*m, 2*m), snd(2*m, 2*m), dsave(29)
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
c       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
c       Quasi-Newton Matrices and their use in Limited Memory Methods'',
c       Mathematical Programming 63 (1994), no. 4, pp. 129-156.
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
 
      logical          prjctd,cnstnd,boxed,updatd,wrk
      character*3      word
      integer          i,k,nintol,itfile,iback,nskip,
     +                 head,col,iter,itail,iupdat,
     +                 nseg,nfgv,info,ifun,
     +                 iword,nfree,nact,ileave,nenter
      double precision theta,fold,ddot,dr,rr,tol,
     +                 xstep,sbgnrm,ddum,dnorm,dtd,epsmch,
     +                 cpu1,cpu2,cachyt,sbtime,lnscht,time1,time2,
     +                 gd,gdold,stp,stpmx,time
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)
      
      if (task .eq. 'START') then

         epsmch = epsilon(one)

         call timer(time1)

c        Initialize counters and scalars when task='START'.

c           for the limited memory BFGS matrices:
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         iback  = 0
         itail  = 0
         iword  = 0
         nact   = 0
         ileave = 0
         nenter = 0
         fold   = zero
         dnorm  = zero
         cpu1   = zero
         gd     = zero
         stpmx  = zero
         sbgnrm = zero
         stp    = zero
         gdold  = zero
         dtd    = zero

c           for operation counts:
         iter   = 0
         nfgv   = 0
         nseg   = 0
         nintol = 0
         nskip  = 0
         nfree  = n
         ifun   = 0
c           for stopping tolerance:
         tol = factr*epsmch

c           for measuring running time:
         cachyt = 0
         sbtime = 0
         lnscht = 0
 
c           'word' records the status of subspace solutions.
         word = '---'

c           'info' records the termination information.
         info = 0

         itfile = 8
         if (iprint .ge. 1) then
c                                open a summary file 'iterate.dat'
            open (itfile, file = 'iterate.dat', status = 'unknown')
         endif            

c        Check the input arguments for errors.

         call errclb(n,m,factr,l,u,nbd,task,info,k)
         if (task(1:5) .eq. 'ERROR') then
            call prn3lb(n,x,f,task,iprint,info,itfile,
     +                  iter,nfgv,nintol,nskip,nact,sbgnrm,
     +                  zero,nseg,word,iback,stp,xstep,k,
     +                  cachyt,sbtime,lnscht)
            return
         endif

         call prn1lb(n,m,l,u,x,iprint,itfile,epsmch)
 
c        Initialize iwhere & project x onto the feasible set.
 
         call active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed) 

c        The end of the initialization.

      else
c          restore local variables.

         prjctd = lsave(1)
         cnstnd = lsave(2)
         boxed  = lsave(3)
         updatd = lsave(4)

         nintol = isave(1)
         itfile = isave(3)
         iback  = isave(4)
         nskip  = isave(5)
         head   = isave(6)
         col    = isave(7)
         itail  = isave(8)
         iter   = isave(9)
         iupdat = isave(10)
         nseg   = isave(12)
         nfgv   = isave(13)
         info   = isave(14)
         ifun   = isave(15)
         iword  = isave(16)
         nfree  = isave(17)
         nact   = isave(18)
         ileave = isave(19)
         nenter = isave(20)

         theta  = dsave(1)
         fold   = dsave(2)
         tol    = dsave(3)
         dnorm  = dsave(4)
         epsmch = dsave(5)
         cpu1   = dsave(6)
         cachyt = dsave(7)
         sbtime = dsave(8)
         lnscht = dsave(9)
         time1  = dsave(10)
         gd     = dsave(11)
         stpmx  = dsave(12)
         sbgnrm = dsave(13)
         stp    = dsave(14)
         gdold  = dsave(15)
         dtd    = dsave(16)
   
c        After returning from the driver go to the point where execution
c        is to resume.

         if (task(1:5) .eq. 'FG_LN') goto 666
         if (task(1:5) .eq. 'NEW_X') goto 777
         if (task(1:5) .eq. 'FG_ST') goto 111
         if (task(1:4) .eq. 'STOP') then
            if (task(7:9) .eq. 'CPU') then
c                                          restore the previous iterate.
               call dcopy(n,t,1,x,1)
               call dcopy(n,r,1,g,1)
               f = fold
            endif
            goto 999
         endif
      endif 

c     Compute f0 and g0.

      task = 'FG_START' 
c          return to the driver to calculate f and g; reenter at 111.
      goto 1000
 111  continue
      nfgv = 1
 
c     Compute the infinity norm of the (-) projected gradient.
 
      call projgr(n,l,u,nbd,x,g,sbgnrm)
  
      if (iprint .ge. 1) then
         write (6,1002) iter,f,sbgnrm
         write (itfile,1003) iter,nfgv,sbgnrm,f
      endif
      if (sbgnrm .le. pgtol) then
c                                terminate the algorithm.
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 999
      endif 
 
c ----------------- the beginning of the loop --------------------------
 
 222  continue
      if (iprint .ge. 99) write (6,1001) iter + 1
      iword = -1
c
      if (.not. cnstnd .and. col .gt. 0) then 
c                                            skip the search for GCP.
         call dcopy(n,x,1,z,1)
         wrk = updatd
         nseg = 0
         goto 333
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the Generalized Cauchy Point (GCP).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call timer(cpu1) 
      call cauchy(n,x,l,u,nbd,g,indx2,iwhere,t,d,z,
     +            m,wy,ws,sy,wt,theta,col,head,
     +            wa(1),wa(2*m+1),wa(4*m+1),wa(6*m+1),nseg,
     +            iprint, sbgnrm, info, epsmch)
      if (info .ne. 0) then 
c         singular triangular system detected; refresh the lbfgs memory.
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         cachyt = cachyt + cpu2 - cpu1
         goto 222
      endif
      call timer(cpu2) 
      cachyt = cachyt + cpu2 - cpu1
      nintol = nintol + nseg

c     Count the entering and leaving variables for iter > 0; 
c     find the index set of free and active variables at the GCP.

      call freev(n,nfree,index,nenter,ileave,indx2,
     +           iwhere,wrk,updatd,cnstnd,iprint,iter)
      nact = n - nfree

 333  continue
 
c     If there are no free variables or B=theta*I, then
c                                        skip the subspace minimization.
 
      if (nfree .eq. 0 .or. col .eq. 0) goto 555
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subspace minimization.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call timer(cpu1) 

c     Form  the LEL^T factorization of the indefinite
c       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c                     [L_a -R_z           theta*S'AA'S ]
c       where     E = [-I  0]
c                     [ 0  I]

      if (wrk) call formk(n,nfree,index,nenter,ileave,indx2,iupdat,
     +                 updatd,wn,snd,m,ws,wy,sy,theta,col,head,info)
      if (info .ne. 0) then
c          nonpositive definiteness in Cholesky factorization;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1006)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         sbtime = sbtime + cpu2 - cpu1 
         goto 222
      endif 

c        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
c                                                   from 'cauchy').
      call cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa,index,
     +           theta,col,head,nfree,cnstnd,info)
      if (info .ne. 0) goto 444

c-jlm-jn   call the direct method. 

      call subsm( n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy,
     +           theta, x, g, col, head, iword, wa, wn, iprint, info)
 444  continue
      if (info .ne. 0) then 
c          singular triangular system detected;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1005)
         info   = 0
         col    = 0
         head   = 1
         theta  = one
         iupdat = 0
         updatd = .false.
         call timer(cpu2) 
         sbtime = sbtime + cpu2 - cpu1 
         goto 222
      endif
 
      call timer(cpu2) 
      sbtime = sbtime + cpu2 - cpu1 
 555  continue
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Line search and optimality tests.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
c     Generate the search direction d:=z-x.

      do 40 i = 1, n
         d(i) = z(i) - x(i)
  40  continue
      call timer(cpu1) 
 666  continue
      call lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm,
     +            dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,
     +            boxed,cnstnd,csave,isave(22),dsave(17))
      if (info .ne. 0 .or. iback .ge. 20) then
c          restore the previous iterate.
         call dcopy(n,t,1,x,1)
         call dcopy(n,r,1,g,1)
         f = fold
         if (col .eq. 0) then
c             abnormal termination.
            if (info .eq. 0) then
               info = -9
c                restore the actual number of f and g evaluations etc.
               nfgv = nfgv - 1
               ifun = ifun - 1
               iback = iback - 1
            endif
            task = 'ABNORMAL_TERMINATION_IN_LNSRCH'
            iter = iter + 1
            goto 999
         else
c             refresh the lbfgs memory and restart the iteration.
            if(iprint .ge. 1) write (6, 1008)
            if (info .eq. 0) nfgv = nfgv - 1
            info   = 0
            col    = 0
            head   = 1
            theta  = one
            iupdat = 0
            updatd = .false.
            task   = 'RESTART_FROM_LNSRCH'
            call timer(cpu2)
            lnscht = lnscht + cpu2 - cpu1
            goto 222
         endif
      else if (task(1:5) .eq. 'FG_LN') then
c          return to the driver for calculating f and g; reenter at 666.
         goto 1000
      else 
c          calculate and print out the quantities related to the new X.
         call timer(cpu2) 
         lnscht = lnscht + cpu2 - cpu1
         iter = iter + 1
 
c        Compute the infinity norm of the projected (-)gradient.
 
         call projgr(n,l,u,nbd,x,g,sbgnrm)
 
c        Print iteration information.

         call prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,
     +               sbgnrm,nseg,word,iword,iback,stp,xstep)
         goto 1000
      endif
 777  continue

c     Test for termination.

      if (sbgnrm .le. pgtol) then
c                                terminate the algorithm.
         task = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
         goto 999
      endif 

      ddum = max(abs(fold), abs(f), one)
      if ((fold - f) .le. tol*ddum) then
c                                        terminate the algorithm.
         task = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
         if (iback .ge. 10) info = -5
c           i.e., to issue a warning if iback>10 in the line search.
         goto 999
      endif 

c     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
 
      do 42 i = 1, n
         r(i) = g(i) - r(i)
  42  continue
      rr = ddot(n,r,1,r,1)
      if (stp .eq. one) then  
         dr = gd - gdold
         ddum = -gdold
      else
         dr = (gd - gdold)*stp
         call dscal(n,stp,d,1)
         ddum = -gdold*stp
      endif
 
      if (dr .le. epsmch*ddum) then
c                            skip the L-BFGS update.
         nskip = nskip + 1
         updatd = .false.
         if (iprint .ge. 1) write (6,1004) dr, ddum
         goto 888
      endif 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update the L-BFGS matrix.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      updatd = .true.
      iupdat = iupdat + 1

c     Update matrices WS and WY and form the middle matrix in B.

      call matupd(n,m,ws,wy,sy,ss,d,r,itail,
     +            iupdat,col,head,theta,rr,dr,stp,dtd)

c     Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
c        Store T in the upper triangular of the array wt;
c        Cholesky factorize T to J*J' with
c           J' stored in the upper triangular of wt.

      call formt(m,wt,sy,ss,col,theta,info)
 
      if (info .ne. 0) then 
c          nonpositive definiteness in Cholesky factorization;
c          refresh the lbfgs memory and restart the iteration.
         if(iprint .ge. 1) write (6, 1007)
         info = 0
         col = 0
         head = 1
         theta = one
         iupdat = 0
         updatd = .false.
         goto 222
      endif

c     Now the inverse of the middle matrix in B is

c       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
c       [ -L*D^(-1/2)   J ] [  0        J'          ]

 888  continue
 
c -------------------- the end of the loop -----------------------------
 
      goto 222
 999  continue
      call timer(time2)
      time = time2 - time1
      call prn3lb(n,x,f,task,iprint,info,itfile,
     +            iter,nfgv,nintol,nskip,nact,sbgnrm,
     +            time,nseg,word,iback,stp,xstep,k,
     +            cachyt,sbtime,lnscht)
 1000 continue

c     Save local variables.

      lsave(1)  = prjctd
      lsave(2)  = cnstnd
      lsave(3)  = boxed
      lsave(4)  = updatd

      isave(1)  = nintol 
      isave(3)  = itfile 
      isave(4)  = iback 
      isave(5)  = nskip 
      isave(6)  = head 
      isave(7)  = col 
      isave(8)  = itail 
      isave(9)  = iter 
      isave(10) = iupdat 
      isave(12) = nseg
      isave(13) = nfgv 
      isave(14) = info 
      isave(15) = ifun 
      isave(16) = iword 
      isave(17) = nfree 
      isave(18) = nact 
      isave(19) = ileave 
      isave(20) = nenter 

      dsave(1)  = theta 
      dsave(2)  = fold 
      dsave(3)  = tol 
      dsave(4)  = dnorm 
      dsave(5)  = epsmch 
      dsave(6)  = cpu1 
      dsave(7)  = cachyt 
      dsave(8)  = sbtime 
      dsave(9)  = lnscht 
      dsave(10) = time1 
      dsave(11) = gd 
      dsave(12) = stpmx 
      dsave(13) = sbgnrm
      dsave(14) = stp
      dsave(15) = gdold
      dsave(16) = dtd  

 1001 format (//,'ITERATION ',i5)
 1002 format
     +  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 1003 format (2(1x,i4),5x,'-',5x,'-',3x,'-',5x,'-',5x,'-',8x,'-',3x,
     +        1p,2(1x,d10.3))
 1004 format ('  ys=',1p,e10.3,'  -gs=',1p,e10.3,' BFGS update SKIPPED')
 1005 format (/, 
     +' Singular triangular system detected;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1006 format (/, 
     +' Nonpositive definiteness in Cholesky factorization in formk;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1007 format (/, 
     +' Nonpositive definiteness in Cholesky factorization in formt;',/,
     +'   refresh the lbfgs memory and restart the iteration.')
 1008 format (/, 
     +' Bad direction in the line search;',/,
     +'   refresh the lbfgs memory and restart the iteration.')

      return   

      end
