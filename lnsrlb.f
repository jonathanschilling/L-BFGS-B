c> \file lnsrlb.f

      subroutine lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t,
     +                  z, stp, dnorm, dtd, xstep, stpmx, iter, ifun,
     +                  iback, nfgv, info, task, boxed, cnstnd, csave,
     +                  isave, dsave)

      character*60     task, csave
      logical          boxed, cnstnd
      integer          n, iter, ifun, iback, nfgv, info,
     +                 nbd(n), isave(2)
      double precision f, fold, gd, gdold, stp, dnorm, dtd, xstep,
     +                 stpmx, x(n), l(n), u(n), g(n), d(n), r(n), t(n),
     +                 z(n), dsave(13)
c     **********
c
c     Subroutine lnsrlb
c
c     This subroutine calls subroutine dcsrch from the Minpack2 library
c       to perform the line search.  Subroutine dscrch is safeguarded so
c       that all trial points lie within the feasible region.
c
c     Subprograms called:
c
c       Minpack2 Library ... dcsrch.
c
c       Linpack ... dtrsl, ddot.
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
c     **********

      integer          i
      double           precision ddot,a1,a2
      double precision one,zero,big
      parameter        (one=1.0d0,zero=0.0d0,big=1.0d+10)
      double precision ftol,gtol,xtol
      parameter        (ftol=1.0d-3,gtol=0.9d0,xtol=0.1d0)

      if (task(1:5) .eq. 'FG_LN') goto 556

      dtd = ddot(n,d,1,d,1)
      dnorm = sqrt(dtd)

c     Determine the maximum step length.

      stpmx = big
      if (cnstnd) then
         if (iter .eq. 0) then
            stpmx = one
         else
            do 43 i = 1, n
               a1 = d(i)
               if (nbd(i) .ne. 0) then
                  if (a1 .lt. zero .and. nbd(i) .le. 2) then
                     a2 = l(i) - x(i)
                     if (a2 .ge. zero) then
                        stpmx = zero
                     else if (a1*stpmx .lt. a2) then
                        stpmx = a2/a1
                     endif
                  else if (a1 .gt. zero .and. nbd(i) .ge. 2) then
                     a2 = u(i) - x(i)
                     if (a2 .le. zero) then
                        stpmx = zero
                     else if (a1*stpmx .gt. a2) then
                        stpmx = a2/a1
                     endif
                  endif
               endif
  43        continue
         endif
      endif
 
      if (iter .eq. 0 .and. .not. boxed) then
         stp = min(one/dnorm, stpmx)
      else
         stp = one
      endif 

      call dcopy(n,x,1,t,1)
      call dcopy(n,g,1,r,1)
      fold = f
      ifun = 0
      iback = 0
      csave = 'START'
 556  continue
      gd = ddot(n,g,1,d,1)
      if (ifun .eq. 0) then
         gdold=gd
         if (gd .ge. zero) then
c                               the directional derivative >=0.
c                               Line search is impossible.
            write(6,*)' ascent direction in projection gd = ', gd
            info = -4
            return
         endif
      endif

      call dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave)

      xstep = stp*dnorm
      if (csave(1:4) .ne. 'CONV' .and. csave(1:4) .ne. 'WARN') then
         task = 'FG_LNSRCH'
         ifun = ifun + 1
         nfgv = nfgv + 1
         iback = ifun - 1 
         if (stp .eq. one) then
            call dcopy(n,z,1,x,1)
         else
            do 41 i = 1, n
               x(i) = stp*d(i) + t(i)
  41        continue
         endif
      else
         task = 'NEW_X'
      endif

      return

      end
