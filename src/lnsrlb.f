c> \file lnsrlb.f

c> \brief This subroutine calls subroutine dcsrch from the Minpack2 library
c>        to perform the line search.  Subroutine dscrch is safeguarded so
c>        that all trial points lie within the feasible region.
c> 
c> @param n number of parameters
c> @param l lower bounds of parameters
c> @param u upper bounds of parameters
c> @param nbd On entry nbd represents the type of bounds imposed on the
c>               variables, and must be specified as follows:
c>               nbd(i)=<ul><li>0 if x(i) is unbounded,</li>
c>                          <li>1 if x(i) has only a lower bound,</li>
c>                          <li>2 if x(i) has both lower and upper bounds, and</li>
c>                          <li>3 if x(i) has only an upper bound.</li></ul>
c>            On exit nbd is unchanged.
c> @param x position
c> @param f function value at x
c> @param fold Function value at the start of this line search (i.e. the
c>              accepted value from the previous iteration). Saved at entry
c>              so the caller can restore x, g, f if the line search fails.
c> @param gd Directional derivative g'd at the current trial step. Computed
c>           on every entry and passed to dcsrch as its `g` argument.
c> @param gdold Directional derivative at stp=0 (i.e. the initial g'd before
c>              any line-search progress this iteration). Saved on the first
c>              call and used by mainlb to test the curvature condition
c>              after the line search returns.
c> @param g Gradient of f at x.
c> @param d Search direction (z - x_current). Length-n vector; the candidate
c>          step is t + stp*d.
c> @param r Workspace: copy of g at the start of this line search. Used
c>          alongside fold to restore the previous iterate on abnormal
c>          line-search termination (mainlb does the restore from r and t).
c> @param t Workspace: copy of x at the start of this line search.
c> @param z Pre-projected candidate (cauchy/subsm output). When stp=1
c>          exactly, lnsrlb sets x := z directly; otherwise it computes
c>          x := t + stp*d.
c> @param stp Current trial step length. On the first entry of a line
c>            search lnsrlb initialises stp; subsequent dcsrch calls
c>            update it.
c> @param dnorm 2-norm of d (||d||).
c> @param dtd Squared 2-norm of d (d'd).
c> @param xstep On exit stp * ||d||, the actual length of the step in x-space.
c> @param stpmx Maximum allowed step. For unconstrained problems set to a
c>              large constant (1e10); for bounded problems lnsrlb scans
c>              the active bounds and tightens stpmx so x + stpmx*d stays
c>              feasible.
c> @param iter Outer iteration number from mainlb.
c> @param ifun On exit number of f/g evaluations performed in this line
c>             search; reset to 0 on each new line search.
c> @param iback On exit number of "backtracks" (ifun - 1). mainlb aborts
c>              the line search if iback >= 20.
c> @param nfgv Cumulative count of f/g evaluations across all iterations;
c>             incremented by 1 per evaluation requested.
c> @param info On exit 0 on success; -4 if the projected directional
c>             derivative gd is non-negative on the first call (no descent
c>             possible).
c> @param task Reverse-comm task. Initial entry: 'START'. While the line
c>             search is running, lnsrlb returns 'FG_LNSRCH' (the user
c>             evaluates f, g at the new x and re-enters with task starting
c>             with 'FG_LN'). On line-search success lnsrlb returns 'NEW_X'.
c> @param boxed .true. if every variable has both lower and upper bounds.
c>              When true, the initial trial step is unit (stp=1) regardless
c>              of d's magnitude.
c> @param cnstnd .true. if the problem has at least one bound. Controls the
c>               stpmx-from-bounds scan.
c> @param csave working array
c> @param isave working array
c> @param dsave working array
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
c     Wolfe-condition tolerances passed to dcsrch:
c       ftol -- sufficient-decrease constant (alpha in eq (2.5) of the
c               algorithm tech report; 1.0d-3 here, vs 1.0d-4 stated in
c               that report). The looser value 1.0d-3 matches the
c               implementation that ships with Algorithm 778; neither
c               the ACM paper (docs/code.pdf) nor the 2011 remark
c               documents the change explicitly.
c       gtol -- curvature constant (beta in eq (2.6) of the tech report,
c               0.9 -- matches the paper).
c       xtol -- relative-step tolerance for dcsrch's bracketing safeguard.
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
