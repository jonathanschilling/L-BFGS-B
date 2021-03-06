c> \file dcsrch.f

c> \brief This subroutine finds a step that satisfies a sufficient
c>        decrease condition and a curvature condition.
c>
c> This subroutine finds a step that satisfies a sufficient
c> decrease condition and a curvature condition.
c>
c> Each call of the subroutine updates an interval with 
c> endpoints stx and sty. The interval is initially chosen 
c> so that it contains a minimizer of the modified function
c>
c>       psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
c>
c> If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c> interval is chosen so that it contains a minimizer of f. 
c>
c> The algorithm is designed to find a step that satisfies 
c> the sufficient decrease condition 
c>
c>       f(stp) <= f(0) + ftol*stp*f'(0),
c>
c> and the curvature condition
c>
c>       abs(f'(stp)) <= gtol*abs(f'(0)).
c>
c> If ftol is less than gtol and if, for example, the function
c> is bounded below, then there is always a step which satisfies
c> both conditions. 
c>
c> If no step can be found that satisfies both conditions, then 
c> the algorithm stops with a warning. In this case stp only 
c> satisfies the sufficient decrease condition.
c>
c> A typical invocation of dcsrch has the following outline:
c>
c> ```Fortran
c>     task = 'START'
c>  10 continue
c>     call dcsrch( ... )
c>     if (task .eq. 'FG') then
c>       Evaluate the function and the gradient at stp 
c>     goto 10
c>     end if
c> ```
c> 
c> NOTE: The user must no alter work arrays between calls.
c>
c> @param f On initial entry f is the value of the function at 0.<br/>
c>          On subsequent entries f is the value of the 
c>             function at stp.<br/>
c>          On exit f is the value of the function at stp.
c> @param g On initial entry g is the derivative of the function at 0.<br/>
c>          On subsequent entries g is the derivative of the 
c>             function at stp.<br/>
c>          On exit g is the derivative of the function at stp.
c> @param stp On entry stp is the current estimate of a satisfactory 
c>               step. On initial entry, a positive initial estimate 
c>               must be provided.<br/>
c>            On exit stp is the current estimate of a satisfactory step
c>               if task = 'FG'. If task = 'CONV' then stp satisfies
c>               the sufficient decrease and curvature condition.
c> @param ftol On entry ftol specifies a nonnegative tolerance for the 
c>                sufficient decrease condition.<br/>
c>             On exit ftol is unchanged.
c> @param gtol On entry gtol specifies a nonnegative tolerance for the 
c>                curvature condition.<br/>
c>             On exit gtol is unchanged.
c> @param xtol On entry xtol specifies a nonnegative relative tolerance
c>                for an acceptable step. The subroutine exits with a
c>                warning if the relative difference between sty and stx
c>                is less than xtol.<br/>
c>             On exit xtol is unchanged.
c> @param stpmin On entry stpmin is a nonnegative lower bound for the step.<br/>
c>               On exit stpmin is unchanged.
c> @param stpmax On entry stpmax is a nonnegative upper bound for the step.<br/>
c>               On exit stpmax is unchanged.
c> @param task On initial entry task must be set to 'START'.<br/>
c>             On exit task indicates the required action:
c>             <ul>
c>             <li>If task(1:2) = 'FG' then evaluate the function and
c>                 derivative at stp and call dcsrch again.</li>
c>             <li>If task(1:4) = 'CONV' then the search is successful.</li>
c>             <li>If task(1:4) = 'WARN' then the subroutine is not able
c>                 to satisfy the convergence conditions. The exit value of
c>                 stp contains the best point found during the search.</li>
c>             <li>If task(1:5) = 'ERROR' then there is an error in the
c>                 input arguments.</li>
c>             </ul>
c>             On exit with convergence, a warning or an error, the
c>                variable task contains additional information.
c> @param isave work array
c> @param dsave work array
      subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
     +                  task,isave,dsave)
      character*(*) task
      integer isave(2)
      double precision f,g,stp,ftol,gtol,xtol,stpmin,stpmax
      double precision dsave(13)
c
c     MINPACK-1 Project. June 1983.
c     Argonne National Laboratory. 
c     Jorge J. More' and David J. Thuente.
c
c     MINPACK-2 Project. October 1993.
c     Argonne National Laboratory and University of Minnesota. 
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'. 
c
c     **********
      double precision zero,p5,p66
      parameter(zero=0.0d0,p5=0.5d0,p66=0.66d0)
      double precision xtrapl,xtrapu
      parameter(xtrapl=1.1d0,xtrapu=4.0d0)

      logical brackt
      integer stage
      double precision finit,ftest,fm,fx,fxm,fy,fym,ginit,gtest,
     +       gm,gx,gxm,gy,gym,stx,sty,stmin,stmax,width,width1

c     Initialization block.

      if (task(1:5) .eq. 'START') then

c        Check the input arguments for errors.

         if (stp .lt. stpmin) task = 'ERROR: STP .LT. STPMIN'
         if (stp .gt. stpmax) task = 'ERROR: STP .GT. STPMAX'
         if (g .ge. zero) task = 'ERROR: INITIAL G .GE. ZERO'
         if (ftol .lt. zero) task = 'ERROR: FTOL .LT. ZERO'
         if (gtol .lt. zero) task = 'ERROR: GTOL .LT. ZERO'
         if (xtol .lt. zero) task = 'ERROR: XTOL .LT. ZERO'
         if (stpmin .lt. zero) task = 'ERROR: STPMIN .LT. ZERO'
         if (stpmax .lt. stpmin) task = 'ERROR: STPMAX .LT. STPMIN'

c        Exit if there are errors on input.

         if (task(1:5) .eq. 'ERROR') return

c        Initialize local variables.

         brackt = .false.
         stage = 1
         finit = f
         ginit = g
         gtest = ftol*ginit
         width = stpmax - stpmin
         width1 = width/p5

c        The variables stx, fx, gx contain the values of the step, 
c        function, and derivative at the best step. 
c        The variables sty, fy, gy contain the value of the step, 
c        function, and derivative at sty.
c        The variables stp, f, g contain the values of the step, 
c        function, and derivative at stp.

         stx = zero
         fx = finit
         gx = ginit
         sty = zero
         fy = finit
         gy = ginit
         stmin = zero
         stmax = stp + xtrapu*stp
         task = 'FG'

         goto 1000

      else

c        Restore local variables.

         if (isave(1) .eq. 1) then
            brackt = .true.
         else
            brackt = .false.
         endif
         stage = isave(2) 
         ginit = dsave(1) 
         gtest = dsave(2) 
         gx = dsave(3) 
         gy = dsave(4) 
         finit = dsave(5) 
         fx = dsave(6) 
         fy = dsave(7) 
         stx = dsave(8) 
         sty = dsave(9) 
         stmin = dsave(10) 
         stmax = dsave(11) 
         width = dsave(12) 
         width1 = dsave(13) 

      endif

c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
c     algorithm enters the second stage.

      ftest = finit + stp*gtest
      if (stage .eq. 1 .and. f .le. ftest .and. g .ge. zero) 
     +   stage = 2

c     Test for warnings.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax))
     +   task = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS'
      if (brackt .and. stmax - stmin .le. xtol*stmax) 
     +   task = 'WARNING: XTOL TEST SATISFIED'
      if (stp .eq. stpmax .and. f .le. ftest .and. g .le. gtest) 
     +   task = 'WARNING: STP = STPMAX'
      if (stp .eq. stpmin .and. (f .gt. ftest .or. g .ge. gtest)) 
     +   task = 'WARNING: STP = STPMIN'

c     Test for convergence.

      if (f .le. ftest .and. abs(g) .le. gtol*(-ginit)) 
     +   task = 'CONVERGENCE'

c     Test for termination.

      if (task(1:4) .eq. 'WARN' .or. task(1:4) .eq. 'CONV') goto 1000

c     A modified function is used to predict the step during the
c     first stage if a lower function value has been obtained but 
c     the decrease is not sufficient.

      if (stage .eq. 1 .and. f .le. fx .and. f .gt. ftest) then

c        Define the modified function and derivative values.

         fm = f - stp*gtest
         fxm = fx - stx*gtest
         fym = fy - sty*gtest
         gm = g - gtest
         gxm = gx - gtest
         gym = gy - gtest

c        Call dcstep to update stx, sty, and to compute the new step.

         call dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,
     +               brackt,stmin,stmax)

c        Reset the function and derivative values for f.

         fx = fxm + stx*gtest
         fy = fym + sty*gtest
         gx = gxm + gtest
         gy = gym + gtest

      else

c       Call dcstep to update stx, sty, and to compute the new step.

        call dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,
     +              brackt,stmin,stmax)

      endif

c     Decide if a bisection step is needed.

      if (brackt) then
         if (abs(sty-stx) .ge. p66*width1) stp = stx + p5*(sty - stx)
         width1 = width
         width = abs(sty-stx)
      endif

c     Set the minimum and maximum steps allowed for stp.

      if (brackt) then
         stmin = min(stx,sty)
         stmax = max(stx,sty)
      else
         stmin = stp + xtrapl*(stp - stx)
         stmax = stp + xtrapu*(stp - stx)
      endif
 
c     Force the step to be within the bounds stpmax and stpmin.
 
      stp = max(stp,stpmin)
      stp = min(stp,stpmax)

c     If further progress is not possible, let stp be the best
c     point obtained during the search.

      if (brackt .and. (stp .le. stmin .or. stp .ge. stmax)
     +   .or. (brackt .and. stmax-stmin .le. xtol*stmax)) stp = stx

c     Obtain another function and derivative.

      task = 'FG'

 1000 continue

c     Save local variables.

      if (brackt) then
         isave(1) = 1
      else
         isave(1) = 0
      endif
      isave(2) = stage
      dsave(1) =  ginit
      dsave(2) =  gtest
      dsave(3) =  gx
      dsave(4) =  gy
      dsave(5) =  finit
      dsave(6) =  fx
      dsave(7) =  fy
      dsave(8) =  stx
      dsave(9) =  sty
      dsave(10) = stmin
      dsave(11) = stmax
      dsave(12) = width
      dsave(13) = width1

      return
      end
