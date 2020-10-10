c> \file active.f

c> \brief This subroutine initializes iwhere and projects the initial x to
c>        the feasible set if necessary.
c> @param n number of parameters
c> @param l lower bounds on parameters
c> @param u upper bounds on parameters
c> @param nbd indicates which bounds are present
c> @param x position
c> @param iwhere
c>      On entry iwhere is unspecified.<br/>
c>      On exit: iwhere(i)=<ul><li>-1  if x(i) has no bounds</li>
c>                             <li> 3   if l(i)=u(i),</li>
c>                             <li> 0   otherwise.</li></ul>
c>      In cauchy, iwhere is given finer gradations.
c> @param iprint console output flag
c> @param prjctd TODO
c> @param cnstnd TODO
c> @param boxed TODO
c 
c                           *  *  *
c 
c     NEOS, November 1994. (Latest revision June 1996.)
c     Optimization Technology Center.
c     Argonne National Laboratory and Northwestern University.
c     Written by
c                        Ciyou Zhu
c     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
      subroutine active(n, l, u, nbd, x, iwhere, iprint,
     +                  prjctd, cnstnd, boxed)

      logical          prjctd, cnstnd, boxed
      integer          n, iprint, nbd(n), iwhere(n)
      double precision x(n), l(n), u(n)

c     ************

      integer          nbdd,i
      double precision zero
      parameter        (zero=0.0d0)

c     Initialize nbdd, prjctd, cnstnd and boxed.

      nbdd = 0
      prjctd = .false.
      cnstnd = .false.
      boxed = .true.

c     Project the initial x to the easible set if necessary.

      do 10 i = 1, n
         if (nbd(i) .gt. 0) then
            if (nbd(i) .le. 2 .and. x(i) .le. l(i)) then
               if (x(i) .lt. l(i)) then
                  prjctd = .true.
                  x(i) = l(i)
               endif
               nbdd = nbdd + 1
            else if (nbd(i) .ge. 2 .and. x(i) .ge. u(i)) then
               if (x(i) .gt. u(i)) then
                  prjctd = .true.
                  x(i) = u(i)
               endif
               nbdd = nbdd + 1
            endif
         endif
  10  continue

c     Initialize iwhere and assign values to cnstnd and boxed.

      do 20 i = 1, n
         if (nbd(i) .ne. 2) boxed = .false.
         if (nbd(i) .eq. 0) then
c                                this variable is always free
            iwhere(i) = -1

c           otherwise set x(i)=mid(x(i), u(i), l(i)).
         else
            cnstnd = .true.
            if (nbd(i) .eq. 2 .and. u(i) - l(i) .le. zero) then
c                   this variable is always fixed
               iwhere(i) = 3
            else 
               iwhere(i) = 0
            endif
         endif
  20  continue

      if (iprint .ge. 0) then
         if (prjctd) write (6,*)
     +   'The initial X is infeasible.  Restart with its projection.'
         if (.not. cnstnd)
     +      write (6,*) 'This problem is unconstrained.'
      endif

      if (iprint .gt. 0) write (6,1001) nbdd

 1001 format (/,'At X0 ',i9,' variables are exactly at the bounds') 

      return

      end
