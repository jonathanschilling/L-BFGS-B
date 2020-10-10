c> \file errclb.f

      subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
 
      character*60     task
      integer          n, m, info, k, nbd(n)
      double precision factr, l(n), u(n)

c     ************
c
c     Subroutine errclb
c
c     This subroutine checks the validity of the input data.
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

      integer          i
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

c     Check the input arguments for errors.

      if (n .le. 0) task = 'ERROR: N .LE. 0'
      if (m .le. 0) task = 'ERROR: M .LE. 0'
      if (factr .lt. zero) task = 'ERROR: FACTR .LT. 0'

c     Check the validity of the arrays nbd(i), u(i), and l(i).

      do 10 i = 1, n
         if (nbd(i) .lt. 0 .or. nbd(i) .gt. 3) then
c                                                   return
            task = 'ERROR: INVALID NBD'
            info = -6
            k = i
         endif
         if (nbd(i) .eq. 2) then
            if (l(i) .gt. u(i)) then
c                                    return
               task = 'ERROR: NO FEASIBLE SOLUTION'
               info = -7
               k = i
            endif
         endif
  10  continue

      return

      end
