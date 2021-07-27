c> \file errclb.f

c> \brief This subroutine checks the validity of the input data.
c>
c> This subroutine checks the validity of the input data.
c>
c> @param n number of parameters
c> @param m history size of approximated Hessian
c> @param factr convergence criterion on function value
c> @param l lower bounds for parameters
c> @param u upper bounds for parameters
c> @param nbd indicates which bounds are present
c> @param task if an error occurs, contains a human-readable error message
c> @param info =0 on success; =-6 if nbd(k) was invalid; =-7 if both limits are given but l(k) > u(k)
c> @param k index of last errournous parameter
      subroutine errclb(n, m, factr, l, u, nbd, task, info, k)
 
      character*60     task
      integer          n, m, info, k, nbd(n)
      double precision factr, l(n), u(n)
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
