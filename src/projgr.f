c> \file projgr.f

c> \brief This subroutine computes the infinity norm of the projected
c>        gradient.
c> 
c> This subroutine computes the infinity norm of the projected
c> gradient.
c>
c> @param n On entry n is the number of variables.<br/>
c>          On exit n is unchanged.
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
c> @param x On entry x is an approximation to the solution.<br/>
c>          On exit x is unchanged.
c>
c> @param g On entry g is the gradient.<br/>
c>          On exit g is unchanged.
c>
c> @param sbgnrm infinity norm of projected gradient
      subroutine projgr(n, l, u, nbd, x, g, sbgnrm)

      integer          n, nbd(n)
      double precision sbgnrm, x(n), l(n), u(n), g(n)

c     ************
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

      integer i
      double precision gi
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

      sbgnrm = zero
      do 15 i = 1, n
        gi = g(i)
        if (nbd(i) .ne. 0) then
           if (gi .lt. zero) then
              if (nbd(i) .ge. 2) gi = max((x(i)-u(i)),gi)
           else
              if (nbd(i) .le. 2) gi = min((x(i)-l(i)),gi)
           endif
        endif
        sbgnrm = max(sbgnrm,abs(gi))
  15  continue

      return

      end
