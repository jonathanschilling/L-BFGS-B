c> \file prn1lb.f

      subroutine prn1lb(n, m, l, u, x, iprint, itfile, epsmch)
 
      integer n, m, iprint, itfile
      double precision epsmch, x(n), l(n), u(n)

c     ************
c
c     Subroutine prn1lb
c
c     This subroutine prints the input data, initial point, upper and
c       lower bounds of each variable, machine precision, as well as 
c       the headings of the output.
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

      integer i

      if (iprint .ge. 0) then
         write (6,7001) epsmch
         write (6,*) 'N = ',n,'    M = ',m
         if (iprint .ge. 1) then
            write (itfile,2001) epsmch
            write (itfile,*)'N = ',n,'    M = ',m
            write (itfile,9001)
            if (iprint .gt. 100) then
               write (6,1004) 'L =',(l(i),i = 1,n)
               write (6,1004) 'X0 =',(x(i),i = 1,n)
               write (6,1004) 'U =',(u(i),i = 1,n)
            endif 
         endif
      endif 

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format ('RUNNING THE L-BFGS-B CODE',/,/,
     + 'it    = iteration number',/,
     + 'nf    = number of function evaluations',/,
     + 'nseg  = number of segments explored during the Cauchy search',/,
     + 'nact  = number of active bounds at the generalized Cauchy point'
     + ,/,
     + 'sub   = manner in which the subspace minimization terminated:'
     + ,/,'        con = converged, bnd = a bound was reached',/,
     + 'itls  = number of iterations performed in the line search',/,
     + 'stepl = step length used',/,
     + 'tstep = norm of the displacement (total step)',/,
     + 'projg = norm of the projected gradient',/,
     + 'f     = function value',/,/,
     + '           * * *',/,/,
     + 'Machine precision =',1p,d10.3)
 7001 format ('RUNNING THE L-BFGS-B CODE',/,/,
     + '           * * *',/,/,
     + 'Machine precision =',1p,d10.3)
 9001 format (/,3x,'it',3x,'nf',2x,'nseg',2x,'nact',2x,'sub',2x,'itls',
     +        2x,'stepl',4x,'tstep',5x,'projg',8x,'f')

      return

      end
