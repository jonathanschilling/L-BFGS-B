c> \file prn1lb.f

c> \brief This subroutine prints the input data, initial point, upper and
c>        lower bounds of each variable, machine precision, as well as 
c>        the headings of the output.
c> 
c> This subroutine prints the input data, initial point, upper and
c>        lower bounds of each variable, machine precision, as well as 
c>        the headings of the output.
c>
c> @param n On entry n is the number of variables.<br/>
c>          On exit n is unchanged.
c>
c> @param m On entry m is the maximum number of variable metric
c>             corrections allowed in the limited memory matrix.<br/>
c>          On exit m is unchanged.
c>
c> @param l On entry l is the lower bound of x.<br/>
c>          On exit l is unchanged.
c>
c> @param u On entry u is the upper bound of x.<br/>
c>          On exit u is unchanged.
c>
c> @param x On entry x is an approximation to the solution.<br/>
c>          On exit x is the current approximation.
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
c> @param itfile unit number of iterate.dat file
c> @param epsmch machine precision epsilon
      subroutine prn1lb(n, m, l, u, x, iprint, itfile, epsmch)
 
      integer n, m, iprint, itfile
      double precision epsmch, x(n), l(n), u(n)

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
