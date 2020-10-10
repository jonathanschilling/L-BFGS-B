c> \file prn3lb.f

      subroutine prn3lb(n, x, f, task, iprint, info, itfile, 
     +                  iter, nfgv, nintol, nskip, nact, sbgnrm, 
     +                  time, nseg, word, iback, stp, xstep, k, 
     +                  cachyt, sbtime, lnscht)
 
      character*60     task
      character*3      word
      integer          n, iprint, info, itfile, iter, nfgv, nintol,
     +                 nskip, nact, nseg, iback, k
      double precision f, sbgnrm, time, stp, xstep, cachyt, sbtime,
     +                 lnscht, x(n)

c     ************
c
c     Subroutine prn3lb
c
c     This subroutine prints out information when either a built-in
c       convergence test is satisfied or when an error message is
c       generated.
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

      if (task(1:5) .eq. 'ERROR') goto 999

      if (iprint .ge. 0) then
         write (6,3003)
         write (6,3004)
         write(6,3005) n,iter,nfgv,nintol,nskip,nact,sbgnrm,f
         if (iprint .ge. 100) then
            write (6,1004) 'X =',(x(i),i = 1,n)
         endif  
         if (iprint .ge. 1) write (6,*) ' F =',f
      endif 
 999  continue
      if (iprint .ge. 0) then
         write (6,3009) task
         if (info .ne. 0) then
            if (info .eq. -1) write (6,9011)
            if (info .eq. -2) write (6,9012)
            if (info .eq. -3) write (6,9013)
            if (info .eq. -4) write (6,9014)
            if (info .eq. -5) write (6,9015)
            if (info .eq. -6) write (6,*)' Input nbd(',k,') is invalid.'
            if (info .eq. -7) 
     +      write (6,*)' l(',k,') > u(',k,').  No feasible solution.'
            if (info .eq. -8) write (6,9018)
            if (info .eq. -9) write (6,9019)
         endif
         if (iprint .ge. 1) write (6,3007) cachyt,sbtime,lnscht
         write (6,3008) time
         if (iprint .ge. 1) then
            if (info .eq. -4 .or. info .eq. -9) then
               write (itfile,3002)
     +             iter,nfgv,nseg,nact,word,iback,stp,xstep
            endif
            write (itfile,3009) task
            if (info .ne. 0) then
               if (info .eq. -1) write (itfile,9011)
               if (info .eq. -2) write (itfile,9012)
               if (info .eq. -3) write (itfile,9013)
               if (info .eq. -4) write (itfile,9014)
               if (info .eq. -5) write (itfile,9015)
               if (info .eq. -8) write (itfile,9018)
               if (info .eq. -9) write (itfile,9019)
            endif
            write (itfile,3008) time
         endif
      endif

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 3002 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6x,'-',10x,'-')
 3003 format (/,
     + '           * * *',/,/,
     + 'Tit   = total number of iterations',/,
     + 'Tnf   = total number of function evaluations',/,
     + 'Tnint = total number of segments explored during',
     +           ' Cauchy searches',/,
     + 'Skip  = number of BFGS updates skipped',/,
     + 'Nact  = number of active bounds at final generalized',
     +          ' Cauchy point',/,
     + 'Projg = norm of the final projected gradient',/,
     + 'F     = final function value',/,/,
     + '           * * *')
 3004 format (/,3x,'N',4x,'Tit',5x,'Tnf',2x,'Tnint',2x,
     +       'Skip',2x,'Nact',5x,'Projg',8x,'F')
 3005 format (i5,2(1x,i6),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d10.3))
 3007 format (/,' Cauchy                time',1p,e10.3,' seconds.',/ 
     +        ' Subspace minimization time',1p,e10.3,' seconds.',/
     +        ' Line search           time',1p,e10.3,' seconds.')
 3008 format (/,' Total User time',1p,e10.3,' seconds.',/)
 3009 format (/,a60)
 9011 format (/,
     +' Matrix in 1st Cholesky factorization in formk is not Pos. Def.')
 9012 format (/,
     +' Matrix in 2st Cholesky factorization in formk is not Pos. Def.')
 9013 format (/,
     +' Matrix in the Cholesky factorization in formt is not Pos. Def.')
 9014 format (/,
     +' Derivative >= 0, backtracking line search impossible.',/,
     +'   Previous x, f and g restored.',/,
     +' Possible causes: 1 error in function or gradient evaluation;',/,
     +'                  2 rounding errors dominate computation.')
 9015 format (/,
     +' Warning:  more than 10 function and gradient',/,
     +'   evaluations in the last line search.  Termination',/,
     +'   may possibly be caused by a bad search direction.')
 9018 format (/,' The triangular system is singular.')
 9019 format (/,
     +' Line search cannot locate an adequate point after 20 function',/
     +,'  and gradient evaluations.  Previous x, f and g restored.',/,
     +' Possible causes: 1 error in function or gradient evaluation;',/,
     +'                  2 rounding error dominate computation.')

      return

      end
