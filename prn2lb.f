c> \file prn2lb.f

      subroutine prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, 
     +                  sbgnrm, nseg, word, iword, iback, stp, xstep)
 
      character*3      word
      integer          n, iprint, itfile, iter, nfgv, nact, nseg,
     +                 iword, iback
      double precision f, sbgnrm, stp, xstep, x(n), g(n)

c     ************
c
c     Subroutine prn2lb
c
c     This subroutine prints out new information after a successful
c       line search. 
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

      integer i,imod

c           'word' records the status of subspace solutions.
      if (iword .eq. 0) then
c                            the subspace minimization converged.
         word = 'con'
      else if (iword .eq. 1) then
c                          the subspace minimization stopped at a bound.
         word = 'bnd'
      else if (iword .eq. 5) then
c                             the truncated Newton step has been used.
         word = 'TNT'
      else
         word = '---'
      endif
      if (iprint .ge. 99) then
         write (6,*) 'LINE SEARCH',iback,' times; norm of step = ',xstep
         write (6,2001) iter,f,sbgnrm
         if (iprint .gt. 100) then      
            write (6,1004) 'X =',(x(i), i = 1, n)
            write (6,1004) 'G =',(g(i), i = 1, n)
         endif
      else if (iprint .gt. 0) then 
         imod = mod(iter,iprint)
         if (imod .eq. 0) write (6,2001) iter,f,sbgnrm
      endif
      if (iprint .ge. 1) write (itfile,3001)
     +          iter,nfgv,nseg,nact,word,iback,stp,xstep,sbgnrm,f

 1004 format (/,a4, 1p, 6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))
 2001 format
     +  (/,'At iterate',i5,4x,'f= ',1p,d12.5,4x,'|proj g|= ',1p,d12.5)
 3001 format(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1p,2(1x,d10.3))

      return

      end
