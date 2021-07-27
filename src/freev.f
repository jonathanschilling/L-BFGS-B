c> \file freev.f

c> \brief This subroutine counts the entering and leaving variables when
c>        iter > 0, and finds the index set of free and active variables
c>        at the GCP.
c>
c> This subroutine counts the entering and leaving variables when
c> iter > 0, and finds the index set of free and active variables
c> at the GCP.
c> 
c> @param n number of parameters
c> @param nfree number of free parameters, i.e., those not at their bounds
c> @param index for i=1,...,nfree, index(i) are the indices of free variables<br/>
c>              for i=nfree+1,...,n, index(i) are the indices of bound variables<br/>
c>              On entry after the first iteration, index gives 
c>                the free variables at the previous iteration.<br/>
c>              On exit it gives the free variables based on the determination
c>                in cauchy using the array iwhere.
c> @param nenter TODO
c> @param ileave TODO
c> @param indx2 On entry indx2 is unspecified.<br/>
c>              On exit with iter>0, indx2 indicates which variables
c>                 have changed status since the previous iteration.<br/>
c>              For i= 1,...,nenter, indx2(i) have changed from bound to free.<br/>
c>              For i= ileave+1,...,n, indx2(i) have changed from free to bound.<br/>
c> @param iwhere TODO
c> @param wrk TODO
c> @param updatd TODO
c> @param cnstnd indicating whether bounds are present
c> @param iprint control screen output
c> @param iter TODO
      subroutine freev(n, nfree, index, nenter, ileave, indx2, 
     +                 iwhere, wrk, updatd, cnstnd, iprint, iter)

      integer n, nfree, nenter, ileave, iprint, iter, 
     +        index(n), indx2(n), iwhere(n)
      logical wrk, updatd, cnstnd
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
 
      integer iact,i,k

      nenter = 0
      ileave = n + 1
      if (iter .gt. 0 .and. cnstnd) then
c                           count the entering and leaving variables.
         do 20 i = 1, nfree
            k = index(i)

c            write(6,*) ' k  = index(i) ', k
c            write(6,*) ' index = ', i

            if (iwhere(k) .gt. 0) then
               ileave = ileave - 1
               indx2(ileave) = k
               if (iprint .ge. 100) write (6,*)
     +             'Variable ',k,' leaves the set of free variables'
            endif
  20     continue
         do 22 i = 1 + nfree, n
            k = index(i)
            if (iwhere(k) .le. 0) then
               nenter = nenter + 1
               indx2(nenter) = k
               if (iprint .ge. 100) write (6,*)
     +             'Variable ',k,' enters the set of free variables'
            endif
  22     continue
         if (iprint .ge. 99) write (6,*)
     +       n+1-ileave,' variables leave; ',nenter,' variables enter'
      endif
      wrk = (ileave .lt. n+1) .or. (nenter .gt. 0) .or. updatd
 
c     Find the index set of free and active variables at the GCP.
 
      nfree = 0 
      iact = n + 1
      do 24 i = 1, n
         if (iwhere(i) .le. 0) then
            nfree = nfree + 1
            index(nfree) = i
         else
            iact = iact - 1
            index(iact) = i
         endif
  24  continue
      if (iprint .ge. 99) write (6,*)
     +      nfree,' variables are free at GCP ',iter + 1  

      return

      end
