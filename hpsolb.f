c> \file hpsolb.f

c> This subroutine sorts out the least element of t, and puts the
c>   remaining elements of t in a heap.
c> 
c> @param n On entry n is the dimension of the arrays t and iorder.<br/>
c>          On exit n is unchanged.
c> 
c> @param t On entry t stores the elements to be sorted.<br/>
c>          On exit t(n) stores the least elements of t, and t(1) to t(n-1)
c>             stores the remaining elements in the form of a heap.
c> 
c> @param iorder On entry iorder(i) is the index of t(i).<br/>
c>               On exit iorder(i) is still the index of t(i), but iorder may be
c>                  permuted in accordance with t.
c> 
c> @param iheap On entry iheap should be set as follows:<ul>
c>                 <li>iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,</li>
c>                 <li>iheap .ne. 0 if otherwise.</li></ul>
c>              On exit iheap is unchanged.
      subroutine hpsolb(n, t, iorder, iheap)
      integer          iheap, n, iorder(n)
      double precision t(n)
c
c     References:
c       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
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
c     ************
  
      integer          i,j,k,indxin,indxou
      double precision ddum,out

      if (iheap .eq. 0) then

c        Rearrange the elements t(1) to t(n) to form a heap.

         do 20 k = 2, n
            ddum  = t(k)
            indxin = iorder(k)

c           Add ddum to the heap.
            i = k
   10       continue
            if (i.gt.1) then
               j = i/2
               if (ddum .lt. t(j)) then
                  t(i) = t(j)
                  iorder(i) = iorder(j)
                  i = j
                  goto 10 
               endif  
            endif  
            t(i) = ddum
            iorder(i) = indxin
   20    continue
      endif
 
c     Assign to 'out' the value of t(1), the least member of the heap,
c        and rearrange the remaining members to form a heap as
c        elements 1 to n-1 of t.
 
      if (n .gt. 1) then
         i = 1
         out = t(1)
         indxou = iorder(1)
         ddum  = t(n)
         indxin  = iorder(n)

c        Restore the heap 
   30    continue
         j = i+i
         if (j .le. n-1) then
            if (t(j+1) .lt. t(j)) j = j+1
            if (t(j) .lt. ddum ) then
               t(i) = t(j)
               iorder(i) = iorder(j)
               i = j
               goto 30
            endif 
         endif 
         t(i) = ddum
         iorder(i) = indxin
 
c     Put the least member in t(n). 

         t(n) = out
         iorder(n) = indxou
      endif 

      return

      end
