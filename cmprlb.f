c> \file cmprlb.f

      subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, 
     +                 theta, col, head, nfree, cnstnd, info)
 
      logical          cnstnd
      integer          n, m, col, head, nfree, info, index(n)
      double precision theta, 
     +                 x(n), g(n), z(n), r(n), wa(4*m), 
     +                 ws(n, m), wy(n, m), sy(m, m), wt(m, m)

c     ************
c
c     Subroutine cmprlb 
c
c       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using 
c         wa(2m+1)=W'(xcp-x) from subroutine cauchy.
c
c     Subprograms called:
c
c       L-BFGS-B Library ... bmv.
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
 
      integer          i,j,k,pointr
      double precision a1,a2

      if (.not. cnstnd .and. col .gt. 0) then 
         do 26 i = 1, n
            r(i) = -g(i)
  26     continue
      else
         do 30 i = 1, nfree
            k = index(i)
            r(i) = -theta*(z(k) - x(k)) - g(k)
  30     continue
         call bmv(m,sy,wt,col,wa(2*m+1),wa(1),info)
         if (info .ne. 0) then
            info = -8
            return
         endif
         pointr = head 
         do 34 j = 1, col
            a1 = wa(j)
            a2 = theta*wa(col + j)
            do 32 i = 1, nfree
               k = index(i)
               r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2
  32        continue
            pointr = mod(pointr,m) + 1
  34     continue
      endif

      return

      end
