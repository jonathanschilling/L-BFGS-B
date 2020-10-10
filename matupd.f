c> \file matupd.f

      subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail, 
     +                  iupdat, col, head, theta, rr, dr, stp, dtd)
 
      integer          n, m, itail, iupdat, col, head
      double precision theta, rr, dr, stp, dtd, d(n), r(n), 
     +                 ws(n, m), wy(n, m), sy(m, m), ss(m, m)

c     ************
c
c     Subroutine matupd
c
c       This subroutine updates matrices WS and WY, and forms the
c         middle matrix in B.
c
c     Subprograms called:
c
c       Linpack ... dcopy, ddot.
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
 
      integer          j,pointr
      double precision ddot
      double precision one
      parameter        (one=1.0d0)

c     Set pointers for matrices WS and WY.
 
      if (iupdat .le. m) then
         col = iupdat
         itail = mod(head+iupdat-2,m) + 1
      else
         itail = mod(itail,m) + 1
         head = mod(head,m) + 1
      endif
 
c     Update matrices WS and WY.

      call dcopy(n,d,1,ws(1,itail),1)
      call dcopy(n,r,1,wy(1,itail),1)
 
c     Set theta=yy/ys.
 
      theta = rr/dr
 
c     Form the middle matrix in B.
 
c        update the upper triangle of SS,
c                                         and the lower triangle of SY:
      if (iupdat .gt. m) then
c                              move old information
         do 50 j = 1, col - 1
            call dcopy(j,ss(2,j+1),1,ss(1,j),1)
            call dcopy(col-j,sy(j+1,j+1),1,sy(j,j),1)
  50     continue
      endif
c        add new information: the last row of SY
c                                             and the last column of SS:
      pointr = head
      do 51 j = 1, col - 1
         sy(col,j) = ddot(n,d,1,wy(1,pointr),1)
         ss(j,col) = ddot(n,ws(1,pointr),1,d,1)
         pointr = mod(pointr,m) + 1
  51  continue
      if (stp .eq. one) then
         ss(col,col) = dtd
      else
         ss(col,col) = stp*stp*dtd
      endif
      sy(col,col) = dr
 
      return

      end
