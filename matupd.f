c> \file matupd.f

c> \brief This subroutine updates matrices WS and WY, and forms the
c>        middle matrix in B.
c>
c> This subroutine updates matrices WS and WY, and forms the
c> middle matrix in B.
c> 
c> @param n On entry n is the number of variables.<br/>
c>          On exit n is unchanged.
c>
c> @param m On entry m is the maximum number of variable metric
c>             corrections allowed in the limited memory matrix.<br/>
c>          On exit m is unchanged.
c>
c> @param ws On entry this stores S, a set of s-vectors, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param wy On entry this stores Y, a set of y-vectors, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param sy On entry this stores S'Y, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c>
c> @param ss On entry this stores S'S, that defines the
c>              limited memory BFGS matrix.<br/>
c>           On exit this array is unchanged.
c> @param d TODO
c> @param r TODO
c> @param itail TODO
c> @param iupdat TODO
c> @param col On entry col is the actual number of variable metric
c>               corrections stored so far.<br/>
c>            On exit col is unchanged.
c>
c> @param head On entry head is the location of the first s-vector (or y-vector)
c>                in S (or Y).<br/>
c>             On exit col is unchanged.
c>
c> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
c>              On exit theta is unchanged.
c> @param rr TODO
c> @param dr TODO
c> @param stp TODO
c> @param dtd TODO
      subroutine matupd(n, m, ws, wy, sy, ss, d, r, itail, 
     +                  iupdat, col, head, theta, rr, dr, stp, dtd)
 
      integer          n, m, itail, iupdat, col, head
      double precision theta, rr, dr, stp, dtd, d(n), r(n), 
     +                 ws(n, m), wy(n, m), sy(m, m), ss(m, m)

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
