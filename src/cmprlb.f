c> \file cmprlb.f

c> \brief This subroutine computes r=-Z'B(xcp-xk)-Z'g by using 
c>        wa(2m+1)=W'(xcp-x) from subroutine cauchy.
c> 
c> This subroutine computes r=-Z'B(xcp-xk)-Z'g by using 
c> wa(2m+1)=W'(xcp-x) from subroutine cauchy.
c> 
c> @param n number of parameters
c> @param m history size of Hessian approximation
c> @param x position
c> @param g gradient
c> @param ws part of L-BFGS matrix
c> @param wy part of L-BFGS matrix
c> @param sy part of L-BFGS matrix
c> @param wt part of L-BFGS matrix
c> @param z The generalized Cauchy point xcp computed by cauchy.
c>          Used here as the linearisation point: r encodes -B(z-x) - g
c>          restricted to the free variables.
c> @param r On exit r(1:nfree) contains -theta*(z(k)-x(k)) - g(k) plus the
c>          W*M^{-1}*W' correction, where k = index(i). Caller uses this as
c>          the residual for the subspace minimisation problem.
c>          When cnstnd=.false. and col>0, the shortcut path sets
c>          r(1:n) = -g(:) without using the index array.
c> @param wa Length-4m workspace shared with cauchy. On entry, the segment
c>           wa(2m+1 : 2m+2col) holds W'(z-x) (filled by cauchy). On exit
c>           wa(1 : 2col) holds M^{-1}*W'(z-x) from the bmv call.
c> @param index Permutation of (1..n): index(1..nfree) lists the indices of
c>              variables that are free at the GCP and are the active
c>              optimisation variables here.
c> @param theta Scaling factor specifying the initial Hessian B_0 = theta*I.
c> @param col Number of stored (s,y) correction pairs (0 on the first
c>            iteration; up to m thereafter).
c> @param head Index in the cyclic WS/WY buffer of the oldest stored
c>             correction. Used to walk the columns in chronological order.
c> @param nfree Number of free variables; size of the subspace problem.
c> @param cnstnd .true. if the problem has bounds; controls the shortcut
c>               path described under @param r.
c>
c> Historical note: this routine used to take an `info` output parameter
c> to forward errors from the embedded `bmv` call. Since `bmv` cannot
c> fail under LAPACK `dtrsm`, the parameter was always 0 on exit and
c> has been removed.
      subroutine cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index,
     +                 theta, col, head, nfree, cnstnd)

      logical          cnstnd
      integer          n, m, col, head, nfree, index(n)
      double precision theta,
     +                 x(n), g(n), z(n), r(n), wa(4*m),
     +                 ws(n, m), wy(n, m), sy(m, m), wt(m, m)
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
         call bmv(m,sy,wt,col,wa(2*m+1),wa(1))
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
