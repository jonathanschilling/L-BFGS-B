c> \file formt.f

      subroutine formt(m, wt, sy, ss, col, theta, info)
 
      integer          m, col, info
      double precision theta, wt(m, m), sy(m, m), ss(m, m)

c     ************
c
c     Subroutine formt
c
c       This subroutine forms the upper half of the pos. def. and symm.
c         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
c         of the array wt, and performs the Cholesky factorization of T
c         to produce J*J', with J' stored in the upper triangle of wt.
c
c     Subprograms called:
c
c       Linpack ... dpofa.
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

      integer          i,j,k,k1
      double precision ddum
      double precision zero
      parameter        (zero=0.0d0)


c     Form the upper half of  T = theta*SS + L*D^(-1)*L',
c        store T in the upper triangle of the array wt.
 
      do 52 j = 1, col
         wt(1,j) = theta*ss(1,j)
  52  continue
      do 55 i = 2, col
         do 54 j = i, col
            k1 = min(i,j) - 1
            ddum  = zero
            do 53 k = 1, k1
               ddum  = ddum + sy(i,k)*sy(j,k)/sy(k,k)
  53        continue
            wt(i,j) = ddum + theta*ss(i,j)
  54     continue
  55  continue
 
c     Cholesky factorize T to J*J' with 
c        J' stored in the upper triangle of wt.
 
      call dpofa(wt,m,col,info)
      if (info .ne. 0) then
         info = -3
      endif

      return

      end
