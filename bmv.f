c> \file bmv.f

c> \brief This subroutine computes the product of the 2m x 2m middle matrix 
c>        in the compact L-BFGS formula of B and a 2m vector v;  
c>        it returns the product in p.
c>
c> @param m On entry m is the maximum number of variable metric corrections
c>             used to define the limited memory matrix.<br/>
c>          On exit m is unchanged.
c>
c> @param sy On entry sy specifies the matrix S'Y.<br/>
c>           On exit sy is unchanged.
c>
c> @param wt On entry wt specifies the upper triangular matrix J' which is 
c>              the Cholesky factor of (thetaS'S+LD^(-1)L').<br/>
c>           On exit wt is unchanged.
c>
c> @param col On entry col specifies the number of s-vectors (or y-vectors)
c>               stored in the compact L-BFGS formula.<br/>
c>            On exit col is unchanged.
c>
c> @param v On entry v specifies vector v.<br/>
c>          On exit v is unchanged.
c>
c> @param p On entry p is unspecified.<br/>
c>          On exit p is the product Mv.
c>
c> @param info On entry info is unspecified.<br/>
c>             On exit info = <ul><li>0       for normal return,</li>
c>                                <li>nonzero for abnormal return when the system
c>                                    to be solved by dtrsl is singular.</li></ul>
      subroutine bmv(m, sy, wt, col, v, p, info)

      integer m, col, info
      double precision sy(m, m), wt(m, m), v(2*col), p(2*col)

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
 
      integer          i,k,i2
      double precision sum
 
      if (col .eq. 0) return
 
c     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
c                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

c       solve Jp2=v2+LD^(-1)v1.
      p(col + 1) = v(col + 1)
      do 20 i = 2, col
         i2 = col + i
         sum = 0.0d0
         do 10 k = 1, i - 1
            sum = sum + sy(i,k)*v(k)/sy(k,k)
  10     continue
         p(i2) = v(i2) + sum
  20  continue  
c     Solve the triangular system
      call dtrsl(wt,m,col,p(col+1),11,info)
      if (info .ne. 0) return
 
c       solve D^(1/2)p1=v1.
      do 30 i = 1, col
         p(i) = v(i)/sqrt(sy(i,i))
  30  continue 
 
c     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
c                    [  0         J'           ] [ p2 ]   [ p2 ]. 
 
c       solve J^Tp2=p2. 
      call dtrsl(wt,m,col,p(col+1),01,info)
      if (info .ne. 0) return
 
c       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
c                 =-D^(-1/2)p1+D^(-1)L'p2.  
      do 40 i = 1, col
         p(i) = -p(i)/sqrt(sy(i,i))
  40  continue
      do 60 i = 1, col
         sum = 0.d0
         do 50 k = i + 1, col
            sum = sum + sy(k,i)*p(col+k)/sy(i,i)
  50     continue
         p(i) = p(i) + sum
  60  continue

      return

      end
