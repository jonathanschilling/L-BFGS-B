c> \file dpofa.f

c> \brief factors a double precision symmetric positive definite matrix.
c>
c> dpofa is usually called by dpoco, but it can be called
c> directly with a saving in time if  rcond  is not needed.
c> (time for dpoco) = (1 + 18/n)*(time for dpofa) .
c>
c> @param a On entry, this is the symmetric matrix to be factored.
c>             Only the diagonal and upper triangle are used.<br/>
c>          On exit this is an upper triangular matrix  r  so that  a = trans(r)*r
c>             where  trans(r)  is the transpose.
c>             The strict lower triangle is unaltered.
c>             If  info .ne. 0 , the factorization is not complete.
c>
c> @param lda On entry, this is the leading dimension of the array  a .<br/>
c>            On exit, this value is unaltered.
c> @param n On entry, this is the order of the matrix  a .<br/>
c>          On exit, this value is unaltered.
c> @param info On exit, this signals success:
c>             <ul>
c>             <li>= 0  for normal return.</li>
c>             <li>= k  signals an error condition.  the leading minor
c>                 of order  k  is not positive definite.</li>
c>             </ul>
      subroutine dpofa(a,lda,n,info)
      integer lda,n,info
      double precision a(lda,*)
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas ddot
c     fortran sqrt
c
c     internal variables
c
      double precision ddot,t
      double precision s
      integer j,jm1,k
c     begin block with ...exits to 40
c
c
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
c     ......exit
            if (s .le. 0.0d0) go to 40
            a(j,j) = sqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
