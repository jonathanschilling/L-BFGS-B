c> \file formk.f

c> \brief Forms the LEL^T factorization of the indefinite matrix K.
c>
c> This subroutine forms the LEL^T factorization of the indefinite matrix
c>
c> K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c>     [L_a -R_z           theta*S'AA'S ]
c>                                    where E = [-I  0]
c>                                              [ 0  I]
c>
c> The matrix K can be shown to be equal to the matrix M^[-1]N
c>   occurring in section 5.1 of [1], as well as to the matrix
c>   Mbar^[-1] Nbar in section 5.3.
c>
c> @param n On entry n is the dimension of the problem.<br/>
c>          On exit n is unchanged.
c>
c> @param nsub On entry nsub is the number of subspace variables in free set.<br/>
c>             On exit nsub is not changed.
c>
c> @param ind On entry ind specifies the indices of subspace variables.<br/>
c>            On exit ind is unchanged.
c>
c> @param nenter On entry nenter is the number of variables entering the
c>                  free set.<br/>
c>               On exit nenter is unchanged.
c>
c> @param ileave On entry indx2(ileave),...,indx2(n) are the variables leaving
c>                  the free set.<br/>
c>               On exit ileave is unchanged.
c>
c> @param indx2 On entry indx2(1),...,indx2(nenter) are the variables entering
c>                 the free set, while indx2(ileave),...,indx2(n) are the
c>                 variables leaving the free set.<br/>
c>              On exit indx2 is unchanged.
c>
c> @param iupdat On entry iupdat is the total number of BFGS updates made so far.<br/>
c>               On exit iupdat is unchanged.
c>
c> @param updatd On entry 'updatd' is true if the L-BFGS matrix is updated.<br/>
c>               On exit 'updatd' is unchanged.
c>
c> @param wn On entry wn is unspecified.<br/>
c>           On exit the upper triangle of wn stores the LEL^T factorization
c>              of the 2*col x 2*col indefinite matrix
c>                         [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
c>                         [L_a -R_z           theta*S'AA'S ]
c>
c> @param wn1 On entry wn1 stores the lower triangular part of
c>                          [Y' ZZ'Y   L_a'+R_z']
c>                          [L_a+R_z   S'AA'S   ]
c>               in the previous iteration.<br/>
c>            On exit wn1 stores the corresponding updated matrices.<br/>
c>            The purpose of wn1 is just to store these inner products
c>            so they can be easily updated and inserted into wn.
c>
c> @param m On entry m is the maximum number of variable metric corrections
c>             used to define the limited memory matrix.<br/>
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
c> @param theta On entry theta is the scaling factor specifying B_0 = theta I.<br/>
c>              On exit theta is unchanged.
c>
c> @param col On entry col is the actual number of variable metric
c>               corrections stored so far.<br/>
c>            On exit col is unchanged.
c>
c> @param head On entry head is the location of the first s-vector (or y-vector)
c>                in S (or Y).<br/>
c>             On exit col is unchanged.
c>
c> @param info On entry info is unspecified.<br/>
c>             On exit info<ul><li>=  0 for normal return;</li>
c>                             <li>= -1 when the 1st Cholesky factorization failed;</li>
c>                             <li>= -2 when the 2st Cholesky factorization failed.</li></ul>
      subroutine formk(n, nsub, ind, nenter, ileave, indx2, iupdat,
     +                 updatd, wn, wn1, m, ws, wy, sy, theta, col,
     +                 head, info)

      integer          n, nsub, m, col, head, nenter, ileave, iupdat,
     +                 info, ind(n), indx2(n)
      double precision theta, wn(2*m, 2*m), wn1(2*m, 2*m),
     +                 ws(n, m), wy(n, m), sy(m, m)
      logical          updatd
c
c     References:
c       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c       memory algorithm for bound constrained optimization'',
c       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
c       limited memory FORTRAN code for solving bound constrained
c       optimization problems'', Tech. Report, NAM-11, EECS Department,
c       Northwestern University, 1994.
c
c       (Postscript files of these papers are available via anonymous
c        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
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

      integer          m2,ipntr,jpntr,iy,is,jy,js,is1,js1,k1,i,k,
     +                 col2,pbegin,pend,dbegin,dend,upcl
      double precision ddot,temp1,temp2,temp3,temp4
      double precision one,zero
      parameter        (one=1.0d0,zero=0.0d0)

c     Form the lower triangular part of
c               WN1 = [Y' ZZ'Y   L_a'+R_z']
c                     [L_a+R_z   S'AA'S   ]
c        where L_a is the strictly lower triangular part of S'AA'Y
c              R_z is the upper triangular part of S'ZZ'Y.

      if (updatd) then
         if (iupdat .gt. m) then
c                                 shift old part of WN1.
            do 10 jy = 1, m - 1
               js = m + jy
               call dcopy(m-jy,wn1(jy+1,jy+1),1,wn1(jy,jy),1)
               call dcopy(m-jy,wn1(js+1,js+1),1,wn1(js,js),1)
               call dcopy(m-1,wn1(m+2,jy+1),1,wn1(m+1,jy),1)
  10        continue
         endif

c          put new rows in blocks (1,1), (2,1) and (2,2).
         pbegin = 1
         pend = nsub
         dbegin = nsub + 1
         dend = n
         iy = col
         is = m + col
         ipntr = head + col - 1
         if (ipntr .gt. m) ipntr = ipntr - m
         jpntr = head
         do 20 jy = 1, col
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
c             compute element jy of row 'col' of Y'ZZ'Y
            do 15 k = pbegin, pend
               k1 = ind(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
  15        continue
c             compute elements jy of row 'col' of L_a and S'AA'S
            do 16 k = dbegin, dend
               k1 = ind(k)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  16        continue
            wn1(iy,jy) = temp1
            wn1(is,js) = temp2
            wn1(is,jy) = temp3
            jpntr = mod(jpntr,m) + 1
  20     continue

c          put new column in block (2,1).
         jy = col
         jpntr = head + col - 1
         if (jpntr .gt. m) jpntr = jpntr - m
         ipntr = head
         do 30 i = 1, col
            is = m + i
            temp3 = zero
c             compute element i of column 'col' of R_z
            do 25 k = pbegin, pend
               k1 = ind(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  25        continue
            ipntr = mod(ipntr,m) + 1
            wn1(is,jy) = temp3
  30     continue
         upcl = col - 1
      else
         upcl = col
      endif

c       modify the old parts in blocks (1,1) and (2,2) due to changes
c       in the set of free variables.
      ipntr = head
      do 45 iy = 1, upcl
         is = m + iy
         jpntr = head
         do 40 jy = 1, iy
            js = m + jy
            temp1 = zero
            temp2 = zero
            temp3 = zero
            temp4 = zero
            do 35 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr)
               temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr)
  35        continue
            do 36 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr)
               temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr)
  36        continue
            wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3
            wn1(is,js) = wn1(is,js) - temp2 + temp4
            jpntr = mod(jpntr,m) + 1
  40     continue
         ipntr = mod(ipntr,m) + 1
  45  continue

c       modify the old parts in block (2,1).
      ipntr = head
      do 60 is = m + 1, m + upcl
         jpntr = head
         do 55 jy = 1, upcl
            temp1 = zero
            temp3 = zero
            do 50 k = 1, nenter
               k1 = indx2(k)
               temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr)
  50        continue
            do 51 k = ileave, n
               k1 = indx2(k)
               temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr)
  51        continue
         if (is .le. jy + m) then
               wn1(is,jy) = wn1(is,jy) + temp1 - temp3
            else
               wn1(is,jy) = wn1(is,jy) - temp1 + temp3
            endif
            jpntr = mod(jpntr,m) + 1
  55     continue
         ipntr = mod(ipntr,m) + 1
  60  continue

c     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
c                                     [-L_a +R_z        S'AA'S*theta]

      m2 = 2*m
      do 70 iy = 1, col
         is = col + iy
         is1 = m + iy
         do 65 jy = 1, iy
            js = col + jy
            js1 = m + jy
            wn(jy,iy) = wn1(iy,jy)/theta
            wn(js,is) = wn1(is1,js1)*theta
  65     continue
         do 66 jy = 1, iy - 1
            wn(jy,is) = -wn1(is1,jy)
  66     continue
         do 67 jy = iy, col
            wn(jy,is) = wn1(is1,jy)
  67     continue
         wn(iy,iy) = wn(iy,iy) + sy(iy,iy)
  70  continue

c     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
c                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ]

c        first Cholesky factor (1,1) block of wn to get LL'
c                          with L' stored in the upper triangle of wn.
      !call dpofa(wn,m2,col,info)
      call dpotrf('U',col,wn,m2,info)

      if (info .ne. 0) then
         info = -1
         return
      endif
c        then form L^-1(-L_a'+R_z') in the (1,2) block.
      col2 = 2*col
      do 71 js = col+1 ,col2
         !call dtrsl(wn,m2,col,wn(1,js),11,info)

         ! TODO: can combine this loop into a single dtrsm call?
         call dtrsm('l','u','t','n',col,1,one,wn,m2,wn(1,js),col)
  71  continue

c     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
c        upper triangle of (2,2) block of wn.


      do 72 is = col+1, col2
         do 74 js = is, col2
               wn(is,js) = wn(is,js) + ddot(col,wn(1,is),1,wn(1,js),1)
  74        continue
  72     continue

c     Cholesky factorization of (2,2) block of wn.

      !call dpofa(wn(col+1,col+1),m2,col,info)
      call dpotrf('U',col,wn(col+1,col+1),m2,info)

      if (info .ne. 0) then
         info = -2
         return
      endif

      return

      end
