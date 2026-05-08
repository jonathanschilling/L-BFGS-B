program test_formk
   use test_assert
   implicit none

   ! formk forms the LEL^T factorisation of the indefinite K matrix used in
   ! subspace minimisation. It has many branches; here we cover the main
   ! paths (no-update vs first-update, both Cholesky-success and the two
   ! Cholesky-failure exits).

   call case_no_update_simple()
   call case_first_update_no_overflow()
   call case_first_cholesky_fail()

   write(*, '("test_formk: PASS")')

contains

   subroutine case_no_update_simple()
      ! updatd = .false. -> skip the "form WN1 update" block. With no entering
      ! or leaving variables, the modify-old loops are also no-ops. wn1 is
      ! consumed as-is to form wn.
      integer, parameter :: n = 1, m = 1
      integer  :: nsub, nenter, ileave, iupdat, col, head, info
      integer  :: ind(n), indx2(n)
      real(dp) :: wn(2*m, 2*m), wn1(2*m, 2*m)
      real(dp) :: ws(n, m), wy(n, m), sy(m, m), theta
      logical  :: updatd
      ws(1,1) = 0.5_dp
      wy(1,1) = 1.0_dp
      sy(1,1) = 0.5_dp
      ! Pre-loaded wn1 (would normally come from the previous iteration).
      wn1 = 0.0_dp
      wn1(1,1) = 1.0_dp; wn1(2,1) = 0.5_dp; wn1(2,2) = 0.25_dp
      wn = 0.0_dp
      nsub = 1; nenter = 0; ileave = n + 1
      iupdat = 0; col = 1; head = 1; theta = 1.0_dp
      ind = (/ 1 /); indx2 = 0
      info = -42
      updatd = .false.
      call formk(n, nsub, ind, nenter, ileave, indx2, iupdat, &
                 updatd, wn, wn1, m, ws, wy, sy, theta, col, head, info)
      call assert_eq_int(info, 0, where="case_no_update_simple info")
   end subroutine case_no_update_simple

   subroutine case_first_update_no_overflow()
      ! updatd = .true., iupdat = 1 <= m: skip the "shift WN1" branch but
      ! enter the "fill new rows/columns of WN1" block.
      integer, parameter :: n = 1, m = 1
      integer  :: nsub, nenter, ileave, iupdat, col, head, info
      integer  :: ind(n), indx2(n)
      real(dp) :: wn(2*m, 2*m), wn1(2*m, 2*m)
      real(dp) :: ws(n, m), wy(n, m), sy(m, m), theta
      logical  :: updatd
      ws(1,1) = 0.5_dp
      wy(1,1) = 1.0_dp
      sy(1,1) = 0.5_dp
      wn1 = 0.0_dp; wn = 0.0_dp
      nsub = 1; nenter = 0; ileave = n + 1
      iupdat = 1; col = 1; head = 1; theta = 1.0_dp
      ind = (/ 1 /); indx2 = 0
      info = -42
      updatd = .true.
      call formk(n, nsub, ind, nenter, ileave, indx2, iupdat, &
                 updatd, wn, wn1, m, ws, wy, sy, theta, col, head, info)
      call assert_eq_int(info, 0, where="case_first_update_no_overflow info")
   end subroutine case_first_update_no_overflow

   subroutine case_first_cholesky_fail()
      ! Make the (1,1) block non-positive-definite so the first dpotrf fails;
      ! formk then sets info = -1 and returns.
      ! sy(1,1) is added to wn(1,1) as the diagonal contribution. Choose
      ! sy(1,1) very negative to overwhelm wy^2 / theta and ensure the
      ! Cholesky's leading minor is non-PD.
      integer, parameter :: n = 1, m = 1
      integer  :: nsub, nenter, ileave, iupdat, col, head, info
      integer  :: ind(n), indx2(n)
      real(dp) :: wn(2*m, 2*m), wn1(2*m, 2*m)
      real(dp) :: ws(n, m), wy(n, m), sy(m, m), theta
      logical  :: updatd
      ws(1,1) = 0.0_dp
      wy(1,1) = 1.0_dp
      sy(1,1) = -10.0_dp
      wn1 = 0.0_dp; wn = 0.0_dp
      nsub = 1; nenter = 0; ileave = n + 1
      iupdat = 1; col = 1; head = 1; theta = 1.0_dp
      ind = (/ 1 /); indx2 = 0
      info = 0
      updatd = .true.
      call formk(n, nsub, ind, nenter, ileave, indx2, iupdat, &
                 updatd, wn, wn1, m, ws, wy, sy, theta, col, head, info)
      call assert_eq_int(info, -1, where="case_first_cholesky_fail info")
   end subroutine case_first_cholesky_fail

end program test_formk
