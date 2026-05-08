program test_bmv
   use test_assert
   implicit none

   call case_col_zero()
   call case_col_one_diagonal()
   call case_col_two_with_L()

   write(*, '("test_bmv: PASS")')

contains

   subroutine case_col_zero()
      ! col=0: bmv returns immediately. p is untouched (caller is responsible
      ! for not relying on its contents; here we just verify no crash).
      ! Allocate sy/wt with m=1 to give addressable storage.
      integer, parameter :: m = 1
      integer  :: col, info
      real(dp) :: sy(m,m), wt(m,m), v(2), p(2)
      sy(1,1) = -1.0_dp; wt(1,1) = -1.0_dp
      v = (/ 7.0_dp, 11.0_dp /)
      p = -999.0_dp
      col = 0; info = -42
      call bmv(m, sy, wt, col, v, p, info)
      ! info is unspecified after early return; bmv doesn't reset it before
      ! `return` on col=0. The caller (cmprlb) doesn't read it in that path.
   end subroutine case_col_zero

   subroutine case_col_one_diagonal()
      ! col=1, m=1. sy(1,1) = s'y > 0. wt(1,1) = sqrt(theta * S'S)
      ! (L is empty for col=1, so wt is just sqrt of the diagonal).
      ! For col=1, M = [[-D, L'], [L, theta*S'S]] reduces to
      !   M = [[-1, 0], [0, theta*s's]].
      ! With sy(1,1)=1, theta*s's=1 -> M = [[-1, 0], [0, 1]],  M^{-1} same.
      ! For v=(3, 2): p = M^{-1} v = (-3, 2).
      integer, parameter :: m = 1
      integer  :: col, info
      real(dp) :: sy(m,m), wt(m,m), v(2), p(2)
      sy(1,1) = 1.0_dp
      wt(1,1) = 1.0_dp
      v = (/ 3.0_dp, 2.0_dp /)
      p = -999.0_dp
      col = 1; info = -42
      call bmv(m, sy, wt, col, v, p, info)
      call assert_close_real(p(1), -3.0_dp, where="case_col_one_diagonal p(1)")
      call assert_close_real(p(2),  2.0_dp, where="case_col_one_diagonal p(2)")
      call assert_eq_int(info, 0, where="case_col_one_diagonal info")
   end subroutine case_col_one_diagonal

   subroutine case_col_two_with_L()
      ! col=2 case constructed to exercise the full block algorithm including
      ! the strict-lower L-driven cross terms.
      !
      ! Design:
      !   D    = diag(1, 2)              (s_i'y_i values)
      !   L    = [[0, 0], [3, 0]]        (strict lower, encoded in sy(2,1))
      !   theta*S'S = diag(2, 5)
      ! So sy must satisfy: sy(1,1)=1, sy(2,2)=2, sy(2,1)=3 (upper part unused).
      !
      ! The Cholesky target J'J = theta*S'S + L*D^{-1}*L'
      !   L*D^{-1} = L (since D^{-1}_{1,1}=1 multiplies the only nonzero col)
      !   L*D^{-1}*L' = [[0,0],[0,9]]
      !   J'J = [[2, 0], [0, 14]]
      ! J = upper Cholesky factor = [[sqrt(2), 0], [0, sqrt(14)]].
      !
      ! For v = (1, 1, 1, 1):
      !   block 2: p2 = (theta*S'S + L*D^{-1}*L')^{-1} (v2 + L*D^{-1}*v1)
      !          = diag(1/2, 1/14) * ((1,1) + (0,3))   = (1/2, 4/14) = (1/2, 2/7)
      !   block 1: p1 = -D^{-1}*v1 + D^{-1}*L'*p2
      !          = -(1, 1/2) + diag(1, 1/2) * [[0,3],[0,0]] * (1/2, 2/7)
      !          = -(1, 1/2) + (6/7, 0)
      !          = (-1/7, -1/2)
      ! Expected p = (-1/7, -1/2, 1/2, 2/7).
      integer, parameter :: m = 2
      integer  :: col, info
      real(dp) :: sy(m,m), wt(m,m), v(4), p(4)
      sy = 0.0_dp
      sy(1,1) = 1.0_dp; sy(2,2) = 2.0_dp; sy(2,1) = 3.0_dp
      wt = 0.0_dp
      wt(1,1) = sqrt(2.0_dp)
      wt(2,2) = sqrt(14.0_dp)
      v = (/ 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp /)
      p = -999.0_dp
      col = 2; info = -42
      call bmv(m, sy, wt, col, v, p, info)
      call assert_close_real(p(1), -1.0_dp/7.0_dp,  where="case_col_two p(1)")
      call assert_close_real(p(2), -1.0_dp/2.0_dp,  where="case_col_two p(2)")
      call assert_close_real(p(3),  1.0_dp/2.0_dp,  where="case_col_two p(3)")
      call assert_close_real(p(4),  2.0_dp/7.0_dp,  where="case_col_two p(4)")
      call assert_eq_int(info, 0, where="case_col_two info")
   end subroutine case_col_two_with_L

end program test_bmv
