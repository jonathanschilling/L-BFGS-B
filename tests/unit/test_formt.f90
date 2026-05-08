program test_formt
   use test_assert
   implicit none

   call case_col_zero()
   call case_col_one_simple()
   call case_col_two_diagonal_T()
   call case_non_positive_definite()

   write(*, '("test_formt: PASS")')

contains

   subroutine case_col_zero()
      ! col=0: all fill loops empty, dpotrf with N=0 returns info=0.
      integer, parameter :: m = 1
      integer  :: col, info
      real(dp) :: wt(m,m), sy(m,m), ss(m,m), theta
      wt = -42.0_dp; sy = 0.0_dp; ss = 0.0_dp
      col = 0; info = -42; theta = 1.0_dp
      call formt(m, wt, sy, ss, col, theta, info)
      call assert_eq_int(info, 0, where="case_col_zero info")
   end subroutine case_col_zero

   subroutine case_col_one_simple()
      ! col=1: only loop 52 fires (j=1). T(1,1) = theta * ss(1,1) = 4.
      ! dpotrf factors -> wt(1,1) = sqrt(4) = 2.
      integer, parameter :: m = 1
      integer  :: col, info
      real(dp) :: wt(m,m), sy(m,m), ss(m,m), theta
      wt = -42.0_dp
      sy(1,1) = 1.0_dp
      ss(1,1) = 4.0_dp
      col = 1; info = -42; theta = 1.0_dp
      call formt(m, wt, sy, ss, col, theta, info)
      call assert_close_real(wt(1,1), 2.0_dp, where="case_col_one_simple wt(1,1)")
      call assert_eq_int(info, 0, where="case_col_one_simple info")
   end subroutine case_col_one_simple

   subroutine case_col_two_diagonal_T()
      ! col=2 case where T = theta*S'S + L D^{-1} L' is diagonal.
      ! sy = [[1, ?], [3, 2]] (sy(1,1)=1, sy(2,1)=3, sy(2,2)=2),
      ! ss = diag(2, 5), theta = 1.
      !   wt(1,1) before = 1 * ss(1,1) = 2
      !   wt(1,2) before = 1 * ss(1,2) = 0
      !   k1 = 1 for (i,j)=(2,2): ddum = sy(2,1)^2 / sy(1,1) = 9
      !                            wt(2,2) before = 9 + 1*5 = 14
      ! T = diag(2, 14). dpotrf gives J' = diag(sqrt(2), sqrt(14)) in upper.
      integer, parameter :: m = 2
      integer  :: col, info
      real(dp) :: wt(m,m), sy(m,m), ss(m,m), theta
      wt = -42.0_dp
      sy = 0.0_dp
      sy(1,1) = 1.0_dp; sy(2,2) = 2.0_dp; sy(2,1) = 3.0_dp
      ss = 0.0_dp
      ss(1,1) = 2.0_dp; ss(2,2) = 5.0_dp
      col = 2; info = -42; theta = 1.0_dp
      call formt(m, wt, sy, ss, col, theta, info)
      call assert_close_real(wt(1,1), sqrt(2.0_dp),  where="case_col_two_diagonal_T wt(1,1)")
      call assert_close_real(wt(1,2), 0.0_dp,        where="case_col_two_diagonal_T wt(1,2)")
      call assert_close_real(wt(2,2), sqrt(14.0_dp), where="case_col_two_diagonal_T wt(2,2)")
      call assert_eq_int(info, 0, where="case_col_two_diagonal_T info")
   end subroutine case_col_two_diagonal_T

   subroutine case_non_positive_definite()
      ! theta < 0 makes T(1,1) = theta * ss(1,1) negative, so dpotrf fails.
      ! formt then maps that to info = -3.
      integer, parameter :: m = 1
      integer  :: col, info
      real(dp) :: wt(m,m), sy(m,m), ss(m,m), theta
      wt = 0.0_dp
      sy(1,1) = 1.0_dp
      ss(1,1) = 1.0_dp
      col = 1; info = 0; theta = -1.0_dp
      call formt(m, wt, sy, ss, col, theta, info)
      call assert_eq_int(info, -3, where="case_non_positive_definite info")
   end subroutine case_non_positive_definite

end program test_formt
