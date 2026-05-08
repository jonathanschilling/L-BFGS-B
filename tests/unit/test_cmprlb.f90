program test_cmprlb
   use test_assert, only: dp, assert_true, assert_eq_int, assert_eq_str, assert_close_real, assert_array_close
   implicit none

   call case_unconstrained_with_history()
   call case_constrained_no_history()
   call case_constrained_zero_step()
   call case_constrained_nontrivial_step()

   write(*, '("test_cmprlb: PASS")')

contains

   subroutine case_unconstrained_with_history()
      ! cnstnd=false AND col>0 -> shortcut: r(i) = -g(i) for all i.
      ! No bmv call, no index lookup; r is filled across the full n-vector.
      integer, parameter :: n = 3, m = 2
      integer  :: index(n), col, head, nfree
      real(dp) :: ws(n,m), wy(n,m), sy(m,m), wt(m,m), wa(4*m)
      real(dp) :: x(n), g(n), z(n), r(n), theta
      logical  :: cnstnd
      ws = 0.0_dp; wy = 0.0_dp; sy = 0.0_dp; wt = 0.0_dp; wa = 0.0_dp
      x = 0.0_dp; z = 0.0_dp
      g = (/ 1.5_dp, -2.5_dp, 0.25_dp /)
      r = -999.0_dp
      theta = 1.0_dp; col = 1; head = 1; nfree = n
      index = (/ 1, 2, 3 /)
      cnstnd = .false.
      call cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, &
                  theta, col, head, nfree, cnstnd)
      call assert_close_real(r(1), -1.5_dp,  where="case_unconstrained_with_history r(1)")
      call assert_close_real(r(2), 2.5_dp,   where="case_unconstrained_with_history r(2)")
      call assert_close_real(r(3), -0.25_dp, where="case_unconstrained_with_history r(3)")
   end subroutine case_unconstrained_with_history

   subroutine case_constrained_no_history()
      ! cnstnd=true AND col=0 -> j-loop empty; r(i) = -theta*(z-x)(k) - g(k).
      integer, parameter :: n = 2, m = 1
      integer  :: index(n), col, head, nfree
      real(dp) :: ws(n,m), wy(n,m), sy(m,m), wt(m,m), wa(4*m)
      real(dp) :: x(n), g(n), z(n), r(n), theta
      logical  :: cnstnd
      ws = 0.0_dp; wy = 0.0_dp; sy = 0.0_dp; wt = 0.0_dp; wa = 0.0_dp
      x = 0.0_dp; z = (/ 2.0_dp, 1.0_dp /)
      g = (/ 1.0_dp, 1.0_dp /)
      r = -999.0_dp
      theta = 1.0_dp; col = 0; head = 1; nfree = 2
      index = (/ 1, 2 /)
      cnstnd = .true.
      call cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, &
                  theta, col, head, nfree, cnstnd)
      ! Expected: r(1) = -1*2 - 1 = -3,  r(2) = -1*1 - 1 = -2.
      call assert_close_real(r(1), -3.0_dp, where="case_constrained_no_history r(1)")
      call assert_close_real(r(2), -2.0_dp, where="case_constrained_no_history r(2)")
   end subroutine case_constrained_no_history

   subroutine case_constrained_zero_step()
      ! z = x (zero step) so wa(2m+1..) = 0; bmv returns zero;
      ! the j-loop adds zero contributions; r reduces to -g(k) for free vars.
      integer, parameter :: n = 2, m = 1
      integer  :: index(n), col, head, nfree
      real(dp) :: ws(n,m), wy(n,m), sy(m,m), wt(m,m), wa(4*m)
      real(dp) :: x(n), g(n), z(n), r(n), theta
      logical  :: cnstnd
      ws(:,1) = (/ 1.0_dp, 0.0_dp /)
      wy(:,1) = (/ 1.0_dp, 1.0_dp /)
      sy(1,1) = 1.0_dp; wt(1,1) = 1.0_dp
      wa = 0.0_dp
      x = (/ 0.5_dp, 0.5_dp /); z = x
      g = (/ 1.0_dp, 1.0_dp /)
      r = -999.0_dp
      theta = 1.0_dp; col = 1; head = 1; nfree = 2
      index = (/ 1, 2 /)
      cnstnd = .true.
      call cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, &
                  theta, col, head, nfree, cnstnd)
      call assert_close_real(r(1), -1.0_dp, where="case_constrained_zero_step r(1)")
      call assert_close_real(r(2), -1.0_dp, where="case_constrained_zero_step r(2)")
   end subroutine case_constrained_zero_step

   subroutine case_constrained_nontrivial_step()
      ! Hand-derived n=2, col=1 problem.
      ! s = (1, 0), y = (1, 1), theta = 1.
      ! sy(1,1) = s'y = 1; wt(1,1) = sqrt(theta * s's) = 1 (col=1, L empty).
      ! z - x = (2, 1).
      ! Caller pre-fills wa(2m+1..) with W'(z-x):
      !   wa(3) = Y'(z-x) = 1*2 + 1*1 = 3
      !   wa(4) = theta * S'(z-x) = 1*(1*2 + 0*1) = 2
      ! M for col=1: [[-D, L'], [L, theta*S'S]] = [[-1, 0], [0, 1]]; M^{-1} same.
      ! So bmv produces p = M^{-1} (3, 2)' = (-3, 2)'.
      ! W = [Y, theta*S] = [[1, 1], [1, 0]].
      ! B = theta*I - W M^{-1} W' = I - W * [[-1,0],[0,1]] * W'
      !   W M^{-1} W' = [[0, -1], [-1, -1]]
      !   B = [[1, 1], [1, 2]]
      ! r = -B(z-x) - g = -[[1,1],[1,2]]*(2,1) - (1,1) = -(3,4) - (1,1) = (-4, -5).
      integer, parameter :: n = 2, m = 1
      integer  :: index(n), col, head, nfree
      real(dp) :: ws(n,m), wy(n,m), sy(m,m), wt(m,m), wa(4*m)
      real(dp) :: x(n), g(n), z(n), r(n), theta
      logical  :: cnstnd
      theta = 1.0_dp; col = 1; head = 1; nfree = 2
      index = (/ 1, 2 /)
      cnstnd = .true.
      ws(:,1) = (/ 1.0_dp, 0.0_dp /)
      wy(:,1) = (/ 1.0_dp, 1.0_dp /)
      sy(1,1) = 1.0_dp; wt(1,1) = 1.0_dp
      x = (/ 0.0_dp, 0.0_dp /); z = (/ 2.0_dp, 1.0_dp /)
      g = (/ 1.0_dp, 1.0_dp /)
      r = -999.0_dp
      ! W'(z-x), evaluated after theta is set:
      wa = 0.0_dp
      wa(2*m + 1) = wy(1,1)*(z(1)-x(1)) + wy(2,1)*(z(2)-x(2))
      wa(2*m + 2) = theta*(ws(1,1)*(z(1)-x(1)) + ws(2,1)*(z(2)-x(2)))
      call cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, &
                  theta, col, head, nfree, cnstnd)
      call assert_close_real(r(1), -4.0_dp, where="case_constrained_nontrivial_step r(1)")
      call assert_close_real(r(2), -5.0_dp, where="case_constrained_nontrivial_step r(2)")
   end subroutine case_constrained_nontrivial_step

end program test_cmprlb
