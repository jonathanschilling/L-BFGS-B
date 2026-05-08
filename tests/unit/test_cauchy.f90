program test_cauchy
   use test_assert
   implicit none

   ! cauchy computes the generalised Cauchy point: the minimiser of the
   ! quadratic model along the projected steepest-descent path. Its body
   ! has many branches; covered here by hand-crafted small cases.

   call case_subgnorm_zero()
   call case_all_fixed_early_return()
   call case_unbounded_no_breakpoints()
   call case_one_breakpoint_col_zero()
   call case_with_history_col_one()

   write(*, '("test_cauchy: PASS")')

contains

   subroutine case_subgnorm_zero()
      ! sbgnrm <= 0 -> early return: xcp = x.
      integer, parameter :: n = 2, m = 1, col = 0
      integer  :: nbd(n), iorder(n), iwhere(n), nseg, info, head
      real(dp) :: x(n), l(n), u(n), g(n), t(n), d(n), xcp(n)
      real(dp) :: wy(n,1), ws(n,1), sy(m,m), wt(m,m)
      real(dp) :: p(2), c(2), wbp(2), v(2)
      real(dp) :: theta, sbgnrm, epsmch
      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 1.5_dp, -2.5_dp /)
      g = 0.0_dp
      iwhere = -1
      iorder = 0
      t = 0.0_dp; d = 0.0_dp; xcp = -42.0_dp
      wy = 0.0_dp; ws = 0.0_dp; sy = 0.0_dp; wt = 0.0_dp
      p = 0.0_dp; c = 0.0_dp; wbp = 0.0_dp; v = 0.0_dp
      theta = 1.0_dp
      sbgnrm = 0.0_dp
      epsmch = 2.22e-16_dp
      info = 0; head = 1; nseg = 0
      call cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, &
                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                  v, nseg, -1, sbgnrm, info, epsmch)
      call assert_close_real(xcp(1),  1.5_dp, where="case_subgnorm_zero xcp(1)")
      call assert_close_real(xcp(2), -2.5_dp, where="case_subgnorm_zero xcp(2)")
      call assert_eq_int(info, 0, where="case_subgnorm_zero info")
   end subroutine case_subgnorm_zero

   subroutine case_all_fixed_early_return()
      ! All variables marked fixed (iwhere=3): d stays 0, nbreak=0,
      ! nfree=n+1 -> the second early-return path fires.
      integer, parameter :: n = 2, m = 1, col = 0
      integer  :: nbd(n), iorder(n), iwhere(n), nseg, info, head
      real(dp) :: x(n), l(n), u(n), g(n), t(n), d(n), xcp(n)
      real(dp) :: wy(n,1), ws(n,1), sy(m,m), wt(m,m)
      real(dp) :: p(2), c(2), wbp(2), v(2)
      real(dp) :: theta, sbgnrm, epsmch
      nbd = (/ 2, 2 /)
      l = (/ 1.0_dp, 1.0_dp /); u = (/ 1.0_dp, 1.0_dp /)
      x = (/ 1.0_dp, 1.0_dp /)
      g = (/ -1.0_dp, -1.0_dp /)
      iwhere = (/ 3, 3 /)
      iorder = 0
      t = 0.0_dp; d = 0.0_dp; xcp = -42.0_dp
      wy = 0.0_dp; ws = 0.0_dp; sy = 0.0_dp; wt = 0.0_dp
      p = 0.0_dp; c = 0.0_dp; wbp = 0.0_dp; v = 0.0_dp
      theta = 1.0_dp
      sbgnrm = 1.0_dp
      epsmch = 2.22e-16_dp
      info = 0; head = 1; nseg = 0
      call cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, &
                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                  v, nseg, -1, sbgnrm, info, epsmch)
      ! xcp should equal x (no movement; all fixed).
      call assert_close_real(xcp(1), 1.0_dp, where="case_all_fixed_early xcp(1)")
      call assert_close_real(xcp(2), 1.0_dp, where="case_all_fixed_early xcp(2)")
      call assert_eq_int(info, 0, where="case_all_fixed_early info")
   end subroutine case_all_fixed_early_return

   subroutine case_unbounded_no_breakpoints()
      ! All variables unbounded (nbd=0, iwhere=-1): nbreak=0 but nfree=1<n+1=3.
      ! After the per-i loop we have d = -g. col=0, so f1=-||g||^2, f2=-theta*f1=theta*||g||^2.
      ! dtm = 1/theta. nbreak=0 -> goto 888. xcp = x + dtm*d = x - g/theta.
      integer, parameter :: n = 2, m = 1, col = 0
      integer  :: nbd(n), iorder(n), iwhere(n), nseg, info, head
      real(dp) :: x(n), l(n), u(n), g(n), t(n), d(n), xcp(n)
      real(dp) :: wy(n,1), ws(n,1), sy(m,m), wt(m,m)
      real(dp) :: p(2), c(2), wbp(2), v(2)
      real(dp) :: theta, sbgnrm, epsmch
      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 0.0_dp, 0.0_dp /)
      g = (/ 1.0_dp, 2.0_dp /)
      iwhere = -1
      iorder = 0
      t = 0.0_dp; d = 0.0_dp; xcp = -42.0_dp
      wy = 0.0_dp; ws = 0.0_dp; sy = 0.0_dp; wt = 0.0_dp
      p = 0.0_dp; c = 0.0_dp; wbp = 0.0_dp; v = 0.0_dp
      theta = 1.0_dp
      sbgnrm = sqrt(5.0_dp)
      epsmch = 2.22e-16_dp
      info = 0; head = 1; nseg = 0
      call cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, &
                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                  v, nseg, -1, sbgnrm, info, epsmch)
      ! Expected: xcp = -g/theta = (-1, -2).
      call assert_close_real(xcp(1), -1.0_dp, where="case_unbounded xcp(1)")
      call assert_close_real(xcp(2), -2.0_dp, where="case_unbounded xcp(2)")
      call assert_eq_int(info, 0, where="case_unbounded info")
   end subroutine case_unbounded_no_breakpoints

   subroutine case_one_breakpoint_col_zero()
      ! n=2, both bounded, x interior, gradient pushes one variable into
      ! its bound first. Exercises one iteration of the breakpoint loop.
      ! Setup:
      !   x = (0.5, 0.5), l = (0, 0), u = (1, 1), nbd = (2, 2)
      !   g = (-1, -2), so -g = (1, 2): both push toward upper.
      ! Breakpoints:
      !   var 1: t1 = (u-x)/(-g) = 0.5 / 1 = 0.5
      !   var 2: t2 = 0.5 / 2 = 0.25       (smaller -> hits first)
      ! col=0: f1 = -||g||^2 = -5, f2 = theta*||g||^2 = 5. dtm = 1.
      ! Iter 1: tj=bkmin=0.25, dt=0.25. dtm(1) >= dt -> fix var 2:
      !   xcp(2) = u(2) = 1, iwhere(2) = 2. Update f1=-0.75, f2=1, dtm=0.75.
      ! Iter 2: tj=0.5, dt=0.25. dtm(0.75) >= dt -> fix var 1:
      !   xcp(1) = u(1) = 1. nleft=0 and nbreak=n -> dtm=dt, goto 999.
      ! Final xcp = (1, 1) (the upper-corner of the box).
      integer, parameter :: n = 2, m = 1, col = 0
      integer  :: nbd(n), iorder(n), iwhere(n), nseg, info, head
      real(dp) :: x(n), l(n), u(n), g(n), t(n), d(n), xcp(n)
      real(dp) :: wy(n,1), ws(n,1), sy(m,m), wt(m,m)
      real(dp) :: p(2), c(2), wbp(2), v(2)
      real(dp) :: theta, sbgnrm, epsmch
      nbd = (/ 2, 2 /)
      l = (/ 0.0_dp, 0.0_dp /); u = (/ 1.0_dp, 1.0_dp /)
      x = (/ 0.5_dp, 0.5_dp /)
      g = (/ -1.0_dp, -2.0_dp /)
      iwhere = (/ 0, 0 /)   ! free, will be classified inside cauchy
      iorder = 0
      t = 0.0_dp; d = 0.0_dp; xcp = -42.0_dp
      wy = 0.0_dp; ws = 0.0_dp; sy = 0.0_dp; wt = 0.0_dp
      p = 0.0_dp; c = 0.0_dp; wbp = 0.0_dp; v = 0.0_dp
      theta = 1.0_dp
      sbgnrm = 2.0_dp
      epsmch = 2.22e-16_dp
      info = 0; head = 1; nseg = 0
      call cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, &
                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                  v, nseg, -1, sbgnrm, info, epsmch)
      call assert_close_real(xcp(1), 1.0_dp, where="case_one_break xcp(1)")
      call assert_close_real(xcp(2), 1.0_dp, where="case_one_break xcp(2)")
      call assert_eq_int(iwhere(1), 2, where="case_one_break iwhere(1)")
      call assert_eq_int(iwhere(2), 2, where="case_one_break iwhere(2)")
      call assert_eq_int(info, 0, where="case_one_break info")
   end subroutine case_one_breakpoint_col_zero

   subroutine case_with_history_col_one()
      ! Exercises the col>0 branch (bmv invocation, p/c updates).
      ! Use the same n=2 problem as above with one stored (s, y) pair.
      integer, parameter :: n = 2, m = 1, col = 1
      integer  :: nbd(n), iorder(n), iwhere(n), nseg, info, head
      real(dp) :: x(n), l(n), u(n), g(n), t(n), d(n), xcp(n)
      real(dp) :: wy(n,col), ws(n,col), sy(m,m), wt(m,m)
      real(dp) :: p(2*col), c(2*col), wbp(2*col), v(2*col)
      real(dp) :: theta, sbgnrm, epsmch
      nbd = (/ 2, 2 /)
      l = (/ 0.0_dp, 0.0_dp /); u = (/ 1.0_dp, 1.0_dp /)
      x = (/ 0.5_dp, 0.5_dp /)
      g = (/ -1.0_dp, -2.0_dp /)
      iwhere = (/ 0, 0 /)
      iorder = 0
      t = 0.0_dp; d = 0.0_dp; xcp = -42.0_dp
      ws(:,1) = (/ 0.5_dp, 0.0_dp /)
      wy(:,1) = (/ 1.0_dp, 1.0_dp /)
      sy(1,1) = 0.5_dp                ! s'y > 0 (curvature)
      wt(1,1) = sqrt(0.25_dp)         ! sqrt(theta * s's) for col=1, L empty
      p = 0.0_dp; c = 0.0_dp; wbp = 0.0_dp; v = 0.0_dp
      theta = 1.0_dp
      sbgnrm = 2.0_dp
      epsmch = 2.22e-16_dp
      info = 0; head = 1; nseg = 0
      call cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, xcp, &
                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                  v, nseg, -1, sbgnrm, info, epsmch)
      ! Don't predict xcp exactly (depends on bmv-driven curvature info);
      ! instead verify info=0 and xcp is feasible.
      call assert_eq_int(info, 0, where="case_with_history info")
      call assert_true(xcp(1) >= l(1) .and. xcp(1) <= u(1), &
                       "case_with_history xcp(1) feasible")
      call assert_true(xcp(2) >= l(2) .and. xcp(2) <= u(2), &
                       "case_with_history xcp(2) feasible")
   end subroutine case_with_history_col_one

end program test_cauchy
