program test_subsm
   use test_assert, only: dp, assert_true, assert_eq_int, assert_eq_str, assert_close_real, assert_array_close
   implicit none

   call case_nsub_zero()
   call case_unconstrained_no_clip()
   call case_constrained_clip_upper()
   call case_lower_only_clip()
   call case_upper_only_clip()
   call case_positive_dir_deriv_backtrack()
   call case_backtrack_lower_clip()
   call case_backtrack_var_at_lower()
   call case_backtrack_var_at_upper()

   write(*, '("test_subsm: PASS")')

contains

   subroutine case_nsub_zero()
      ! nsub <= 0 -> immediate return.
      integer, parameter :: n = 1, m = 1, nsub = 0
      integer  :: nbd(n), ind(max(nsub,1))
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(max(nsub,1)), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = 0; ind = 0
      l = 0.0_dp; u = 0.0_dp; x = (/ 1.5_dp /)
      d = 0.0_dp; xp = 0.0_dp; xx = 0.0_dp; gg = 0.0_dp
      ws = 0.0_dp; wy = 0.0_dp; wn = 0.0_dp; wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
   end subroutine case_nsub_zero

   subroutine case_unconstrained_no_clip()
      ! All variables free (nbd=0): the projection loop just does x(k) += d(i).
      ! iword stays 0 -> early return via 911.
      integer, parameter :: n = 2, m = 1, nsub = 2
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = 0
      ind = (/ 1, 2 /)
      l = 0.0_dp; u = 0.0_dp
      x = (/ 0.0_dp, 0.0_dp /)
      ! d incoming is the residual r. After subsm transforms it,
      ! the new x is x + transformed_d. With wn = identity and ws/wy = 0,
      ! the K^{-1} wv computation doesn't add anything, so d unchanged
      ! until the (1/theta) scaling.
      d = (/ 1.0_dp, 2.0_dp /)
      xp = -42.0_dp; xx = 0.0_dp; gg = 0.0_dp
      ws = 0.0_dp; wy = 0.0_dp
      ! wn must be invertible (dtrsm with diagonal=0 would fail).
      wn = 0.0_dp
      wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
      ! With ws=wy=0, wv stays 0, so d is just (1/theta)*d_in = d_in.
      ! Then x(k) += d(i): x = (1, 2). iword stays 0.
      call assert_close_real(x(1), 1.0_dp, where="case_unconstrained x(1)")
      call assert_close_real(x(2), 2.0_dp, where="case_unconstrained x(2)")
      call assert_eq_int(iword, 0, where="case_unconstrained iword")
   end subroutine case_unconstrained_no_clip

   subroutine case_constrained_clip_upper()
      ! Bounded variable, step pushes x past upper bound -> clipped, iword=1.
      ! n=1, nbd=2 (both bounds), x=0.5, l=0, u=1, d (transformed)=2.
      ! Projection: x(1) = max(l, xk + dk) = max(0, 2.5) = 2.5; min(u, 2.5) = 1.
      ! x(1) = u(1) -> iword = 1.
      ! After projection, the dir-deriv check uses (x - xx)*gg. Set gg, xx to
      ! ensure dd_p < 0 so we go to 911 (no backtracking).
      integer, parameter :: n = 1, m = 1, nsub = 1
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = (/ 2 /)
      ind = (/ 1 /)
      l = (/ 0.0_dp /); u = (/ 1.0_dp /)
      x = (/ 0.5_dp /)
      d = (/ 2.0_dp /)
      xp = -42.0_dp
      xx = (/ 0.5_dp /)         ! reference x for dd_p
      gg = (/ -1.0_dp /)        ! reference g; dd_p = (x_new - 0.5)*(-1) negative
      ws = 0.0_dp; wy = 0.0_dp
      wn = 0.0_dp
      wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
      call assert_close_real(x(1), 1.0_dp, where="case_clip_upper x(1)")
      call assert_eq_int(iword, 1, where="case_clip_upper iword")
   end subroutine case_constrained_clip_upper

   subroutine case_lower_only_clip()
      ! nbd=1 (lower bound only): exercises L247-249.
      ! x=0.5, l=0, dk=-2 -> x_new = max(0, -1.5) = 0. Hits bound -> iword=1.
      integer, parameter :: n = 1, m = 1, nsub = 1
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = (/ 1 /)
      ind = (/ 1 /)
      l = (/ 0.0_dp /); u = 0.0_dp
      x = (/ 0.5_dp /)
      d = (/ -2.0_dp /)
      xp = -42.0_dp
      xx = (/ 0.5_dp /); gg = (/ 1.0_dp /)         ! dd_p = (0-0.5)*1 = -0.5 < 0
      ws = 0.0_dp; wy = 0.0_dp
      wn = 0.0_dp; wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
      call assert_close_real(x(1), 0.0_dp, where="case_lower_only_clip x(1)")
      call assert_eq_int(iword, 1, where="case_lower_only_clip iword")
   end subroutine case_lower_only_clip

   subroutine case_upper_only_clip()
      ! nbd=3 (upper bound only): exercises L258-260.
      ! x=0.5, u=1, dk=2 -> x_new = min(1, 2.5) = 1. iword=1.
      integer, parameter :: n = 1, m = 1, nsub = 1
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = (/ 3 /)
      ind = (/ 1 /)
      l = 0.0_dp; u = (/ 1.0_dp /)
      x = (/ 0.5_dp /)
      d = (/ 2.0_dp /)
      xp = -42.0_dp
      xx = (/ 0.5_dp /); gg = (/ -1.0_dp /)         ! dd_p = (1-0.5)*(-1) < 0
      ws = 0.0_dp; wy = 0.0_dp
      wn = 0.0_dp; wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
      call assert_close_real(x(1), 1.0_dp, where="case_upper_only_clip x(1)")
      call assert_eq_int(iword, 1, where="case_upper_only_clip iword")
   end subroutine case_upper_only_clip

   subroutine case_positive_dir_deriv_backtrack()
      ! Positive directional derivative dd_p > 0 -> backtracking branch.
      ! Construct: x bumped against upper, but gg points uphill so the
      ! projected step actually increased the objective. The line-search
      ! clipping loop (L290+) computes alpha from the tightest active
      ! bound, then x = xp + alpha*d. Diagnostic prints to unit 6 are
      ! expected.
      ! n=2 with nbd=2 both. d projects var 1 to upper bound; gg makes
      ! dd_p > 0 so backtracking kicks in.
      integer, parameter :: n = 2, m = 1, nsub = 2
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = (/ 2, 2 /)
      ind = (/ 1, 2 /)
      l = (/ 0.0_dp, 0.0_dp /); u = (/ 1.0_dp, 1.0_dp /)
      x = (/ 0.5_dp, 0.5_dp /)
      d = (/ 2.0_dp, 0.5_dp /)            ! var 1 will hit upper, var 2 won't
      xp = -42.0_dp                        ! will get filled with x on entry
      xx = (/ 0.5_dp, 0.5_dp /)
      ! gg crafted so dd_p = (x_new(1) - 0.5)*1 + (x_new(2) - 0.5)*1 > 0
      ! (var 1 moves to 1.0, var 2 moves to 1.0 after subspace step;
      ! gg=(+1, +1) makes total positive)
      gg = (/ 1.0_dp, 1.0_dp /)
      ws = 0.0_dp; wy = 0.0_dp
      wn = 0.0_dp; wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
      ! After backtracking, x must be feasible.
      call assert_true(x(1) >= l(1) .and. x(1) <= u(1), &
                       "case_backtrack x(1) feasible")
      call assert_true(x(2) >= l(2) .and. x(2) <= u(2), &
                       "case_backtrack x(2) feasible")
   end subroutine case_positive_dir_deriv_backtrack

   subroutine case_backtrack_lower_clip()
      ! Backtracking with d<0 -> exercises the dk<0 lower-bound path in
      ! the alpha-clipping loop (L297-302) and the dk<0 branch in the
      ! "alpha < 1" final adjustment (L325-327).
      ! n=1, nbd=2, x=0.5, d=-2, l=0, u=1.
      ! Initial projection clips x to 0; gg=(-1) gives dd_p>0 so we
      ! enter the backtrack. alpha = (l-xp)/dk = -0.5/-2 = 0.25.
      integer, parameter :: n = 1, m = 1, nsub = 1
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = (/ 2 /)
      ind = (/ 1 /)
      l = (/ 0.0_dp /); u = (/ 1.0_dp /)
      x = (/ 0.5_dp /)
      d = (/ -2.0_dp /)
      xp = -42.0_dp
      xx = (/ 0.5_dp /)
      gg = (/ -1.0_dp /)         ! dd_p = (0-0.5)*(-1) = 0.5 > 0
      ws = 0.0_dp; wy = 0.0_dp
      wn = 0.0_dp; wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
      ! Final x should be at lower bound after backtrack.
      call assert_close_real(x(1), 0.0_dp, where="case_backtrack_lower_clip x(1)")
   end subroutine case_backtrack_lower_clip

   subroutine case_backtrack_var_at_lower()
      ! Backtrack with one variable already AT its lower bound on entry.
      ! Hits L300: temp2 = l(k) - x(k) = 0 (>= 0) -> temp1 = 0 -> alpha = 0.
      ! Setup: var 1 starts at lower, d(1)<0 wants to push it past;
      !        var 2 contributes a positive dd_p so we enter the backtrack.
      integer, parameter :: n = 2, m = 1, nsub = 2
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = (/ 2, 2 /)
      ind = (/ 1, 2 /)
      l = (/ 0.0_dp, 0.0_dp /); u = (/ 1.0_dp, 1.0_dp /)
      x = (/ 0.0_dp, 0.5_dp /)             ! var 1 already at lower
      d = (/ -1.0_dp, 1.0_dp /)
      xp = -42.0_dp
      xx = (/ 0.0_dp, 0.5_dp /)
      gg = (/ 0.0_dp, 1.0_dp /)            ! dd_p = (1-0.5)*1 = 0.5 > 0
      ws = 0.0_dp; wy = 0.0_dp
      wn = 0.0_dp; wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
   end subroutine case_backtrack_var_at_lower

   subroutine case_backtrack_var_at_upper()
      ! Symmetric: var already at upper, d wants to push it past.
      ! Hits L307: temp2 = u(k) - x(k) = 0 (<= 0) -> temp1 = 0.
      integer, parameter :: n = 2, m = 1, nsub = 2
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword
      real(dp) :: l(n), u(n), x(n), d(nsub), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = (/ 2, 2 /)
      ind = (/ 1, 2 /)
      l = (/ 0.0_dp, 0.0_dp /); u = (/ 1.0_dp, 1.0_dp /)
      x = (/ 1.0_dp, 0.5_dp /)             ! var 1 already at upper
      d = (/ 1.0_dp, -1.0_dp /)
      xp = -42.0_dp
      xx = (/ 1.0_dp, 0.5_dp /)
      gg = (/ 0.0_dp, -1.0_dp /)           ! dd_p = (0-0.5)*(-1) = 0.5 > 0
      ws = 0.0_dp; wy = 0.0_dp
      wn = 0.0_dp; wn(1,1) = 1.0_dp; wn(2,2) = 1.0_dp
      wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1)
   end subroutine case_backtrack_var_at_upper

end program test_subsm
