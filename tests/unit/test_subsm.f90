program test_subsm
   use test_assert
   implicit none

   call case_nsub_zero()
   call case_unconstrained_no_clip()
   call case_constrained_clip_upper()

   write(*, '("test_subsm: PASS")')

contains

   subroutine case_nsub_zero()
      ! nsub <= 0 -> immediate return.
      integer, parameter :: n = 1, m = 1, nsub = 0
      integer  :: nbd(n), ind(max(nsub,1))
      integer  :: col, head, iword, info
      real(dp) :: l(n), u(n), x(n), d(max(nsub,1)), xp(n), xx(n), gg(n)
      real(dp) :: ws(n,m), wy(n,m), wn(2*m, 2*m)
      real(dp) :: wv(2*m), theta
      nbd = 0; ind = 0
      l = 0.0_dp; u = 0.0_dp; x = (/ 1.5_dp /)
      d = 0.0_dp; xp = 0.0_dp; xx = 0.0_dp; gg = 0.0_dp
      ws = 0.0_dp; wy = 0.0_dp; wn = 0.0_dp; wv = 0.0_dp
      theta = 1.0_dp; col = 1; head = 1; iword = -1; info = 0
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1, info)
      ! No assertions on x (nsub=0 path doesn't touch it); just verify info=0.
      call assert_eq_int(info, 0, where="case_nsub_zero info")
   end subroutine case_nsub_zero

   subroutine case_unconstrained_no_clip()
      ! All variables free (nbd=0): the projection loop just does x(k) += d(i).
      ! iword stays 0 -> early return via 911.
      integer, parameter :: n = 2, m = 1, nsub = 2
      integer  :: nbd(n), ind(nsub)
      integer  :: col, head, iword, info
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
      theta = 1.0_dp; col = 1; head = 1; iword = -1; info = -42
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1, info)
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
      integer  :: col, head, iword, info
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
      theta = 1.0_dp; col = 1; head = 1; iword = -1; info = -42
      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, &
                 ws, wy, theta, xx, gg, col, head, iword, wv, wn, -1, info)
      call assert_close_real(x(1), 1.0_dp, where="case_clip_upper x(1)")
      call assert_eq_int(iword, 1, where="case_clip_upper iword")
   end subroutine case_constrained_clip_upper

end program test_subsm
