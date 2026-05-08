program test_lnsrlb
   use test_assert
   implicit none

   call case_initial_unbounded_iter0()
   call case_initial_bounded_iter0()
   call case_bounded_iter_pos_lower()
   call case_bounded_iter_pos_upper()
   call case_ascent_direction()
   call case_continuation_path()

   write(*, '("test_lnsrlb: PASS")')

contains

   subroutine case_initial_unbounded_iter0()
      ! cnstnd=F (no bounds), iter=0, .not. boxed -> stp = min(one/dnorm, big).
      integer, parameter :: n = 1
      integer  :: nbd(n), iter, ifun, iback, nfgv, info, isave(2)
      real(dp) :: l(n), u(n), x(n), g(n), d(n), r(n), t(n), z(n)
      real(dp) :: f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
      character(len=60) :: task, csave
      logical  :: boxed, cnstnd
      real(dp) :: dsave(13)

      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 0.0_dp /)
      g = (/ -1.0_dp /)
      d = (/ 1.0_dp /)
      r = 0.0_dp; t = 0.0_dp
      z = (/ 1.0_dp /)
      f = 0.0_dp
      stp = -42.0_dp; dnorm = 0.0_dp; dtd = 0.0_dp; xstep = 0.0_dp; stpmx = 0.0_dp
      iter = 0; ifun = -1; iback = -1; nfgv = 0; info = 0
      task = 'START'; csave = ''
      boxed = .false.; cnstnd = .false.
      isave = 0; dsave = 0.0_dp

      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
                  stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, &
                  info, task, boxed, cnstnd, csave, isave, dsave)
      call assert_close_real(dnorm, 1.0_dp, where="case_initial_unbounded_iter0 dnorm")
      call assert_close_real(stp,   1.0_dp, where="case_initial_unbounded_iter0 stp")
      call assert_eq_int(info, 0, where="case_initial_unbounded_iter0 info")
      ! After dcsrch returns 'FG', lnsrlb sets task='FG_LNSRCH' and increments ifun.
      call assert_true(task(1:9) == 'FG_LNSRCH', "case_initial_unbounded_iter0 task")
      call assert_eq_int(ifun, 1, where="case_initial_unbounded_iter0 ifun")
   end subroutine case_initial_unbounded_iter0

   subroutine case_initial_bounded_iter0()
      ! cnstnd=T, iter=0: stpmx = 1.
      integer, parameter :: n = 1
      integer  :: nbd(n), iter, ifun, iback, nfgv, info, isave(2)
      real(dp) :: l(n), u(n), x(n), g(n), d(n), r(n), t(n), z(n)
      real(dp) :: f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
      character(len=60) :: task, csave
      logical  :: boxed, cnstnd
      real(dp) :: dsave(13)

      nbd = (/ 2 /)
      l = (/ -1.0_dp /); u = (/ 1.0_dp /)
      x = (/ 0.0_dp /); g = (/ -1.0_dp /); d = (/ 0.5_dp /)
      r = 0.0_dp; t = 0.0_dp; z = (/ 0.5_dp /)
      f = 0.0_dp
      iter = 0; ifun = 0; iback = 0; nfgv = 0; info = 0
      task = 'START'; csave = ''
      boxed = .true.; cnstnd = .true.
      isave = 0; dsave = 0.0_dp

      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
                  stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, &
                  info, task, boxed, cnstnd, csave, isave, dsave)
      call assert_close_real(stpmx, 1.0_dp, where="case_initial_bounded_iter0 stpmx")
      call assert_eq_int(info, 0, where="case_initial_bounded_iter0 info")
   end subroutine case_initial_bounded_iter0

   subroutine case_bounded_iter_pos_lower()
      ! cnstnd=T, iter>0, d<0 with lower bound: stpmx limited by (l-x)/d.
      integer, parameter :: n = 1
      integer  :: nbd(n), iter, ifun, iback, nfgv, info, isave(2)
      real(dp) :: l(n), u(n), x(n), g(n), d(n), r(n), t(n), z(n)
      real(dp) :: f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
      character(len=60) :: task, csave
      logical  :: boxed, cnstnd
      real(dp) :: dsave(13)

      nbd = (/ 2 /)
      l = (/ -1.0_dp /); u = (/ 5.0_dp /)
      x = (/ 0.0_dp /); g = (/ 1.0_dp /); d = (/ -2.0_dp /)
      r = 0.0_dp; t = 0.0_dp; z = (/ -1.0_dp /)
      f = 0.0_dp
      iter = 1; ifun = 0; iback = 0; nfgv = 0; info = 0
      task = 'START'; csave = ''
      boxed = .true.; cnstnd = .true.
      isave = 0; dsave = 0.0_dp

      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
                  stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, &
                  info, task, boxed, cnstnd, csave, isave, dsave)
      ! a1=d=-2 < 0 with nbd=2; a2 = l - x = -1 < 0; a1*stpmx = -2*big < a2 -> stpmx = a2/a1 = 0.5.
      call assert_close_real(stpmx, 0.5_dp, where="case_bounded_iter_pos_lower stpmx")
   end subroutine case_bounded_iter_pos_lower

   subroutine case_bounded_iter_pos_upper()
      ! cnstnd=T, iter>0, d>0 with upper bound: stpmx limited by (u-x)/d.
      integer, parameter :: n = 1
      integer  :: nbd(n), iter, ifun, iback, nfgv, info, isave(2)
      real(dp) :: l(n), u(n), x(n), g(n), d(n), r(n), t(n), z(n)
      real(dp) :: f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
      character(len=60) :: task, csave
      logical  :: boxed, cnstnd
      real(dp) :: dsave(13)

      nbd = (/ 2 /)
      l = (/ -5.0_dp /); u = (/ 1.0_dp /)
      x = (/ 0.0_dp /); g = (/ -1.0_dp /); d = (/ 2.0_dp /)
      r = 0.0_dp; t = 0.0_dp; z = (/ 1.0_dp /)
      f = 0.0_dp
      iter = 1; ifun = 0; iback = 0; nfgv = 0; info = 0
      task = 'START'; csave = ''
      boxed = .true.; cnstnd = .true.
      isave = 0; dsave = 0.0_dp

      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
                  stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, &
                  info, task, boxed, cnstnd, csave, isave, dsave)
      ! a1=d=2 > 0 with nbd=2; a2 = u - x = 1 > 0; a1*stpmx > a2 -> stpmx = a2/a1 = 0.5.
      call assert_close_real(stpmx, 0.5_dp, where="case_bounded_iter_pos_upper stpmx")
   end subroutine case_bounded_iter_pos_upper

   subroutine case_ascent_direction()
      ! gd = g'd >= 0: not a descent direction. info = -4 set.
      ! lnsrlb prints a diagnostic line ("ascent direction in projection") to
      ! unit 6; that's expected and harmless for the test.
      integer, parameter :: n = 1
      integer  :: nbd(n), iter, ifun, iback, nfgv, info, isave(2)
      real(dp) :: l(n), u(n), x(n), g(n), d(n), r(n), t(n), z(n)
      real(dp) :: f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
      character(len=60) :: task, csave
      logical  :: boxed, cnstnd
      real(dp) :: dsave(13)

      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 0.0_dp /); g = (/ 1.0_dp /); d = (/ 1.0_dp /)   ! g'd = 1 > 0
      r = 0.0_dp; t = 0.0_dp; z = 0.0_dp
      f = 0.0_dp
      iter = 0; ifun = 0; iback = 0; nfgv = 0; info = 0
      task = 'START'; csave = ''
      boxed = .false.; cnstnd = .false.
      isave = 0; dsave = 0.0_dp

      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
                  stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, &
                  info, task, boxed, cnstnd, csave, isave, dsave)
      call assert_eq_int(info, -4, where="case_ascent_direction info")
   end subroutine case_ascent_direction

   subroutine case_continuation_path()
      ! Two-call sequence: first call hits the setup path, second call starts
      ! with task='FG_LN' and skips to the dcsrch reentry (L75 -> 556).
      integer, parameter :: n = 1
      integer  :: nbd(n), iter, ifun, iback, nfgv, info, isave(2)
      real(dp) :: l(n), u(n), x(n), g(n), d(n), r(n), t(n), z(n)
      real(dp) :: f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
      character(len=60) :: task, csave
      logical  :: boxed, cnstnd
      real(dp) :: dsave(13)

      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 0.0_dp /); g = (/ -2.0_dp /); d = (/ 1.0_dp /)
      r = 0.0_dp; t = 0.0_dp; z = (/ 1.0_dp /)
      f = 0.0_dp
      iter = 0; ifun = 0; iback = 0; nfgv = 0; info = 0
      task = 'START'; csave = ''
      boxed = .false.; cnstnd = .false.
      isave = 0; dsave = 0.0_dp

      ! First call: runs setup, dcsrch returns 'FG', lnsrlb sets task='FG_LNSRCH'.
      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
                  stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, &
                  info, task, boxed, cnstnd, csave, isave, dsave)

      ! User would now compute f and g at the new x. Simulate that.
      f = (x(1) - 1.0_dp)**2 - 1.0_dp
      g = (/ 2.0_dp * (x(1) - 1.0_dp) /)

      ! Second call: task starts with 'FG_LN' so lnsrlb skips setup.
      task = 'FG_LN'
      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, &
                  stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, &
                  info, task, boxed, cnstnd, csave, isave, dsave)
      call assert_eq_int(info, 0, where="case_continuation_path info")
   end subroutine case_continuation_path

end program test_lnsrlb
