program test_setulb
   use test_assert, only: dp, assert_true, assert_eq_int, assert_eq_str, assert_close_real, assert_array_close
   implicit none

   call case_start_initialises_workspace_offsets()
   call case_continuation_uses_existing_offsets()

   write(*, '("test_setulb: PASS")')

contains

   subroutine case_start_initialises_workspace_offsets()
      ! task='START': setulb should populate isave(1..16) with workspace
      ! sizes and offsets, then forward to mainlb. mainlb returns 'FG'
      ! (asking for f, g at the initial point).
      integer, parameter :: n = 2, m = 3
      integer  :: nbd(n), iwa(3*n), iprint, isave(44)
      real(dp) :: x(n), l(n), u(n), g(n), f, factr, pgtol
      real(dp) :: dsave(29)
      real(dp) :: wa(2*m*n + 5*n + 11*m*m + 8*m)
      logical  :: lsave(4)
      character(len=60) :: task, csave

      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 0.5_dp, -0.5_dp /)
      g = 0.0_dp
      f = 0.0_dp
      factr = 1.0e7_dp; pgtol = 1.0e-5_dp
      iprint = -1
      task = 'START'; csave = ''
      isave = 0; dsave = 0.0_dp; lsave = .false.

      call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                  task, iprint, csave, lsave, isave, dsave)
      ! Workspace offsets: see setulb.f.
      call assert_eq_int(isave(1), m*n,   where="case_start isave(1) m*n")
      call assert_eq_int(isave(2), m*m,   where="case_start isave(2) m^2")
      call assert_eq_int(isave(3), 4*m*m, where="case_start isave(3) 4m^2")
      call assert_eq_int(isave(4), 1,     where="case_start isave(4) ws offset")
      call assert_eq_int(isave(5), 1 + m*n, where="case_start isave(5) wy offset")
      ! ws (m*n) + wy (m*n) + sy (m^2) + ss (m^2) + wt (m^2) + wn (4*m^2) +
      ! snd (4*m^2) + z, r, d, t, xp (5*n) -> 2*m*n + 11*m^2 + 5*n.
      call assert_eq_int(isave(16), 1 + 2*m*n + 11*m*m + 5*n, &
                         where="case_start isave(16) wa offset")
      ! First setulb call returns 'FG' (request function/gradient).
      call assert_true(task(1:2) == 'FG', "case_start task=FG on return")
   end subroutine case_start_initialises_workspace_offsets

   subroutine case_continuation_uses_existing_offsets()
      ! After a START call, subsequent calls feed back f, g and progress
      ! the iteration. This case drives a few iterations of the n=2
      ! bound-Rosenbrock to verify the workspace plumbing is consistent.
      integer, parameter :: n = 2, m = 3
      integer  :: nbd(n), iwa(3*n), iprint, isave(44)
      real(dp) :: x(n), l(n), u(n), g(n), f, factr, pgtol
      real(dp) :: dsave(29)
      real(dp) :: wa(2*m*n + 5*n + 11*m*m + 8*m)
      logical  :: lsave(4)
      character(len=60) :: task, csave
      integer  :: iter

      nbd = (/ 2, 2 /)
      l = (/ -10.0_dp, -10.0_dp /); u = (/ 10.0_dp, 10.0_dp /)
      x = (/ -1.2_dp, 1.0_dp /)         ! Rosenbrock starting point
      g = 0.0_dp; f = 0.0_dp
      factr = 1.0e7_dp; pgtol = 1.0e-5_dp
      iprint = -1
      task = 'START'; csave = ''
      isave = 0; dsave = 0.0_dp; lsave = .false.

      do iter = 1, 200
         call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                     task, iprint, csave, lsave, isave, dsave)
         if (task(1:2) == 'FG') then
            ! Rosenbrock f and gradient
            f = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2
            g(1) = -2.0_dp * (1.0_dp - x(1)) - 400.0_dp * x(1) * (x(2) - x(1)**2)
            g(2) = 200.0_dp * (x(2) - x(1)**2)
            cycle
         end if
         if (task(1:5) == 'NEW_X') cycle
         exit
      end do
      ! Should reach convergence and return 'CONV...' on a feasible problem.
      call assert_true(task(1:4) == 'CONV', "case_continuation final task=CONV")
      ! Optimum is at x = (1, 1) with f = 0.
      call assert_close_real(f, 0.0_dp, rtol=1.0e-3_dp, where="case_continuation f")
   end subroutine case_continuation_uses_existing_offsets

end program test_setulb
