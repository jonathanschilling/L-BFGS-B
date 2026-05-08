program test_mainlb
   use test_assert
   implicit none

   ! mainlb is the state-machine driver. Exercising it directly requires
   ! correctly partitioned workspaces; we go through setulb instead.

   call case_invalid_input_returns_error()
   call case_pgtol_convergence()
   call case_user_signals_stop()
   call case_factr_convergence()

   write(*, '("test_mainlb: PASS")')

contains

   subroutine case_invalid_input_returns_error()
      ! n <= 0 is rejected by errclb inside mainlb. task should come back
      ! starting with 'ERROR'. n=0 still requires the wa/iwa arrays to be
      ! addressable, so set m=1 and use minimum-size buffers.
      integer, parameter :: m = 1
      integer  :: n_bad
      integer  :: nbd_dummy(1), iwa(3), iprint, isave(44)
      real(dp) :: x(1), l(1), u(1), g(1), f, factr, pgtol
      real(dp) :: dsave(29), wa(2*m + 5 + 11*m*m + 8*m)
      logical  :: lsave(4)
      character(len=60) :: task, csave
      n_bad = 0
      nbd_dummy = 0; l = 0.0_dp; u = 0.0_dp
      x = 0.0_dp; g = 0.0_dp; f = 0.0_dp
      factr = 1.0e7_dp; pgtol = 1.0e-5_dp
      iprint = -1
      task = 'START'; csave = ''
      isave = 0; dsave = 0.0_dp; lsave = .false.
      call setulb(n_bad, m, x, l, u, nbd_dummy, f, g, factr, pgtol, wa, iwa, &
                  task, iprint, csave, lsave, isave, dsave)
      call assert_true(task(1:5) == 'ERROR', "case_invalid_input task=ERROR")
   end subroutine case_invalid_input_returns_error

   subroutine case_pgtol_convergence()
      ! Small gradient stopping criterion: pgtol = 1e-3 so the algorithm
      ! exits early via |proj g| <= pgtol rather than via the f-reduction
      ! criterion. This exercises mainlb's task='CONVERGENCE: NORM_...' path.
      integer, parameter :: n = 2, m = 3
      integer  :: nbd(n), iwa(3*n), iprint, isave(44)
      real(dp) :: x(n), l(n), u(n), g(n), f, factr, pgtol
      real(dp) :: dsave(29), wa(2*m*n + 5*n + 11*m*m + 8*m)
      logical  :: lsave(4)
      character(len=60) :: task, csave
      integer  :: iter
      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 5.0_dp, -3.0_dp /)
      g = 0.0_dp; f = 0.0_dp
      factr = 0.0_dp                         ! suppress f-reduction stop
      pgtol = 1.0e-3_dp
      iprint = -1
      task = 'START'; csave = ''
      isave = 0; dsave = 0.0_dp; lsave = .false.
      do iter = 1, 200
         call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                     task, iprint, csave, lsave, isave, dsave)
         if (task(1:2) == 'FG') then
            ! phi(x) = sum_i (x_i - i)^2 ; minimum at x = (1, 2).
            f = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
            g(1) = 2.0_dp * (x(1) - 1.0_dp)
            g(2) = 2.0_dp * (x(2) - 2.0_dp)
            cycle
         end if
         if (task(1:5) == 'NEW_X') cycle
         exit
      end do
      call assert_true(task(1:25) == 'CONVERGENCE: NORM_OF_PROJ', &
                       "case_pgtol_convergence task")
      call assert_close_real(x(1), 1.0_dp, rtol=1.0e-3_dp, where="case_pgtol x(1)")
      call assert_close_real(x(2), 2.0_dp, rtol=1.0e-3_dp, where="case_pgtol x(2)")
   end subroutine case_pgtol_convergence

   subroutine case_user_signals_stop()
      ! User assigns task='STOP' on a NEW_X return: mainlb should honor it
      ! and exit cleanly. Drives the goto-999 path through the STOP branch.
      integer, parameter :: n = 2, m = 3
      integer  :: nbd(n), iwa(3*n), iprint, isave(44)
      real(dp) :: x(n), l(n), u(n), g(n), f, factr, pgtol
      real(dp) :: dsave(29), wa(2*m*n + 5*n + 11*m*m + 8*m)
      logical  :: lsave(4)
      character(len=60) :: task, csave
      integer  :: iter
      logical  :: stopped
      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 0.0_dp, 0.0_dp /)
      g = 0.0_dp; f = 0.0_dp
      factr = 1.0e7_dp; pgtol = 1.0e-10_dp     ! prevent fast convergence
      iprint = -1
      task = 'START'; csave = ''
      isave = 0; dsave = 0.0_dp; lsave = .false.
      stopped = .false.
      do iter = 1, 100
         call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                     task, iprint, csave, lsave, isave, dsave)
         if (task(1:2) == 'FG') then
            f = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
            g(1) = 2.0_dp * (x(1) - 1.0_dp)
            g(2) = 2.0_dp * (x(2) - 2.0_dp)
            cycle
         end if
         if (task(1:5) == 'NEW_X') then
            ! Tell mainlb to stop on the next entry.
            task = 'STOP: USER REQUESTED'
            stopped = .true.
            cycle
         end if
         exit
      end do
      call assert_true(stopped, "case_user_signals_stop user did issue STOP")
      call assert_true(task(1:4) == 'STOP', "case_user_signals_stop task=STOP")
   end subroutine case_user_signals_stop

   subroutine case_factr_convergence()
      ! Slow gradient convergence on the same quadratic, with factr large
      ! enough that the rel-reduction-of-f criterion fires first.
      integer, parameter :: n = 2, m = 3
      integer  :: nbd(n), iwa(3*n), iprint, isave(44)
      real(dp) :: x(n), l(n), u(n), g(n), f, factr, pgtol
      real(dp) :: dsave(29), wa(2*m*n + 5*n + 11*m*m + 8*m)
      logical  :: lsave(4)
      character(len=60) :: task, csave
      integer  :: iter
      nbd = 0
      l = 0.0_dp; u = 0.0_dp
      x = (/ 5.0_dp, -3.0_dp /)
      g = 0.0_dp; f = 0.0_dp
      factr = 1.0e12_dp                          ! loose f-stop
      pgtol = 0.0_dp                             ! suppress g-stop
      iprint = -1
      task = 'START'; csave = ''
      isave = 0; dsave = 0.0_dp; lsave = .false.
      do iter = 1, 200
         call setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, &
                     task, iprint, csave, lsave, isave, dsave)
         if (task(1:2) == 'FG') then
            f = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
            g(1) = 2.0_dp * (x(1) - 1.0_dp)
            g(2) = 2.0_dp * (x(2) - 2.0_dp)
            cycle
         end if
         if (task(1:5) == 'NEW_X') cycle
         exit
      end do
      ! For a strict quadratic with factr=1e12, either the f-reduction
      ! criterion fires first or, if gradient hits machine zero, the
      ! pgtol=0 path also catches it. Either way task starts with 'CONV'.
      call assert_true(task(1:4) == 'CONV', "case_factr_convergence task")
   end subroutine case_factr_convergence

end program test_mainlb
