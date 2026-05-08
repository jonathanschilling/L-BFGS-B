program test_dcsrch
   use test_assert, only: dp, assert_true, assert_eq_int, assert_eq_str, assert_close_real, assert_array_close
   implicit none

   ! dcsrch is the Moré–Thuente line search driver, called via reverse
   ! communication: the user evaluates phi(stp) and phi'(stp), then calls
   ! again. Test cases below cover input-validation, normal convergence on
   ! a quadratic, and the STPMAX hit warning.

   call case_err_stp_below_stpmin()
   call case_err_stp_above_stpmax()
   call case_err_initial_g_nonneg()
   call case_err_ftol_negative()
   call case_err_gtol_negative()
   call case_err_xtol_negative()
   call case_err_stpmin_negative()
   call case_err_stpmax_below_stpmin()
   call case_quadratic_converges()
   call case_warning_stp_at_stpmax()
   call case_phi1_more_thuente()
   call case_modified_function_path()

   write(*, '("test_dcsrch: PASS")')

contains

   subroutine call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      real(dp), intent(inout) :: f, g, stp
      real(dp), intent(in)    :: ftol, gtol, xtol, stpmin, stpmax
      character(len=60), intent(inout) :: task
      integer :: isave(2)
      real(dp) :: dsave(13)
      isave = 0; dsave = 0.0_dp
      call dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task, isave, dsave)
   end subroutine call_dcsrch

   subroutine case_err_stp_below_stpmin()
      ! stp < stpmin, all other validation checks pass.
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = -1.0_dp; stp = 0.05_dp
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = 0.1_dp; stpmax = 1.0_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: STP .LT. STPMIN', where="case_err_stp_below_stpmin")
   end subroutine case_err_stp_below_stpmin

   subroutine case_err_stp_above_stpmax()
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = -1.0_dp; stp = 5.0_dp
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = 0.1_dp; stpmax = 1.0_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: STP .GT. STPMAX', where="case_err_stp_above_stpmax")
   end subroutine case_err_stp_above_stpmax

   subroutine case_err_initial_g_nonneg()
      ! g >= 0 at t=0 means the search direction isn't a descent direction.
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = 0.5_dp; stp = 0.5_dp
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 1.0_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: INITIAL G .GE. ZERO', where="case_err_initial_g_nonneg")
   end subroutine case_err_initial_g_nonneg

   subroutine case_err_ftol_negative()
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = -1.0_dp; stp = 0.5_dp
      ftol = -1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 1.0_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: FTOL .LT. ZERO', where="case_err_ftol_negative")
   end subroutine case_err_ftol_negative

   subroutine case_err_gtol_negative()
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = -1.0_dp; stp = 0.5_dp
      ftol = 1.0e-3_dp; gtol = -0.1_dp; xtol = 1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 1.0_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: GTOL .LT. ZERO', where="case_err_gtol_negative")
   end subroutine case_err_gtol_negative

   subroutine case_err_xtol_negative()
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = -1.0_dp; stp = 0.5_dp
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = -1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 1.0_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: XTOL .LT. ZERO', where="case_err_xtol_negative")
   end subroutine case_err_xtol_negative

   subroutine case_err_stpmin_negative()
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = -1.0_dp; stp = 0.5_dp
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = -0.1_dp; stpmax = 1.0_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: STPMIN .LT. ZERO', where="case_err_stpmin_negative")
   end subroutine case_err_stpmin_negative

   subroutine case_err_stpmax_below_stpmin()
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      f = 0.0_dp; g = -1.0_dp; stp = 0.5_dp
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = 1.0_dp; stpmax = 0.1_dp
      task = 'START'
      call call_dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task)
      call assert_eq_str(task, 'ERROR: STPMAX .LT. STPMIN', where="case_err_stpmax_below_stpmin")
   end subroutine case_err_stpmax_below_stpmin

   subroutine case_quadratic_converges()
      ! phi(t) = (t-1)^2 - 1: minimum at t=1 with phi(1)=-1, phi'(1)=0.
      ! At t=0: phi=0, phi'=-2 (descent direction).
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      integer :: isave(2), iter
      real(dp) :: dsave(13)

      f = 0.0_dp; g = -2.0_dp; stp = 0.5_dp
      ftol = 1.0e-3_dp; gtol = 0.1_dp; xtol = 1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 5.0_dp
      task = 'START'
      isave = 0; dsave = 0.0_dp
      do iter = 1, 50
         call dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task, isave, dsave)
         if (task(1:2) /= 'FG') exit
         f = (stp - 1.0_dp)**2 - 1.0_dp
         g = 2.0_dp * (stp - 1.0_dp)
      end do
      call assert_true(task(1:4) == 'CONV', "case_quadratic_converges task=CONV")
      ! At convergence, |phi'(stp)| should be small. For a quadratic with
      ! phi'(1) = 0, the minimiser is at t=1; allow small numerical slack.
      call assert_close_real(stp, 1.0_dp, rtol=1.0e-3_dp, where="case_quadratic_converges stp")
   end subroutine case_quadratic_converges

   subroutine case_warning_stp_at_stpmax()
      ! Use phi(t) = -t (linear, monotone decreasing): no minimum exists.
      ! Search will push stp up to stpmax. Algorithm exits with
      ! 'WARNING: STP = STPMAX'.
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      integer :: isave(2), iter
      real(dp) :: dsave(13)
      f = 0.0_dp; g = -1.0_dp; stp = 0.5_dp
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 2.0_dp
      task = 'START'
      isave = 0; dsave = 0.0_dp
      do iter = 1, 50
         call dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task, isave, dsave)
         if (task(1:2) /= 'FG') exit
         f = -stp
         g = -1.0_dp
      end do
      ! For a strictly-descending phi, the search hits stpmax and warns.
      call assert_true(task(1:4) == 'WARN', "case_warning_stp_at_stpmax task=WARN")
      call assert_close_real(stp, stpmax, where="case_warning_stp_at_stpmax stp")
   end subroutine case_warning_stp_at_stpmax

   subroutine case_phi1_more_thuente()
      ! Moré–Thuente test function phi1(t) = -t/(t^2 + b), b=2.
      ! Minimum at t = sqrt(b) = sqrt(2) ≈ 1.4142.
      ! With a small initial step (0.001) the line search has to travel
      ! far and ends up bracketing once it overshoots, which exercises
      ! the brackt=true paths (state restore L172, bisection L259-260,
      ! brackt-stmin/stmax block L267).
      real(dp), parameter :: b = 2.0_dp
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      integer :: isave(2), iter
      real(dp) :: dsave(13)
      f = 0.0_dp
      g = -1.0_dp / b                            ! phi1'(0)
      stp = 1.0e-3_dp
      ftol = 1.0e-3_dp; gtol = 0.1_dp; xtol = 1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 1.0e3_dp
      task = 'START'
      isave = 0; dsave = 0.0_dp
      do iter = 1, 100
         call dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task, isave, dsave)
         if (task(1:2) /= 'FG') exit
         f = -stp / (stp*stp + b)
         g = (stp*stp - b) / ((stp*stp + b)**2)
      end do
      call assert_true(task(1:4) == 'CONV', "case_phi1_more_thuente task")
      ! dcsrch only guarantees the strong-Wolfe curvature condition
      !   |phi'(stp)| <= gtol * |phi'(0)|
      ! not stp exactly at the minimum. Verify that curvature condition.
      call assert_true(abs(g) <= gtol * (1.0_dp / b) + 1.0e-14_dp, &
                       "case_phi1 curvature condition satisfied")
      ! Also assert stp is in the right neighbourhood of sqrt(b).
      call assert_true(stp > 1.0_dp .and. stp < 2.0_dp, &
                       "case_phi1 stp in (1, 2)")
   end subroutine case_phi1_more_thuente

   subroutine case_modified_function_path()
      ! Drives the line search through the "psi" / modified-function branch
      ! at L225-246. To reach it we need stage=1, f<=fx, f>ftest -- meaning
      ! the step decreased f but by less than the sufficient-decrease line
      ! ftest = finit + stp*ftol*ginit predicts. Starting from a step that
      ! lands AT the minimum of phi(t) = -t + 0.5*t^2 (where phi' = 0)
      ! gives exactly that condition.
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      integer :: isave(2), iter
      real(dp) :: dsave(13)
      f = 0.0_dp
      g = -1.0_dp
      stp = 2.0_dp                   ! overshoots the t=1 minimum
      ftol = 1.0e-3_dp; gtol = 0.9_dp; xtol = 1.0e-10_dp
      stpmin = 0.0_dp; stpmax = 10.0_dp
      task = 'START'
      isave = 0; dsave = 0.0_dp
      do iter = 1, 50
         call dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task, isave, dsave)
         if (task(1:2) /= 'FG') exit
         f = -stp + 0.5_dp*stp*stp     ! phi(t)  = -t + t^2 / 2
         g = -1.0_dp + stp             ! phi'(t) = -1 + t
      end do
      call assert_true(task(1:4) == 'CONV', "case_modified_function_path task")
   end subroutine case_modified_function_path

end program test_dcsrch
