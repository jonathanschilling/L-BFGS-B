program test_errclb
   use test_assert, only: dp, assert_true, assert_eq_int, assert_eq_str, assert_close_real, assert_array_close
   implicit none

   call case_valid_input()
   call case_n_le_zero()
   call case_m_le_zero()
   call case_factr_negative()
   call case_invalid_nbd_below()
   call case_invalid_nbd_above()
   call case_no_feasible_solution()
   call case_nbd2_feasible_no_error()

   write(*, '("test_errclb: PASS")')

contains

   subroutine call_errclb(n, m, factr, l, u, nbd, task, info, k)
      integer, intent(in)            :: n, m
      real(dp), intent(in)           :: factr, l(n), u(n)
      integer, intent(in)            :: nbd(n)
      character(len=60), intent(inout) :: task
      integer, intent(inout)         :: info, k
      call errclb(n, m, factr, l, u, nbd, task, info, k)
   end subroutine call_errclb

   subroutine case_valid_input()
      ! All inputs valid, n=2, mixed nbd values that all pass validation.
      ! errclb does not write task or info on success, so they retain
      ! their entry values.
      integer, parameter :: n = 2
      integer  :: nbd(n), info, k
      real(dp) :: l(n), u(n), factr
      character(len=60) :: task
      nbd = (/ 0, 1 /)
      l = (/ -1.0_dp, 0.0_dp /)
      u = (/ 1.0_dp, 0.0_dp /)
      factr = 1.0e7_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(n, 5, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'START', where="case_valid_input task")
      call assert_eq_int(info, 0, where="case_valid_input info")
   end subroutine case_valid_input

   subroutine case_n_le_zero()
      ! n=0 triggers the n<=0 check; nbd loop is empty.
      integer  :: nbd(1), info, k
      real(dp) :: l(1), u(1), factr
      character(len=60) :: task
      nbd = 0; l = 0.0_dp; u = 0.0_dp; factr = 1.0e7_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(0, 5, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'ERROR: N .LE. 0', where="case_n_le_zero")
   end subroutine case_n_le_zero

   subroutine case_m_le_zero()
      ! m=0 triggers the m<=0 check.
      integer, parameter :: n = 1
      integer  :: nbd(n), info, k
      real(dp) :: l(n), u(n), factr
      character(len=60) :: task
      nbd = 0; l = 0.0_dp; u = 0.0_dp; factr = 1.0e7_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(n, 0, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'ERROR: M .LE. 0', where="case_m_le_zero")
   end subroutine case_m_le_zero

   subroutine case_factr_negative()
      ! factr<0 triggers the factr check (last sequential overwrite wins).
      integer, parameter :: n = 1
      integer  :: nbd(n), info, k
      real(dp) :: l(n), u(n), factr
      character(len=60) :: task
      nbd = 0; l = 0.0_dp; u = 0.0_dp; factr = -1.0_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(n, 5, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'ERROR: FACTR .LT. 0', where="case_factr_negative")
   end subroutine case_factr_negative

   subroutine case_invalid_nbd_below()
      ! nbd(i)<0 trips the nbd-validity check; info=-6, k records the index.
      integer, parameter :: n = 2
      integer  :: nbd(n), info, k
      real(dp) :: l(n), u(n), factr
      character(len=60) :: task
      nbd = (/ 0, -1 /)
      l = 0.0_dp; u = 0.0_dp; factr = 1.0e7_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(n, 5, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'ERROR: INVALID NBD', where="case_invalid_nbd_below task")
      call assert_eq_int(info, -6, where="case_invalid_nbd_below info")
      call assert_eq_int(k, 2, where="case_invalid_nbd_below k")
   end subroutine case_invalid_nbd_below

   subroutine case_invalid_nbd_above()
      ! nbd(i)>3 trips the nbd-validity check (other side of the range).
      integer, parameter :: n = 1
      integer  :: nbd(n), info, k
      real(dp) :: l(n), u(n), factr
      character(len=60) :: task
      nbd = (/ 4 /)
      l = 0.0_dp; u = 0.0_dp; factr = 1.0e7_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(n, 5, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'ERROR: INVALID NBD', where="case_invalid_nbd_above task")
      call assert_eq_int(info, -6, where="case_invalid_nbd_above info")
      call assert_eq_int(k, 1, where="case_invalid_nbd_above k")
   end subroutine case_invalid_nbd_above

   subroutine case_no_feasible_solution()
      ! nbd=2 with l>u: infeasible, info=-7.
      integer, parameter :: n = 1
      integer  :: nbd(n), info, k
      real(dp) :: l(n), u(n), factr
      character(len=60) :: task
      nbd = (/ 2 /)
      l = (/ 5.0_dp /); u = (/ 3.0_dp /); factr = 1.0e7_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(n, 5, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'ERROR: NO FEASIBLE SOLUTION', where="case_no_feasible_solution task")
      call assert_eq_int(info, -7, where="case_no_feasible_solution info")
      call assert_eq_int(k, 1, where="case_no_feasible_solution k")
   end subroutine case_no_feasible_solution

   subroutine case_nbd2_feasible_no_error()
      ! nbd=2 with l<=u: covers L53-true / L54-false branch combination.
      integer, parameter :: n = 1
      integer  :: nbd(n), info, k
      real(dp) :: l(n), u(n), factr
      character(len=60) :: task
      nbd = (/ 2 /)
      l = (/ 1.0_dp /); u = (/ 2.0_dp /); factr = 1.0e7_dp
      task = 'START' ; info = 0 ; k = 0
      call call_errclb(n, 5, factr, l, u, nbd, task, info, k)
      call assert_eq_str(task, 'START', where="case_nbd2_feasible_no_error task")
      call assert_eq_int(info, 0, where="case_nbd2_feasible_no_error info")
   end subroutine case_nbd2_feasible_no_error

end program test_errclb
