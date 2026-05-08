program test_active
   use test_assert
   implicit none

   call case_all_unbounded()
   call case_nbd2_interior()
   call case_nbd2_at_lower_no_project()
   call case_nbd2_below_lower_projects()
   call case_nbd2_at_upper_no_project()
   call case_nbd2_above_upper_projects()
   call case_nbd1_below_lower_projects()
   call case_nbd3_above_upper_projects()
   call case_nbd2_fixed_variable()
   call case_iprint_unconstrained()
   call case_iprint_with_bounded_variable()

   write(*, '("test_active: PASS")')

contains

   subroutine case_all_unbounded()
      ! nbd=0 for all variables: cnstnd=false, boxed=false, iwhere=-1.
      ! No projection. iprint<0 silences output.
      integer, parameter :: n = 2
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 0; l = 0.0_dp; u = 0.0_dp
      x = (/ 1.5_dp, -2.0_dp /)
      iwhere = 0
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_true(.not. prjctd, "case_all_unbounded prjctd")
      call assert_true(.not. cnstnd, "case_all_unbounded cnstnd")
      call assert_true(.not. boxed,  "case_all_unbounded boxed")
      call assert_eq_int(iwhere(1), -1, where="case_all_unbounded iwhere(1)")
      call assert_eq_int(iwhere(2), -1, where="case_all_unbounded iwhere(2)")
      call assert_close_real(x(1), 1.5_dp, where="case_all_unbounded x(1)")
      call assert_close_real(x(2), -2.0_dp, where="case_all_unbounded x(2)")
   end subroutine case_all_unbounded

   subroutine case_nbd2_interior()
      ! All nbd=2 with x strictly inside (l, u): boxed=true, no projection.
      integer, parameter :: n = 2
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 2
      l = (/ 0.0_dp, -1.0_dp /); u = (/ 10.0_dp, 1.0_dp /)
      x = (/ 5.0_dp, 0.0_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_true(.not. prjctd, "case_nbd2_interior prjctd")
      call assert_true(cnstnd, "case_nbd2_interior cnstnd")
      call assert_true(boxed,  "case_nbd2_interior boxed")
      call assert_eq_int(iwhere(1), 0, where="case_nbd2_interior iwhere(1)")
      call assert_eq_int(iwhere(2), 0, where="case_nbd2_interior iwhere(2)")
   end subroutine case_nbd2_interior

   subroutine case_nbd2_at_lower_no_project()
      ! nbd=2 with x exactly at l: enters L56 (x<=l) but skips L57 (x<l strict).
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 2; l = (/ 1.0_dp /); u = (/ 5.0_dp /); x = l
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_true(.not. prjctd, "case_nbd2_at_lower_no_project prjctd")
      call assert_close_real(x(1), 1.0_dp, where="case_nbd2_at_lower_no_project x")
   end subroutine case_nbd2_at_lower_no_project

   subroutine case_nbd2_below_lower_projects()
      ! x strictly below l: projects x up, sets prjctd=true.
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 2; l = (/ 1.0_dp /); u = (/ 5.0_dp /); x = (/ -3.0_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, 0, prjctd, cnstnd, boxed)
      call assert_true(prjctd, "case_nbd2_below_lower_projects prjctd")
      call assert_close_real(x(1), 1.0_dp, where="case_nbd2_below_lower_projects x")
   end subroutine case_nbd2_below_lower_projects

   subroutine case_nbd2_at_upper_no_project()
      ! nbd=2 with x exactly at u: enters L62 (x>=u) but skips L63 (x>u strict).
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 2; l = (/ 0.0_dp /); u = (/ 5.0_dp /); x = u
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_true(.not. prjctd, "case_nbd2_at_upper_no_project prjctd")
      call assert_close_real(x(1), 5.0_dp, where="case_nbd2_at_upper_no_project x")
   end subroutine case_nbd2_at_upper_no_project

   subroutine case_nbd2_above_upper_projects()
      ! x strictly above u: projects x down to u.
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 2; l = (/ 0.0_dp /); u = (/ 5.0_dp /); x = (/ 9.0_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_true(prjctd, "case_nbd2_above_upper_projects prjctd")
      call assert_close_real(x(1), 5.0_dp, where="case_nbd2_above_upper_projects x")
   end subroutine case_nbd2_above_upper_projects

   subroutine case_nbd1_below_lower_projects()
      ! nbd=1 (lower bound only) with x<l: projects up.
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 1; l = (/ 2.0_dp /); u = 0.0_dp; x = (/ -1.0_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_true(prjctd, "case_nbd1_below_lower_projects prjctd")
      call assert_true(.not. boxed, "case_nbd1_below_lower_projects boxed")
      call assert_close_real(x(1), 2.0_dp, where="case_nbd1_below_lower_projects x")
   end subroutine case_nbd1_below_lower_projects

   subroutine case_nbd3_above_upper_projects()
      ! nbd=3 (upper bound only) with x>u: projects down.
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 3; l = 0.0_dp; u = (/ 2.0_dp /); x = (/ 7.0_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_true(prjctd, "case_nbd3_above_upper_projects prjctd")
      call assert_close_real(x(1), 2.0_dp, where="case_nbd3_above_upper_projects x")
   end subroutine case_nbd3_above_upper_projects

   subroutine case_nbd2_fixed_variable()
      ! nbd=2 with l == u: variable is fixed; iwhere=3.
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 2; l = (/ 4.0_dp /); u = (/ 4.0_dp /); x = (/ 4.0_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, -1, prjctd, cnstnd, boxed)
      call assert_eq_int(iwhere(1), 3, where="case_nbd2_fixed_variable iwhere")
   end subroutine case_nbd2_fixed_variable

   subroutine case_iprint_unconstrained()
      ! iprint>=0 with cnstnd=false (nbd=0 only): triggers the "unconstrained"
      ! diagnostic line. Also covers L93 false (no projection happened).
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 0; l = 0.0_dp; u = 0.0_dp; x = (/ 1.0_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, 0, prjctd, cnstnd, boxed)
      call assert_true(.not. cnstnd, "case_iprint_unconstrained cnstnd")
   end subroutine case_iprint_unconstrained

   subroutine case_iprint_with_bounded_variable()
      ! iprint>0 with cnstnd=true: triggers the nbdd-count diagnostic.
      ! No projection so L93 is false; cnstnd is true so L95 is false.
      integer, parameter :: n = 1
      integer :: nbd(n), iwhere(n)
      real(dp) :: l(n), u(n), x(n)
      logical :: prjctd, cnstnd, boxed
      nbd = 2; l = (/ 0.0_dp /); u = (/ 1.0_dp /); x = (/ 0.5_dp /)
      iwhere = -99
      call active(n, l, u, nbd, x, iwhere, 1, prjctd, cnstnd, boxed)
      call assert_true(.not. prjctd, "case_iprint_with_bounded_variable prjctd")
      call assert_true(cnstnd, "case_iprint_with_bounded_variable cnstnd")
   end subroutine case_iprint_with_bounded_variable

end program test_active
