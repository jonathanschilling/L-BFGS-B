program test_projgr
   use test_assert
   implicit none

   call case_unbounded()
   call case_unbounded_negative_g()
   call case_nbd2_interior_positive_g()
   call case_nbd2_lower_clipped()
   call case_nbd2_upper_clipped()
   call case_nbd1_negative_g_no_upper()
   call case_nbd3_positive_g_no_lower()
   call case_max_aggregation()

   write(*, '("test_projgr: PASS")')

contains

   subroutine case_unbounded()
      ! nbd=0 -> outer if FALSE: variable contributes |g| with no projection.
      integer, parameter :: n = 1
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd(1) = 0
      l = 0.0_dp; u = 0.0_dp; x = 0.0_dp
      g = (/ 3.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 3.0_dp, where="case_unbounded")
   end subroutine case_unbounded

   subroutine case_unbounded_negative_g()
      ! nbd=0 with g<0: same path, |g| accumulates absolute value.
      integer, parameter :: n = 1
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd(1) = 0
      l = 0.0_dp; u = 0.0_dp; x = 0.0_dp
      g = (/ -3.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 3.0_dp, where="case_unbounded_negative_g")
   end subroutine case_unbounded_negative_g

   subroutine case_nbd2_interior_positive_g()
      ! nbd=2, x in interior, g>=0: hits "nbd<=2" branch.
      ! min(x-l, g) = min(4, 2) = 2 -> sbgnrm = 2.
      integer, parameter :: n = 1
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd(1) = 2
      l = (/ 1.0_dp /); u = (/ 10.0_dp /); x = (/ 5.0_dp /)
      g = (/ 2.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 2.0_dp, where="case_nbd2_interior_positive_g")
   end subroutine case_nbd2_interior_positive_g

   subroutine case_nbd2_lower_clipped()
      ! nbd=2, x near lower bound, g>=0: gradient pushes toward lower.
      ! gi = min(x-l, g) = min(0.5, 2) = 0.5.
      integer, parameter :: n = 1
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd(1) = 2
      l = (/ 1.0_dp /); u = (/ 10.0_dp /); x = (/ 1.5_dp /)
      g = (/ 2.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 0.5_dp, where="case_nbd2_lower_clipped")
   end subroutine case_nbd2_lower_clipped

   subroutine case_nbd2_upper_clipped()
      ! nbd=2, x near upper bound, g<0: gradient pushes toward upper.
      ! gi = max(x-u, g) = max(-1, -2) = -1, |gi| = 1.
      integer, parameter :: n = 1
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd(1) = 2
      l = (/ 1.0_dp /); u = (/ 10.0_dp /); x = (/ 9.0_dp /)
      g = (/ -2.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 1.0_dp, where="case_nbd2_upper_clipped")
   end subroutine case_nbd2_upper_clipped

   subroutine case_nbd1_negative_g_no_upper()
      ! nbd=1 (lower only), g<0: variable can move up freely.
      ! Inner branch "nbd>=2" is FALSE -> no clip, gi=g, sbgnrm=|g|.
      integer, parameter :: n = 1
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd(1) = 1
      l = (/ 1.0_dp /); u = 0.0_dp; x = (/ 5.0_dp /)
      g = (/ -7.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 7.0_dp, where="case_nbd1_negative_g_no_upper")
   end subroutine case_nbd1_negative_g_no_upper

   subroutine case_nbd3_positive_g_no_lower()
      ! nbd=3 (upper only), g>=0: variable can move down freely.
      ! Inner branch "nbd<=2" is FALSE -> no clip, gi=g, sbgnrm=g.
      integer, parameter :: n = 1
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd(1) = 3
      l = 0.0_dp; u = (/ 10.0_dp /); x = (/ 5.0_dp /)
      g = (/ 4.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 4.0_dp, where="case_nbd3_positive_g_no_lower")
   end subroutine case_nbd3_positive_g_no_lower

   subroutine case_max_aggregation()
      ! Multi-component: verify max() across components, mixed bound types.
      ! i=1 nbd=0, g=1.5  -> contributes 1.5
      ! i=2 nbd=2 interior, g=-2.0 in (l=0,u=10), x=5: max(-5,-2)=-2 -> 2.0
      ! i=3 nbd=2 upper near, g=-3, x=8, u=10: max(-2,-3)=-2 -> 2.0
      ! i=4 nbd=1 lower near, g=4, x=1.1, l=1: min(0.1,4)=0.1 -> 0.1
      ! Expected sbgnrm = 2.0 (from i=2 or i=3).
      integer, parameter :: n = 4
      integer :: nbd(n)
      real(dp) :: l(n), u(n), x(n), g(n), sbgnrm
      nbd = (/ 0, 2, 2, 1 /)
      l = (/ 0.0_dp, 0.0_dp,  0.0_dp, 1.0_dp /)
      u = (/ 0.0_dp, 10.0_dp, 10.0_dp, 0.0_dp /)
      x = (/ 0.0_dp, 5.0_dp,  8.0_dp, 1.1_dp /)
      g = (/ 1.5_dp, -2.0_dp, -3.0_dp, 4.0_dp /)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call assert_close_real(sbgnrm, 2.0_dp, where="case_max_aggregation")
   end subroutine case_max_aggregation

end program test_projgr
