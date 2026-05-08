program test_hpsolb
   use test_assert
   implicit none

   call case_n1_iheap0()
   call case_n1_iheap1()
   call case_unsorted_iheap0()
   call case_prebuilt_iheap1_right_child_smaller()
   call case_iheap0_min_at_first()
   call case_sift_down_stops_early()

   write(*, '("test_hpsolb: PASS")')

contains

   subroutine assert_min_at_n(t_in, t_out, iorder_out, n, where)
      ! Verify (a) t_out(n) = min(t_in), (b) iorder_out(n) points at the
      ! original location of that minimum, (c) the multiset
      ! { (t_in(iorder_out(i)), iorder_out(i)) }_i equals the input multiset.
      integer, intent(in) :: n
      real(dp), intent(in) :: t_in(n), t_out(n)
      integer, intent(in) :: iorder_out(n)
      character(len=*), intent(in) :: where
      real(dp) :: tmin, expected
      integer :: i, idx
      tmin = minval(t_in)
      call assert_close_real(t_out(n), tmin, where=where // " t_out(n)")
      ! Check iorder_out(n) is some index of the min in t_in.
      idx = iorder_out(n)
      call assert_close_real(t_in(idx), tmin, where=where // " iorder_out(n) idx")
      ! Check each surviving slot's t-value matches t_in at that index.
      do i = 1, n
         expected = t_in(iorder_out(i))
         call assert_close_real(t_out(i), expected, where=where // " t_out(i) consistency")
      end do
   end subroutine assert_min_at_n

   subroutine case_n1_iheap0()
      ! Trivial: single-element with iheap=0. Heap-build loop is empty;
      ! extract phase guard (n>1) is false so nothing happens.
      integer, parameter :: n = 1
      real(dp) :: t(n)
      integer :: iorder(n)
      t = (/ 7.5_dp /); iorder = (/ 1 /)
      call hpsolb(n, t, iorder, 0)
      call assert_close_real(t(1), 7.5_dp, where="case_n1_iheap0 t")
      call assert_eq_int(iorder(1), 1, where="case_n1_iheap0 iorder")
   end subroutine case_n1_iheap0

   subroutine case_n1_iheap1()
      ! n=1 with iheap=1: skips heap-build (L42 false), n>1 guard also false.
      integer, parameter :: n = 1
      real(dp) :: t(n)
      integer :: iorder(n)
      t = (/ -3.25_dp /); iorder = (/ 5 /)
      call hpsolb(n, t, iorder, 1)
      call assert_close_real(t(1), -3.25_dp, where="case_n1_iheap1 t")
      call assert_eq_int(iorder(1), 5, where="case_n1_iheap1 iorder")
   end subroutine case_n1_iheap1

   subroutine case_unsorted_iheap0()
      ! Build heap from arbitrary input, then extract min.
      ! Input  t = [3, 1, 4, 2] -> exercises sift-up swap (L55 T) and
      ! no-swap (L55 F at k=3 inserting 4). Subsequent extract exercises
      ! L82 F (right child not smaller) and L81 F (sift-down reaches leaf).
      integer, parameter :: n = 4
      real(dp) :: t(n)
      integer :: iorder(n)
      real(dp) :: t_in(n)
      integer :: i
      t_in = (/ 3.0_dp, 1.0_dp, 4.0_dp, 2.0_dp /)
      t = t_in
      do i = 1, n
         iorder(i) = i
      end do
      call hpsolb(n, t, iorder, 0)
      call assert_min_at_n(t_in, t, iorder, n, "case_unsorted_iheap0")
      call assert_eq_int(iorder(n), 2, where="case_unsorted_iheap0 iorder(n)")
   end subroutine case_unsorted_iheap0

   subroutine case_prebuilt_iheap1_right_child_smaller()
      ! Pre-built heap, exercises L42 F (skip build) and L82 T (right child
      ! smaller than left child during sift-down). Heap form: [1, 3, 2, 5].
      integer, parameter :: n = 4
      real(dp) :: t(n), t_in(n)
      integer :: iorder(n), iorder_in(n)
      integer :: i
      t_in = (/ 1.0_dp, 3.0_dp, 2.0_dp, 5.0_dp /)
      iorder_in = (/ 1, 2, 3, 4 /)
      t = t_in; iorder = iorder_in
      call hpsolb(n, t, iorder, 1)
      call assert_min_at_n(t_in, t, iorder, n, "case_prebuilt_iheap1")
      call assert_eq_int(iorder(n), 1, where="case_prebuilt_iheap1 iorder(n)")
   end subroutine case_prebuilt_iheap1_right_child_smaller

   subroutine case_sift_down_stops_early()
      ! Heap [1, 3, 2] with iheap=1: ddum = t(3) = 2 equals the chosen
      ! sift child after L82's swap, so the sift-down inequality
      ! "t(j) < ddum" is false at the first iteration (L83 F branch).
      integer, parameter :: n = 3
      real(dp) :: t(n), t_in(n)
      integer :: iorder(n)
      t_in = (/ 1.0_dp, 3.0_dp, 2.0_dp /)
      t = t_in; iorder = (/ 1, 2, 3 /)
      call hpsolb(n, t, iorder, 1)
      call assert_min_at_n(t_in, t, iorder, n, "case_sift_down_stops_early")
      call assert_eq_int(iorder(n), 1, where="case_sift_down_stops_early iorder(n)")
   end subroutine case_sift_down_stops_early

   subroutine case_iheap0_min_at_first()
      ! Min already at position 1 on entry: heap-build is a no-op for k=2..n
      ! that just confirms structure (each newly-added element passes L55 F).
      integer, parameter :: n = 3
      real(dp) :: t(n), t_in(n)
      integer :: iorder(n)
      t_in = (/ 1.0_dp, 2.0_dp, 3.0_dp /)
      t = t_in; iorder = (/ 1, 2, 3 /)
      call hpsolb(n, t, iorder, 0)
      call assert_min_at_n(t_in, t, iorder, n, "case_iheap0_min_at_first")
      call assert_eq_int(iorder(n), 1, where="case_iheap0_min_at_first iorder(n)")
   end subroutine case_iheap0_min_at_first

end program test_hpsolb
