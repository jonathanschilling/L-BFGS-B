program test_freev
   use test_assert, only: assert_true, assert_eq_int
   implicit none

   call case_iter_zero_mixed_iwhere()
   call case_iter_pos_unconstrained()
   call case_iter_pos_no_changes_updatd_true()
   call case_iter_pos_one_leaves()
   call case_iter_pos_one_enters()
   call case_iter_pos_both_events_iprint_100()

   write(*, '("test_freev: PASS")')

contains

   subroutine case_iter_zero_mixed_iwhere()
      ! iter=0: skips the count block (L55 F via iter check).
      ! Second loop (L86-) always runs and partitions index by current iwhere.
      ! iwhere = (-1, 2) : var 1 free, var 2 bound -> nfree=1, index=(1,2).
      integer, parameter :: n = 2
      integer :: index(n), indx2(n), iwhere(n)
      integer :: nfree, nenter, ileave
      logical :: wrk, updatd, cnstnd
      iwhere = (/ -1, 2 /)
      index = 0; indx2 = 0
      nfree = 0; nenter = 99; ileave = 99
      updatd = .false.; cnstnd = .true.
      call freev(n, nfree, index, nenter, ileave, indx2, &
                 iwhere, wrk, updatd, cnstnd, -1, 0)
      call assert_eq_int(nfree, 1, where="case_iter_zero_mixed_iwhere nfree")
      call assert_eq_int(index(1), 1, where="case_iter_zero_mixed_iwhere index(1)")
      call assert_eq_int(index(2), 2, where="case_iter_zero_mixed_iwhere index(2)")
      call assert_eq_int(nenter, 0, where="case_iter_zero_mixed_iwhere nenter")
      call assert_eq_int(ileave, n+1, where="case_iter_zero_mixed_iwhere ileave")
      call assert_true(.not. wrk, "case_iter_zero_mixed_iwhere wrk")
   end subroutine case_iter_zero_mixed_iwhere

   subroutine case_iter_pos_unconstrained()
      ! iter>0 but cnstnd=false: also skips the count block (L55 F via cnstnd).
      integer, parameter :: n = 2
      integer :: index(n), indx2(n), iwhere(n)
      integer :: nfree, nenter, ileave
      logical :: wrk, updatd, cnstnd
      iwhere = (/ -1, -1 /)
      index = 0; indx2 = 0
      nfree = 0; nenter = 99; ileave = 99
      updatd = .false.; cnstnd = .false.
      call freev(n, nfree, index, nenter, ileave, indx2, &
                 iwhere, wrk, updatd, cnstnd, -1, 1)
      call assert_eq_int(nfree, 2, where="case_iter_pos_unconstrained nfree")
      call assert_true(.not. wrk, "case_iter_pos_unconstrained wrk")
   end subroutine case_iter_pos_unconstrained

   subroutine case_iter_pos_no_changes_updatd_true()
      ! iter>0, cnstnd=true: enters count block. Previous active set matches
      ! current iwhere -> no leaving, no entering. updatd=true forces wrk=true
      ! via the OR-clause in L82.
      ! Previous: nfree=1, index(1)=2, index(2)=1 (var 2 free, var 1 bound).
      ! Current : iwhere = (1, -1) -> var 2 still free, var 1 still bound.
      integer, parameter :: n = 2
      integer :: index(n), indx2(n), iwhere(n)
      integer :: nfree, nenter, ileave
      logical :: wrk, updatd, cnstnd
      iwhere = (/ 1, -1 /)
      index = (/ 2, 1 /); indx2 = 0
      nfree = 1; nenter = 99; ileave = 99
      updatd = .true.; cnstnd = .true.
      call freev(n, nfree, index, nenter, ileave, indx2, &
                 iwhere, wrk, updatd, cnstnd, -1, 1)
      call assert_eq_int(nenter, 0, where="case_iter_pos_no_changes_updatd_true nenter")
      call assert_eq_int(ileave, n+1, where="case_iter_pos_no_changes_updatd_true ileave")
      call assert_true(wrk, "case_iter_pos_no_changes_updatd_true wrk")
   end subroutine case_iter_pos_no_changes_updatd_true

   subroutine case_iter_pos_one_leaves()
      ! Var 2 was free, now bound -> "leaves" the free set.
      ! Previous: nfree=2, index=(1,2). Current: iwhere=(-1, 1).
      integer, parameter :: n = 2
      integer :: index(n), indx2(n), iwhere(n)
      integer :: nfree, nenter, ileave
      logical :: wrk, updatd, cnstnd
      iwhere = (/ -1, 1 /)
      index = (/ 1, 2 /); indx2 = 0
      nfree = 2; nenter = 99; ileave = 99
      updatd = .false.; cnstnd = .true.
      call freev(n, nfree, index, nenter, ileave, indx2, &
                 iwhere, wrk, updatd, cnstnd, -1, 1)
      call assert_eq_int(ileave, 2, where="case_iter_pos_one_leaves ileave")
      call assert_eq_int(indx2(2), 2, where="case_iter_pos_one_leaves indx2(2)")
      call assert_eq_int(nenter, 0, where="case_iter_pos_one_leaves nenter")
      call assert_true(wrk, "case_iter_pos_one_leaves wrk")
   end subroutine case_iter_pos_one_leaves

   subroutine case_iter_pos_one_enters()
      ! Var 1 was bound, now free -> "enters" the free set.
      ! Previous: nfree=1, index=(2, 1). Current: iwhere=(-1, -1).
      integer, parameter :: n = 2
      integer :: index(n), indx2(n), iwhere(n)
      integer :: nfree, nenter, ileave
      logical :: wrk, updatd, cnstnd
      iwhere = (/ -1, -1 /)
      index = (/ 2, 1 /); indx2 = 0
      nfree = 1; nenter = 99; ileave = 99
      updatd = .false.; cnstnd = .true.
      call freev(n, nfree, index, nenter, ileave, indx2, &
                 iwhere, wrk, updatd, cnstnd, -1, 1)
      call assert_eq_int(nenter, 1, where="case_iter_pos_one_enters nenter")
      call assert_eq_int(indx2(1), 1, where="case_iter_pos_one_enters indx2(1)")
      call assert_eq_int(ileave, n+1, where="case_iter_pos_one_enters ileave")
      call assert_true(wrk, "case_iter_pos_one_enters wrk")
   end subroutine case_iter_pos_one_enters

   subroutine case_iter_pos_both_events_iprint_100()
      ! Var 1 was free now bound (leaves), var 2 was bound now free (enters).
      ! iprint=100 exercises the "Variable X enters/leaves" diagnostic prints
      ! and the iprint>=99 summary prints.
      ! Previous: nfree=1, index=(1, 2). Current: iwhere=(2, -1).
      integer, parameter :: n = 2
      integer :: index(n), indx2(n), iwhere(n)
      integer :: nfree, nenter, ileave
      logical :: wrk, updatd, cnstnd
      iwhere = (/ 2, -1 /)
      index = (/ 1, 2 /); indx2 = 0
      nfree = 1; nenter = 99; ileave = 99
      updatd = .false.; cnstnd = .true.
      call freev(n, nfree, index, nenter, ileave, indx2, &
                 iwhere, wrk, updatd, cnstnd, 100, 1)
      call assert_eq_int(nenter, 1, where="case_iter_pos_both_events nenter")
      call assert_true(ileave .lt. n+1, "case_iter_pos_both_events ileave decreased")
   end subroutine case_iter_pos_both_events_iprint_100

end program test_freev
