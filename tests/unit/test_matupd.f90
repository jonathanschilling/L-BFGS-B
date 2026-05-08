program test_matupd
   use test_assert, only: dp, assert_true, assert_eq_int, assert_eq_str, assert_close_real, assert_array_close
   implicit none

   call case_walkthrough_three_updates_with_overflow()
   call case_stp_equals_one()

   write(*, '("test_matupd: PASS")')

contains

   subroutine case_walkthrough_three_updates_with_overflow()
      ! n=2, m=2 walkthrough: iupdat=1, 2, 3. The third update wraps the
      ! cyclic buffer (iupdat > m), exercising the shift-old-info loop and
      ! the head/itail wrap. Throughout, stp != 1 so we hit L112's else.
      integer, parameter :: n = 2, m = 2
      integer  :: itail, iupdat, col, head
      real(dp) :: ws(n,m), wy(n,m), sy(m,m), ss(m,m)
      real(dp) :: d(n), r(n), theta, rr, dr, stp, dtd

      ws = 0.0_dp; wy = 0.0_dp; sy = 0.0_dp; ss = 0.0_dp
      itail = 0; col = 0; head = 1
      stp = 0.5_dp

      ! ---- iupdat=1: first update, slot 1 ----
      iupdat = 1
      d = (/ 1.0_dp, 0.0_dp /)
      r = (/ 1.0_dp, 1.0_dp /)
      rr = 2.0_dp; dr = 1.0_dp; dtd = 1.0_dp
      call matupd(n, m, ws, wy, sy, ss, d, r, itail, &
                  iupdat, col, head, theta, rr, dr, stp, dtd)
      call assert_eq_int(col,   1, where="iupdat=1 col")
      call assert_eq_int(itail, 1, where="iupdat=1 itail")
      call assert_eq_int(head,  1, where="iupdat=1 head")
      call assert_close_real(theta, 2.0_dp, where="iupdat=1 theta")
      call assert_close_real(ws(1,1), 1.0_dp, where="iupdat=1 ws(1,1)")
      call assert_close_real(ws(2,1), 0.0_dp, where="iupdat=1 ws(2,1)")
      call assert_close_real(wy(1,1), 1.0_dp, where="iupdat=1 wy(1,1)")
      call assert_close_real(wy(2,1), 1.0_dp, where="iupdat=1 wy(2,1)")
      ! stp=0.5, dtd=1 -> ss(1,1) = 0.25.
      call assert_close_real(ss(1,1), 0.25_dp, where="iupdat=1 ss(1,1)")
      call assert_close_real(sy(1,1), 1.0_dp,  where="iupdat=1 sy(1,1)")

      ! ---- iupdat=2: second update, slot 2 ----
      iupdat = 2
      d = (/ 0.0_dp, 1.0_dp /)
      r = (/ 0.0_dp, 1.0_dp /)
      rr = 1.0_dp; dr = 1.0_dp; dtd = 1.0_dp
      stp = 0.5_dp
      call matupd(n, m, ws, wy, sy, ss, d, r, itail, &
                  iupdat, col, head, theta, rr, dr, stp, dtd)
      call assert_eq_int(col,   2, where="iupdat=2 col")
      call assert_eq_int(itail, 2, where="iupdat=2 itail")
      call assert_eq_int(head,  1, where="iupdat=2 head")
      call assert_close_real(theta, 1.0_dp, where="iupdat=2 theta")
      ! ws(:,1) preserved, ws(:,2) = new d.
      call assert_close_real(ws(1,1), 1.0_dp, where="iupdat=2 ws(1,1)")
      call assert_close_real(ws(1,2), 0.0_dp, where="iupdat=2 ws(1,2)")
      call assert_close_real(ws(2,2), 1.0_dp, where="iupdat=2 ws(2,2)")
      ! Loop 51 fills sy(2,1) = d'wy(:,1) = (0,1)'(1,1) = 1.
      call assert_close_real(sy(2,1), 1.0_dp,  where="iupdat=2 sy(2,1)")
      ! Loop 51 fills ss(1,2) = ws(:,1)'d = (1,0)'(0,1) = 0.
      call assert_close_real(ss(1,2), 0.0_dp,  where="iupdat=2 ss(1,2)")
      ! Final-row diagonal: ss(2,2) = stp^2 * dtd = 0.25; sy(2,2) = dr = 1.
      call assert_close_real(ss(2,2), 0.25_dp, where="iupdat=2 ss(2,2)")
      call assert_close_real(sy(2,2), 1.0_dp,  where="iupdat=2 sy(2,2)")

      ! ---- iupdat=3: overflow. itail wraps to 1, head wraps to 2. ----
      iupdat = 3
      d = (/ 1.0_dp, 1.0_dp /)
      r = (/ 2.0_dp, 0.0_dp /)
      rr = 4.0_dp; dr = 2.0_dp; dtd = 2.0_dp
      stp = 0.5_dp
      call matupd(n, m, ws, wy, sy, ss, d, r, itail, &
                  iupdat, col, head, theta, rr, dr, stp, dtd)
      call assert_eq_int(col,   2, where="iupdat=3 col")
      call assert_eq_int(itail, 1, where="iupdat=3 itail")
      call assert_eq_int(head,  2, where="iupdat=3 head")
      call assert_close_real(theta, 2.0_dp, where="iupdat=3 theta")
      ! Old slot 1 (iupdat=1) overwritten; ws(:,1) = new d, wy(:,1) = new r.
      call assert_close_real(ws(1,1), 1.0_dp, where="iupdat=3 ws(1,1)")
      call assert_close_real(ws(2,1), 1.0_dp, where="iupdat=3 ws(2,1)")
      call assert_close_real(wy(1,1), 2.0_dp, where="iupdat=3 wy(1,1)")
      call assert_close_real(wy(2,1), 0.0_dp, where="iupdat=3 wy(2,1)")
      ! Slot 2 from iupdat=2 still preserved.
      call assert_close_real(ws(2,2), 1.0_dp, where="iupdat=3 ws(2,2) preserved")
      ! Shift loop pulled old ss(2,2)=0.25 into ss(1,1), sy(2,2)=1 into sy(1,1).
      call assert_close_real(ss(1,1), 0.25_dp, where="iupdat=3 ss(1,1) shifted")
      call assert_close_real(sy(1,1), 1.0_dp,  where="iupdat=3 sy(1,1) shifted")
      ! New row/col: pointr=head=2.
      ! sy(2,1) = d'wy(:,2) = (1,1)'(0,1) = 1.
      ! ss(1,2) = ws(:,2)'d = (0,1)'(1,1) = 1.
      call assert_close_real(sy(2,1), 1.0_dp, where="iupdat=3 sy(2,1) new")
      call assert_close_real(ss(1,2), 1.0_dp, where="iupdat=3 ss(1,2) new")
      ! Diagonal: stp=0.5, dtd=2 -> ss(2,2) = 0.25 * 2 = 0.5; sy(2,2) = dr = 2.
      call assert_close_real(ss(2,2), 0.5_dp, where="iupdat=3 ss(2,2)")
      call assert_close_real(sy(2,2), 2.0_dp, where="iupdat=3 sy(2,2)")
   end subroutine case_walkthrough_three_updates_with_overflow

   subroutine case_stp_equals_one()
      ! Single update with stp = 1.0 exactly: covers the L112 true branch
      ! (ss(col,col) = dtd, no stp^2 multiply).
      integer, parameter :: n = 2, m = 2
      integer  :: itail, iupdat, col, head
      real(dp) :: ws(n,m), wy(n,m), sy(m,m), ss(m,m)
      real(dp) :: d(n), r(n), theta, rr, dr, stp, dtd
      ws = 0.0_dp; wy = 0.0_dp; sy = 0.0_dp; ss = 0.0_dp
      itail = 0; col = 0; head = 1
      iupdat = 1
      d = (/ 2.0_dp, 0.0_dp /)
      r = (/ 1.0_dp, 1.0_dp /)
      rr = 2.0_dp; dr = 2.0_dp; dtd = 4.0_dp
      stp = 1.0_dp
      call matupd(n, m, ws, wy, sy, ss, d, r, itail, &
                  iupdat, col, head, theta, rr, dr, stp, dtd)
      call assert_close_real(ss(1,1), 4.0_dp, where="case_stp_equals_one ss(1,1)")
   end subroutine case_stp_equals_one

end program test_matupd
