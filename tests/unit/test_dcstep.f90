program test_dcstep
   use test_assert
   implicit none

   ! dcstep is the More-Thuente safeguarded cubic-step generator. Its outer
   ! IF/ELSEIF chain has four cases:
   !   Case 1:  fp > fx                                  (function increased)
   !   Case 2:  fp <= fx, sgnd < 0                       (derivatives flipped)
   !   Case 3:  fp <= fx, sgnd >= 0, |dpval| < |dx|         (descent, |g| shrunk)
   !   Case 4:  fp <= fx, sgnd >= 0, |dpval| >= |dx|        (descent, |g| held)
   ! Each is accompanied by an interval-update block at the bottom of the
   ! routine (sty/stx/fy/fx/dy/dx). Test cases below each pin one branch
   ! of the case decomposition and assert the structural updates plus
   ! that the new stp lands in [stpmin, stpmax].

   call case1_fp_higher()
   call case1_stp_below_stx()
   call case2_opposite_sign_derivs()
   call case3_decreasing_grad_brackt_true()
   call case3_decreasing_grad_brackt_false()
   call case3_stp_above_stx_no_bracket()
   call case4_brackt_true()
   call case4_brackt_false_stp_above()
   call case4_brackt_false_stp_below()

   write(*, '("test_dcstep: PASS")')

contains

   subroutine assert_in_range(stp, stpmin, stpmax, where)
      real(dp), intent(in) :: stp, stpmin, stpmax
      character(len=*), intent(in) :: where
      call assert_true(stp >= stpmin .and. stp <= stpmax, &
                       where // " stp in [stpmin, stpmax]")
   end subroutine assert_in_range

   subroutine case1_fp_higher()
      ! fp > fx with stp > stx: Case 1 path. Bracket gets set, sty <- stp.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp; fx = -1.0_dp; dx = -1.0_dp
      sty = 0.0_dp; fy =  0.0_dp; dy =  0.0_dp
      stp = 1.0_dp; fp =  1.0_dp; dpval =  4.0_dp
      brackt = .false.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_true(brackt, "case1_fp_higher brackt")
      call assert_close_real(sty, 1.0_dp, where="case1_fp_higher sty<-stp")
      call assert_close_real(fy,  1.0_dp, where="case1_fp_higher fy<-fp")
      call assert_close_real(dy,  4.0_dp, where="case1_fp_higher dy<-dpval")
      ! stx, fx, dx unchanged.
      call assert_close_real(stx, 0.0_dp,  where="case1_fp_higher stx unchanged")
      call assert_in_range(stp, stpmin, stpmax, "case1_fp_higher")
   end subroutine case1_fp_higher

   subroutine case1_stp_below_stx()
      ! Case 1 with stp < stx: triggers the gamma-sign flip "if (stp < stx)
      ! gamma = -gamma". Take a quadratic translated so f(0)=0.16 with
      ! f'(0)=-0.8; choose stp=-1 with fp manufactured to be > fx=0.16.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp;  fx = 0.16_dp; dx = -0.8_dp
      sty = 0.0_dp;  fy = 0.0_dp;  dy =  0.0_dp
      stp = -1.0_dp; fp = 5.0_dp;  dpval = -3.0_dp
      brackt = .false.; stpmin = -10.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_true(brackt, "case1_stp_below_stx brackt")
      call assert_close_real(sty, -1.0_dp, where="case1_stp_below_stx sty<-stp")
   end subroutine case1_stp_below_stx

   subroutine case2_opposite_sign_derivs()
      ! fp <= fx with sgnd < 0: Case 2. sty <- old stx, stx <- stp.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp; fx = 1.0_dp;  dx = -1.0_dp
      sty = 0.0_dp; fy = 0.0_dp;  dy =  0.0_dp
      stp = 1.0_dp; fp = 0.5_dp;  dpval =  2.0_dp
      brackt = .false.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_true(brackt, "case2 brackt")
      ! Update block: fp <= fx, sgnd<0 -> sty<-old stx; then stx<-stp.
      call assert_close_real(sty, 0.0_dp, where="case2 sty<-old stx")
      call assert_close_real(fy,  1.0_dp, where="case2 fy<-old fx")
      call assert_close_real(dy, -1.0_dp, where="case2 dy<-old dx")
      call assert_close_real(stx, 1.0_dp, where="case2 stx<-stp")
      call assert_close_real(fx,  0.5_dp, where="case2 fx<-fp")
      call assert_close_real(dx,  2.0_dp, where="case2 dx<-dpval")
   end subroutine case2_opposite_sign_derivs

   subroutine case3_decreasing_grad_brackt_true()
      ! Case 3 with brackt=T: |dpval| < |dx|, sgnd >= 0, fp <= fx.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp; fx = 1.0_dp;  dx = -2.0_dp
      sty = 5.0_dp; fy = 2.0_dp;  dy =  1.0_dp
      stp = 1.0_dp; fp = 0.0_dp;  dpval = -0.5_dp
      brackt = .true.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      ! Case 3 does not toggle brackt; it stays true.
      call assert_true(brackt, "case3_brackt_true brackt unchanged")
      ! Update block: fp <= fx, sgnd >= 0 -> stx <- stp, sty unchanged.
      call assert_close_real(stx, 1.0_dp,  where="case3_brackt_true stx<-stp")
      call assert_close_real(sty, 5.0_dp,  where="case3_brackt_true sty unchanged")
      call assert_in_range(stp, stpmin, stpmax, "case3_brackt_true")
   end subroutine case3_decreasing_grad_brackt_true

   subroutine case3_decreasing_grad_brackt_false()
      ! Case 3 with brackt=F: cubic step compared against secant; clamped to
      ! [stpmin, stpmax] explicitly at the end of the case-3 else block.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp; fx = 1.0_dp;  dx = -2.0_dp
      sty = 0.0_dp; fy = 0.0_dp;  dy =  0.0_dp
      stp = 1.0_dp; fp = 0.0_dp;  dpval = -0.5_dp
      brackt = .false.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_in_range(stp, stpmin, stpmax, "case3_brackt_false")
      call assert_close_real(stx, 1.0_dp, where="case3_brackt_false stx<-stp")
   end subroutine case3_decreasing_grad_brackt_false

   subroutine case3_stp_above_stx_no_bracket()
      ! Case 3 with stp > stx and the cubic-step branch where r >= 0
      ! (forces the "stpc = stpmax" branch when stp > stx).
      ! Choose dx, dpval same sign and |dpval| < |dx|, but with q yielding r >= 0.
      ! Achievable when theta is small relative to gamma.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp; fx = 1.0_dp;  dx = -2.0_dp
      sty = 0.0_dp; fy = 0.0_dp;  dy =  0.0_dp
      stp = 1.0_dp; fp = -1.0_dp; dpval = -1.0_dp
      brackt = .false.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_in_range(stp, stpmin, stpmax, "case3_stp_above_stx_no_bracket")
   end subroutine case3_stp_above_stx_no_bracket

   subroutine case4_brackt_true()
      ! Case 4 with brackt=T: cubic step computed using sty data instead of stx.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp; fx = 2.0_dp;  dx = -1.0_dp
      sty = 5.0_dp; fy = 1.0_dp;  dy = -2.0_dp
      stp = 1.0_dp; fp = 1.5_dp;  dpval = -1.5_dp
      brackt = .true.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_close_real(stx, 1.0_dp, where="case4_brackt_true stx<-stp")
      call assert_in_range(stp, stpmin, stpmax, "case4_brackt_true")
   end subroutine case4_brackt_true

   subroutine case4_brackt_false_stp_above()
      ! Case 4 brackt=F with stp > stx: stpf = stpmax.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 0.0_dp; fx = 1.0_dp; dx = -1.0_dp
      sty = 0.0_dp; fy = 0.0_dp; dy =  0.0_dp
      stp = 1.0_dp; fp = 0.5_dp; dpval = -2.0_dp
      brackt = .false.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_close_real(stp, stpmax, where="case4_above stp=stpmax")
      call assert_close_real(stx, 1.0_dp, where="case4_above stx<-stp_in")
   end subroutine case4_brackt_false_stp_above

   subroutine case4_brackt_false_stp_below()
      ! Case 4 brackt=F with stp < stx: stpf = stpmin.
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dpval, stpmin, stpmax
      logical  :: brackt
      stx = 2.0_dp; fx = 1.0_dp; dx =  1.0_dp
      sty = 0.0_dp; fy = 0.0_dp; dy =  0.0_dp
      stp = 1.0_dp; fp = 0.5_dp; dpval =  2.0_dp
      brackt = .false.; stpmin = 0.0_dp; stpmax = 10.0_dp
      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dpval, brackt, stpmin, stpmax)
      call assert_close_real(stp, stpmin, where="case4_below stp=stpmin")
      call assert_close_real(stx, 1.0_dp, where="case4_below stx<-stp_in")
   end subroutine case4_brackt_false_stp_below

end program test_dcstep
