module test_assert
   implicit none
   private

   public :: assert_true
   public :: assert_eq_int
   public :: assert_eq_str
   public :: assert_close_real
   public :: assert_array_close

   integer, parameter, public :: dp = kind(0.0d0)
   real(dp), parameter :: default_rtol = 1.0e-12_dp

contains

   subroutine fail(msg, where)
      character(len=*), intent(in) :: msg
      character(len=*), intent(in), optional :: where
      if (present(where)) then
         write(0, '(A,": ",A)') trim(where), trim(msg)
      else
         write(0, '("FAIL: ",A)') trim(msg)
      end if
      error stop 1
   end subroutine fail

   subroutine assert_true(cond, where)
      logical, intent(in) :: cond
      character(len=*), intent(in), optional :: where
      if (.not. cond) call fail("assertion is false", where)
   end subroutine assert_true

   subroutine assert_eq_int(actual, expected, where)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in), optional :: where
      character(len=128) :: msg
      if (actual /= expected) then
         write(msg, '("expected ",I0,", got ",I0)') expected, actual
         call fail(trim(msg), where)
      end if
   end subroutine assert_eq_int

   subroutine assert_eq_str(actual, expected, where)
      character(len=*), intent(in) :: actual, expected
      character(len=*), intent(in), optional :: where
      character(len=512) :: msg
      if (actual /= expected) then
         write(msg, '("expected """,A,""", got """,A,"""")') &
            trim(expected), trim(actual)
         call fail(trim(msg), where)
      end if
   end subroutine assert_eq_str

   subroutine assert_close_real(actual, expected, rtol, where)
      real(dp), intent(in) :: actual, expected
      real(dp), intent(in), optional :: rtol
      character(len=*), intent(in), optional :: where
      real(dp) :: tol, dev
      character(len=256) :: msg
      if (present(rtol)) then
         tol = rtol
      else
         tol = default_rtol
      end if
      dev = abs(actual - expected) / (1.0_dp + abs(expected))
      if (dev >= tol) then
         write(msg, '("|a-b|/(1+|b|) = ",ES12.5," >= rtol = ",ES12.5, &
                     &" (a = ",ES23.16,", b = ",ES23.16,")")') &
            dev, tol, actual, expected
         call fail(trim(msg), where)
      end if
   end subroutine assert_close_real

   subroutine assert_array_close(actual, expected, rtol, where)
      real(dp), intent(in) :: actual(:), expected(:)
      real(dp), intent(in), optional :: rtol
      character(len=*), intent(in), optional :: where
      real(dp) :: tol, dev, worst
      integer :: i, iworst
      character(len=256) :: msg
      if (present(rtol)) then
         tol = rtol
      else
         tol = default_rtol
      end if
      if (size(actual) /= size(expected)) then
         write(msg, '("size mismatch: expected ",I0,", got ",I0)') &
            size(expected), size(actual)
         call fail(trim(msg), where)
      end if
      worst = 0.0_dp
      iworst = 0
      do i = 1, size(actual)
         dev = abs(actual(i) - expected(i)) / (1.0_dp + abs(expected(i)))
         if (dev > worst) then
            worst = dev
            iworst = i
         end if
      end do
      if (worst >= tol) then
         write(msg, '("worst element [",I0,"]: |a-b|/(1+|b|) = ",ES12.5, &
                     &" >= rtol = ",ES12.5," (a = ",ES23.16,", b = ",ES23.16,")")') &
            iworst, worst, tol, actual(iworst), expected(iworst)
         call fail(trim(msg), where)
      end if
   end subroutine assert_array_close

end module test_assert
