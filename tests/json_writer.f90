!> \file json_writer.f90
!> Test instrumentation for the L-BFGS-B drivers: writes the aggregate
!> end-of-run state (final task, final f, final |proj g|) to a JSON file.
!> Activated by setting the LBFGSB_JSON_OUTPUT environment variable.
!>
!> Per-iteration data is intentionally not captured: the L-BFGS-B
!> trajectory is sensitive to FP rounding and diverges across toolchains
!> after enough iterations. Aggregate quantities are robust enough to
!> compare across compiler/BLAS versions while still catching real
!> algorithmic regressions (different stopping reason, different f).

module json_writer
   implicit none
   private
   public :: json_write_aggregate

   integer, parameter :: dp = kind(0.0d0)

contains

   subroutine json_write_aggregate(path, task, f, projg)
      character(len=*), intent(in) :: path, task
      real(dp),         intent(in) :: f, projg
      integer :: u

      open(newunit=u, file=path, status='replace', action='write')
      write(u,'(A)')           '{'
      write(u,'(A,A,A)')       '  "final_task": "',  trim(task), '",'
      write(u,'(A,ES23.16,A)') '  "final_f": ',      f,          ','
      write(u,'(A,ES23.16)')   '  "final_projg": ',  projg
      write(u,'(A)')           '}'
      close(u)
   end subroutine json_write_aggregate

end module json_writer
