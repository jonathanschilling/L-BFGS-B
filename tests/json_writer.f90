!> \file json_writer.f90
!> Test instrumentation for the L-BFGS-B drivers: writes per-iteration
!> optimization data to a JSON file in full double precision.
!> Activated by setting the LBFGSB_JSON_OUTPUT environment variable.

module json_writer
   implicit none
   private
   public :: json_open, json_iter, json_close

   integer, parameter :: dp = kind(0.0d0)
   integer, parameter :: U  = 99
   logical, save :: pending_close = .false.

contains

   subroutine json_open(path)
      character(len=*), intent(in) :: path
      open(unit=U, file=path, status='replace', action='write')
      write(U,'(A)') '{'
      write(U,'(A)') '  "iterations": ['
      pending_close = .false.
   end subroutine json_open

   subroutine json_iter(n, x, g, f, projg, isave)
      integer,  intent(in) :: n
      real(dp), intent(in) :: x(n), g(n), f, projg
      integer,  intent(in) :: isave(*)
      integer :: i

      if (pending_close) then
         write(U,'(A)') '    },'
      end if
      pending_close = .true.

      write(U,'(A)') '    {'
      write(U,'(A,I0,A)')      '      "iter": ',  isave(30), ','
      write(U,'(A,I0,A)')      '      "nfg": ',   isave(34), ','
      write(U,'(A,ES23.16,A)') '      "f": ',     f,         ','
      write(U,'(A,ES23.16,A)') '      "projg": ', projg,     ','
      call write_array('x', x, n, .true.)
      call write_array('g', g, n, .false.)
   end subroutine json_iter

   subroutine json_close(n, x, f, task)
      integer,  intent(in) :: n
      real(dp), intent(in) :: x(n), f
      character(len=*), intent(in) :: task

      if (pending_close) then
         write(U,'(A)') '    }'
         pending_close = .false.
      end if
      write(U,'(A)') '  ],'
      write(U,'(A,A,A)')       '  "final_task": "', trim(task), '",'
      write(U,'(A,ES23.16,A)') '  "final_f": ',     f,          ','
      call write_array_top('final_x', x, n, .false.)
      write(U,'(A)') '}'
      close(U)
   end subroutine json_close

   subroutine write_array(key, v, n, trailing_comma)
      character(len=*), intent(in) :: key
      integer,  intent(in) :: n
      real(dp), intent(in) :: v(n)
      logical,  intent(in) :: trailing_comma
      integer :: i

      write(U,'(A,A,A)', advance='no') '      "', key, '": ['
      do i = 1, n
         if (i > 1) write(U,'(A)', advance='no') ', '
         write(U,'(ES23.16)', advance='no') v(i)
      end do
      if (trailing_comma) then
         write(U,'(A)') '],'
      else
         write(U,'(A)') ']'
      end if
   end subroutine write_array

   subroutine write_array_top(key, v, n, trailing_comma)
      character(len=*), intent(in) :: key
      integer,  intent(in) :: n
      real(dp), intent(in) :: v(n)
      logical,  intent(in) :: trailing_comma
      integer :: i

      write(U,'(A,A,A)', advance='no') '  "', key, '": ['
      do i = 1, n
         if (i > 1) write(U,'(A)', advance='no') ', '
         write(U,'(ES23.16)', advance='no') v(i)
      end do
      if (trailing_comma) then
         write(U,'(A)') '],'
      else
         write(U,'(A)') ']'
      end if
   end subroutine write_array_top

end module json_writer
