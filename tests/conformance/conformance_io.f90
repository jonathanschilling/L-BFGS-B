!> \file conformance_io.f90
!>
!> \brief Tiny key-value text I/O for the F77 conformance driver.
!>
!> Reads the input file produced by docs/spec/runner/conformance.py
!> (--engine fortran) into a fixed-size store of key/values, and writes
!> output in the same format. The format is:
!>
!>   sub <name>
!>   <type> <key> [<count>|<rows> <cols>] <values...>
!>   ...
!>   end
!>
!> with types: i, r, s, b, iv, rv, im, rm.
!> Matrices are column-major (Fortran natural).

module conformance_io
   implicit none
   private

   integer, parameter :: dp = kind(1.0d0)
   integer, parameter, public :: dp_kind = dp

   integer, parameter :: MAX_ENTRIES = 64
   integer, parameter :: KEY_LEN = 64
   integer, parameter :: VAL_LEN = 32768       ! large for big matrices

   type, public :: kv_store
      integer                       :: n = 0
      character(len=KEY_LEN)        :: keys(MAX_ENTRIES)  = ''
      character(len=VAL_LEN)        :: vals(MAX_ENTRIES)  = ''
   end type kv_store

   public :: read_kvs, get_int, get_real, get_str, get_bool
   public :: get_int_vec, get_real_vec, get_int_mat, get_real_mat
   public :: write_int, write_real, write_str, write_bool
   public :: write_int_vec, write_real_vec, write_int_mat, write_real_mat
   public :: write_end

contains

   !> Read key-value pairs from a unit until "end" is hit.
   subroutine read_kvs(unit, store)
      integer,         intent(in)    :: unit
      type(kv_store),  intent(inout) :: store
      character(len=VAL_LEN) :: line
      character(len=KEY_LEN) :: typ, key
      integer                :: ios, sp1, sp2

      store%n = 0
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:3) == 'end') exit
         if (line(1:4) == 'sub ') cycle

         ! Parse "type key rest" by finding two space separators.
         sp1 = index(line, ' ')
         if (sp1 == 0) cycle
         typ = line(1:sp1-1)
         line = adjustl(line(sp1+1:))
         sp2 = index(line, ' ')
         if (sp2 == 0) then
            key = line
            line = ''
         else
            key = line(1:sp2-1)
            line = line(sp2+1:)
         end if

         store%n = store%n + 1
         if (store%n > MAX_ENTRIES) then
            write(*,*) 'conformance_io: too many key-value entries (max ', MAX_ENTRIES, ')'
            stop 1
         end if
         store%keys(store%n) = trim(typ) // ':' // trim(key)
         store%vals(store%n) = line
      end do
   end subroutine read_kvs

   !> Look up the value string for "<typ>:<key>". Sets ios=-1 if not found.
   subroutine lookup(store, typ, key, val, ios)
      type(kv_store),         intent(in)  :: store
      character(len=*),       intent(in)  :: typ, key
      character(len=VAL_LEN), intent(out) :: val
      integer,                intent(out) :: ios
      character(len=KEY_LEN) :: needle
      integer :: i

      needle = trim(typ) // ':' // trim(key)
      ios = -1
      val = ''
      do i = 1, store%n
         if (trim(store%keys(i)) == trim(needle)) then
            val = store%vals(i)
            ios = 0
            return
         end if
      end do
   end subroutine lookup

   subroutine fail_missing(typ, key)
      character(len=*), intent(in) :: typ, key
      write(*,*) 'conformance_io: missing key ', trim(typ), ':', trim(key)
      stop 1
   end subroutine fail_missing

   subroutine get_int(store, key, val)
      type(kv_store), intent(in)  :: store
      character(len=*), intent(in)  :: key
      integer,        intent(out) :: val
      character(len=VAL_LEN) :: s
      integer :: ios
      call lookup(store, 'i', key, s, ios)
      if (ios /= 0) call fail_missing('i', key)
      read(s, *) val
   end subroutine get_int

   subroutine get_real(store, key, val)
      type(kv_store), intent(in)  :: store
      character(len=*), intent(in)  :: key
      real(dp),       intent(out) :: val
      character(len=VAL_LEN) :: s
      integer :: ios
      call lookup(store, 'r', key, s, ios)
      if (ios /= 0) call fail_missing('r', key)
      read(s, *) val
   end subroutine get_real

   subroutine get_str(store, key, val)
      type(kv_store), intent(in)  :: store
      character(len=*), intent(in)  :: key
      character(len=*), intent(out) :: val
      character(len=VAL_LEN) :: s
      integer :: ios
      call lookup(store, 's', key, s, ios)
      if (ios /= 0) call fail_missing('s', key)
      val = s
   end subroutine get_str

   subroutine get_bool(store, key, val)
      type(kv_store), intent(in)  :: store
      character(len=*), intent(in)  :: key
      logical,        intent(out) :: val
      character(len=VAL_LEN) :: s
      integer :: ios
      call lookup(store, 'b', key, s, ios)
      if (ios /= 0) call fail_missing('b', key)
      s = adjustl(s)
      val = (s(1:1) == 'T' .or. s(1:1) == 't')
   end subroutine get_bool

   subroutine get_int_vec(store, key, vec, n)
      type(kv_store),     intent(in)  :: store
      character(len=*),    intent(in)  :: key
      integer,             intent(out) :: n
      integer, allocatable, intent(out) :: vec(:)
      character(len=VAL_LEN) :: s
      integer :: ios, i
      call lookup(store, 'iv', key, s, ios)
      if (ios /= 0) call fail_missing('iv', key)
      read(s, *) n
      allocate(vec(n))
      ! Re-parse: skip past the leading n.
      read(s, *) i, vec(:)         ! 'i' is a throwaway holding n
   end subroutine get_int_vec

   subroutine get_real_vec(store, key, vec, n)
      type(kv_store),     intent(in)  :: store
      character(len=*),    intent(in)  :: key
      integer,             intent(out) :: n
      real(dp), allocatable, intent(out) :: vec(:)
      character(len=VAL_LEN) :: s
      integer :: ios, dummy_n
      call lookup(store, 'rv', key, s, ios)
      if (ios /= 0) call fail_missing('rv', key)
      read(s, *) dummy_n
      n = dummy_n
      allocate(vec(n))
      read(s, *) dummy_n, vec(:)
   end subroutine get_real_vec

   subroutine get_int_mat(store, key, mat, rows, cols)
      type(kv_store),       intent(in)  :: store
      character(len=*),      intent(in)  :: key
      integer,               intent(out) :: rows, cols
      integer, allocatable,  intent(out) :: mat(:,:)
      character(len=VAL_LEN) :: s
      integer :: ios, dr, dc
      call lookup(store, 'im', key, s, ios)
      if (ios /= 0) call fail_missing('im', key)
      read(s, *) dr, dc
      rows = dr; cols = dc
      allocate(mat(rows, cols))
      read(s, *) dr, dc, mat(:,:)
   end subroutine get_int_mat

   subroutine get_real_mat(store, key, mat, rows, cols)
      type(kv_store),         intent(in)  :: store
      character(len=*),        intent(in)  :: key
      integer,                 intent(out) :: rows, cols
      real(dp), allocatable,   intent(out) :: mat(:,:)
      character(len=VAL_LEN) :: s
      integer :: ios, dr, dc
      call lookup(store, 'rm', key, s, ios)
      if (ios /= 0) call fail_missing('rm', key)
      read(s, *) dr, dc
      rows = dr; cols = dc
      allocate(mat(rows, cols))
      read(s, *) dr, dc, mat(:,:)
   end subroutine get_real_mat

   ! ----------------------- writers --------------------------

   subroutine write_int(unit, key, val)
      integer,           intent(in) :: unit, val
      character(len=*),  intent(in) :: key
      write(unit, '("i ",A," ",I0)') trim(key), val
   end subroutine write_int

   subroutine write_real(unit, key, val)
      integer,           intent(in) :: unit
      character(len=*),  intent(in) :: key
      real(dp),          intent(in) :: val
      ! ES23.16 writes 17 sig figs; sufficient for IEEE-754 round-trip.
      write(unit, '("r ",A," ",ES26.17E3)') trim(key), val
   end subroutine write_real

   subroutine write_str(unit, key, val)
      integer,           intent(in) :: unit
      character(len=*),  intent(in) :: key, val
      write(unit, '("s ",A," ",A)') trim(key), trim(val)
   end subroutine write_str

   subroutine write_bool(unit, key, val)
      integer,           intent(in) :: unit
      character(len=*),  intent(in) :: key
      logical,           intent(in) :: val
      character(len=1) :: c
      c = merge('T', 'F', val)
      write(unit, '("b ",A," ",A)') trim(key), c
   end subroutine write_bool

   subroutine write_int_vec(unit, key, vec)
      integer,           intent(in) :: unit
      character(len=*),  intent(in) :: key
      integer,           intent(in) :: vec(:)
      integer :: i
      write(unit, '("iv ",A," ",I0)', advance='no') trim(key), size(vec)
      do i = 1, size(vec)
         write(unit, '(" ",I0)', advance='no') vec(i)
      end do
      write(unit, '(A)') ''
   end subroutine write_int_vec

   subroutine write_real_vec(unit, key, vec)
      integer,           intent(in) :: unit
      character(len=*),  intent(in) :: key
      real(dp),          intent(in) :: vec(:)
      integer :: i
      write(unit, '("rv ",A," ",I0)', advance='no') trim(key), size(vec)
      do i = 1, size(vec)
         write(unit, '(" ",ES26.17E3)', advance='no') vec(i)
      end do
      write(unit, '(A)') ''
   end subroutine write_real_vec

   subroutine write_int_mat(unit, key, mat)
      integer,           intent(in) :: unit
      character(len=*),  intent(in) :: key
      integer,           intent(in) :: mat(:,:)
      integer :: i, j
      write(unit, '("im ",A," ",I0," ",I0)', advance='no') &
            trim(key), size(mat,1), size(mat,2)
      do j = 1, size(mat, 2)
         do i = 1, size(mat, 1)
            write(unit, '(" ",I0)', advance='no') mat(i, j)
         end do
      end do
      write(unit, '(A)') ''
   end subroutine write_int_mat

   subroutine write_real_mat(unit, key, mat)
      integer,           intent(in) :: unit
      character(len=*),  intent(in) :: key
      real(dp),          intent(in) :: mat(:,:)
      integer :: i, j
      write(unit, '("rm ",A," ",I0," ",I0)', advance='no') &
            trim(key), size(mat,1), size(mat,2)
      do j = 1, size(mat, 2)
         do i = 1, size(mat, 1)
            write(unit, '(" ",ES26.17E3)', advance='no') mat(i, j)
         end do
      end do
      write(unit, '(A)') ''
   end subroutine write_real_mat

   subroutine write_end(unit)
      integer, intent(in) :: unit
      write(unit, '("end")')
   end subroutine write_end

end module conformance_io
