!> \file conformance_driver.f90
!>
!> \brief Driver that runs F77 subroutines on JSON-derived test vectors.
!>
!> Reads a structured text input file (one key-value entry per line; see
!> conformance_io.f90 for the format), dispatches by subroutine name to
!> the appropriate F77 call, and writes outputs in the same format.
!>
!> Used by docs/spec/runner/conformance.py with --engine fortran. The
!> Python runner converts JSON test vectors to the text format, invokes
!> this driver per case, and parses the output.
!>
!> Build target: conformance_driver (built when ctest is enabled). Links
!> against liblbfgsb so it can call any F77 routine in src/.
!>
!> Usage:
!>   conformance_driver <input_file> <output_file>

program conformance_driver
   use conformance_io
   implicit none

   integer, parameter :: dp = kind(1.0d0)

   character(len=1024) :: input_path, output_path
   character(len=1024) :: line
   character(len=64)   :: sub_name
   integer             :: in_unit, out_unit, ios
   type(kv_store)      :: store

   if (command_argument_count() /= 2) then
      write(*,*) 'usage: conformance_driver <input_file> <output_file>'
      stop 2
   end if

   call get_command_argument(1, input_path)
   call get_command_argument(2, output_path)

   open(newunit=in_unit, file=trim(input_path), action='read', status='old')
   open(newunit=out_unit, file=trim(output_path), action='write', status='replace')

   ! First non-blank line is "sub <name>".
   do
      read(in_unit, '(A)', iostat=ios) line
      if (ios /= 0) then
         write(*,*) 'conformance_driver: missing "sub <name>" header'
         stop 1
      end if
      line = adjustl(line)
      if (len_trim(line) > 0) exit
   end do
   if (line(1:4) /= 'sub ') then
      write(*,*) 'conformance_driver: first line must be "sub <name>", got: ', trim(line)
      stop 1
   end if
   sub_name = adjustl(line(5:))

   call read_kvs(in_unit, store)

   write(out_unit, '("sub ",A)') trim(sub_name)

   select case (trim(sub_name))
   case ('projgr');  call run_projgr(store, out_unit)
   case ('errclb');  call run_errclb(store, out_unit)
   case ('active');  call run_active(store, out_unit)
   case ('freev');   call run_freev(store, out_unit)
   case ('hpsolb');  call run_hpsolb(store, out_unit)
   case ('cmprlb');  call run_cmprlb(store, out_unit)
   case ('bmv');     call run_bmv(store, out_unit)
   case ('formt');   call run_formt(store, out_unit)
   case ('matupd');  call run_matupd(store, out_unit)
   case ('dcstep');  call run_dcstep(store, out_unit)
   case ('formk');   call run_formk(store, out_unit)
   case ('dcsrch');  call run_dcsrch(store, out_unit)
   case ('lnsrlb');  call run_lnsrlb(store, out_unit)
   case ('cauchy');  call run_cauchy(store, out_unit)
   case ('subsm');   call run_subsm(store, out_unit)
   case default
      write(*,*) 'conformance_driver: unknown subroutine: ', trim(sub_name)
      stop 1
   end select

   call write_end(out_unit)
   close(in_unit)
   close(out_unit)

contains

   ! ===================================================================
   ! projgr
   ! ===================================================================
   subroutine run_projgr(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, nlen
      integer,  allocatable :: nbd(:)
      real(dp), allocatable :: x(:), l(:), u(:), g(:)
      real(dp) :: sbgnrm

      call get_int(s, 'n', n)
      call get_real_vec(s, 'x', x, nlen)
      call get_real_vec(s, 'l', l, nlen)
      call get_real_vec(s, 'u', u, nlen)
      call get_int_vec(s, 'nbd', nbd, nlen)
      call get_real_vec(s, 'g', g, nlen)
      call projgr(n, l, u, nbd, x, g, sbgnrm)
      call write_real(ou, 'sbgnrm', sbgnrm)
   end subroutine run_projgr

   ! ===================================================================
   ! errclb
   ! ===================================================================
   subroutine run_errclb(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, m, info, k, nlen
      integer,  allocatable :: nbd(:)
      real(dp), allocatable :: l(:), u(:)
      real(dp) :: factr
      character(len=60) :: task

      call get_int(s, 'n', n)
      call get_int(s, 'm', m)
      call get_real(s, 'factr', factr)
      call get_real_vec(s, 'l', l, nlen)
      call get_real_vec(s, 'u', u, nlen)
      call get_int_vec(s, 'nbd', nbd, nlen)
      call get_str(s, 'task_in', task)
      call get_int(s, 'info_in', info)
      call get_int(s, 'k_in', k)

      call errclb(n, m, factr, l, u, nbd, task, info, k)

      call write_str(ou, 'task', task)
      call write_int(ou, 'info', info)
      call write_int(ou, 'k', k)
   end subroutine run_errclb

   ! ===================================================================
   ! active
   ! ===================================================================
   subroutine run_active(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, iprint, nlen
      integer,  allocatable :: nbd(:), iwhere(:)
      real(dp), allocatable :: x(:), l(:), u(:)
      logical  :: prjctd, cnstnd, boxed

      call get_int(s, 'n', n)
      call get_real_vec(s, 'l', l, nlen)
      call get_real_vec(s, 'u', u, nlen)
      call get_int_vec(s, 'nbd', nbd, nlen)
      call get_real_vec(s, 'x_in', x, nlen)
      call get_int(s, 'iprint', iprint)
      allocate(iwhere(n))
      iwhere = 0

      call active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed)

      call write_real_vec(ou, 'x', x)
      call write_int_vec(ou, 'iwhere', iwhere)
      call write_bool(ou, 'prjctd', prjctd)
      call write_bool(ou, 'cnstnd', cnstnd)
      call write_bool(ou, 'boxed', boxed)
   end subroutine run_active

   ! ===================================================================
   ! freev
   ! ===================================================================
   subroutine run_freev(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, nfree_in, iter, iprint, nenter, ileave, nfree, nlen
      integer,  allocatable :: index_arr(:), indx2(:), iwhere(:)
      logical  :: updatd, cnstnd, wrk

      call get_int(s, 'n', n)
      call get_int(s, 'nfree_in', nfree_in)
      call get_int_vec(s, 'index_in', index_arr, nlen)
      call get_int_vec(s, 'indx2_in', indx2, nlen)
      call get_int_vec(s, 'iwhere', iwhere, nlen)
      call get_bool(s, 'updatd', updatd)
      call get_bool(s, 'cnstnd', cnstnd)
      call get_int(s, 'iter', iter)
      call get_int(s, 'iprint', iprint)
      nfree = nfree_in

      call freev(n, nfree, index_arr, nenter, ileave, indx2, &
                 iwhere, wrk, updatd, cnstnd, iprint, iter)

      call write_int(ou, 'nfree', nfree)
      call write_int_vec(ou, 'index', index_arr)
      call write_int_vec(ou, 'indx2', indx2)
      call write_int(ou, 'nenter', nenter)
      call write_int(ou, 'ileave', ileave)
      call write_bool(ou, 'wrk', wrk)
   end subroutine run_freev

   ! ===================================================================
   ! hpsolb
   ! ===================================================================
   subroutine run_hpsolb(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, iheap, nlen
      integer,  allocatable :: iorder(:)
      real(dp), allocatable :: t(:)

      call get_int(s, 'n', n)
      call get_real_vec(s, 't_in', t, nlen)
      call get_int_vec(s, 'iorder_in', iorder, nlen)
      call get_int(s, 'iheap', iheap)

      call hpsolb(n, t, iorder, iheap)

      call write_real_vec(ou, 't', t)
      call write_int_vec(ou, 'iorder', iorder)
   end subroutine run_hpsolb

   ! ===================================================================
   ! bmv
   ! ===================================================================
   subroutine run_bmv(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: m, col, info, rows, cols, vlen
      real(dp), allocatable :: sy(:,:), wt(:,:), vv(:), p(:)

      call get_int(s, 'm', m)
      call get_int(s, 'col', col)
      call get_real_mat(s, 'sy', sy, rows, cols)
      call get_real_mat(s, 'wt', wt, rows, cols)
      call get_real_vec(s, 'v', vv, vlen)
      call get_real_vec(s, 'p_in', p, vlen)
      ! bmv only sets info on the col>0 path; init to 0 so the col=0
      ! early-return produces a defined output.
      info = 0

      call bmv(m, sy, wt, col, vv, p, info)

      call write_real_vec(ou, 'p', p)
      call write_int(ou, 'info', info)
   end subroutine run_bmv

   ! ===================================================================
   ! formt
   ! ===================================================================
   subroutine run_formt(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: m, col, info, rows, cols
      real(dp) :: theta
      real(dp), allocatable :: wt(:,:), sy(:,:), ss(:,:)

      call get_int(s, 'm', m)
      call get_int(s, 'col', col)
      call get_real(s, 'theta', theta)
      call get_real_mat(s, 'sy', sy, rows, cols)
      call get_real_mat(s, 'ss', ss, rows, cols)
      call get_real_mat(s, 'wt_in', wt, rows, cols)

      call formt(m, wt, sy, ss, col, theta, info)

      call write_real_mat(ou, 'wt_upper', wt)
      call write_int(ou, 'info', info)
   end subroutine run_formt

   ! ===================================================================
   ! matupd
   ! ===================================================================
   subroutine run_matupd(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, m, iupdat, head, itail, col, rows, cols, nlen
      real(dp) :: theta, rr, dr, stp, dtd
      real(dp), allocatable :: ws(:,:), wy(:,:), sy(:,:), ss(:,:), d(:), r(:)

      call get_int(s, 'n', n);   call get_int(s, 'm', m)
      call get_int(s, 'iupdat', iupdat)
      call get_int(s, 'head_in', head); call get_int(s, 'itail_in', itail)
      call get_int(s, 'col_in', col)
      call get_real_mat(s, 'ws_in', ws, rows, cols)
      call get_real_mat(s, 'wy_in', wy, rows, cols)
      call get_real_mat(s, 'sy_in', sy, rows, cols)
      call get_real_mat(s, 'ss_in', ss, rows, cols)
      call get_real_vec(s, 'd', d, nlen);  call get_real_vec(s, 'r', r, nlen)
      call get_real(s, 'rr', rr); call get_real(s, 'dr', dr)
      call get_real(s, 'stp', stp); call get_real(s, 'dtd', dtd)
      ! theta is overwritten by matupd; provide any starting value.
      theta = 1.0_dp

      call matupd(n, m, ws, wy, sy, ss, d, r, itail, &
                  iupdat, col, head, theta, rr, dr, stp, dtd)

      call write_int(ou, 'head', head)
      call write_int(ou, 'itail', itail)
      call write_int(ou, 'col', col)
      call write_real(ou, 'theta', theta)
      call write_real_mat(ou, 'ws', ws)
      call write_real_mat(ou, 'wy', wy)
      call write_real_mat(ou, 'sy', sy)
      call write_real_mat(ou, 'ss', ss)
   end subroutine run_matupd

   ! ===================================================================
   ! dcstep
   ! ===================================================================
   subroutine run_dcstep(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      real(dp) :: stx, fx, dx, sty, fy, dy, stp, fp, dp_, stpmin, stpmax
      logical  :: brackt

      call get_real(s, 'stx_in', stx); call get_real(s, 'fx_in', fx)
      call get_real(s, 'dx_in', dx)
      call get_real(s, 'sty_in', sty); call get_real(s, 'fy_in', fy)
      call get_real(s, 'dy_in', dy)
      call get_real(s, 'stp_in', stp); call get_real(s, 'fp', fp)
      call get_real(s, 'dp', dp_)
      call get_bool(s, 'brackt_in', brackt)
      call get_real(s, 'stpmin', stpmin); call get_real(s, 'stpmax', stpmax)

      call dcstep(stx, fx, dx, sty, fy, dy, stp, fp, dp_, brackt, stpmin, stpmax)

      call write_real(ou, 'stx', stx); call write_real(ou, 'fx', fx)
      call write_real(ou, 'dx', dx)
      call write_real(ou, 'sty', sty); call write_real(ou, 'fy', fy)
      call write_real(ou, 'dy', dy)
      call write_real(ou, 'stp', stp); call write_bool(ou, 'brackt', brackt)
   end subroutine run_dcstep

   ! ===================================================================
   ! formk
   ! ===================================================================
   subroutine run_formk(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, m, nsub, nenter, ileave, iupdat, col, head, info, &
                  rows, cols, nlen
      logical  :: updatd
      real(dp) :: theta
      integer,  allocatable :: ind(:), indx2(:)
      real(dp), allocatable :: wn(:,:), wn1(:,:), ws(:,:), wy(:,:), sy(:,:)

      call get_int(s, 'n', n);  call get_int(s, 'm', m)
      call get_int(s, 'nsub', nsub)
      call get_int(s, 'nenter', nenter); call get_int(s, 'ileave', ileave)
      call get_int(s, 'iupdat', iupdat)
      call get_bool(s, 'updatd', updatd)
      call get_int(s, 'col', col); call get_int(s, 'head', head)
      call get_real(s, 'theta', theta)
      call get_int_vec(s, 'ind', ind, nlen)
      call get_int_vec(s, 'indx2', indx2, nlen)
      call get_real_mat(s, 'ws', ws, rows, cols)
      call get_real_mat(s, 'wy', wy, rows, cols)
      call get_real_mat(s, 'sy', sy, rows, cols)
      call get_real_mat(s, 'wn_in', wn, rows, cols)
      call get_real_mat(s, 'wn1_in', wn1, rows, cols)

      call formk(n, nsub, ind, nenter, ileave, indx2, iupdat, &
                 updatd, wn, wn1, m, ws, wy, sy, theta, col, head, info)

      call write_int(ou, 'info', info)
      call write_real_mat(ou, 'wn', wn)
      call write_real_mat(ou, 'wn1', wn1)
   end subroutine run_formk

   ! ===================================================================
   ! cmprlb
   ! ===================================================================
   subroutine run_cmprlb(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, m, col, head, nfree, info, rows, cols, nlen
      logical  :: cnstnd
      real(dp) :: theta
      integer,  allocatable :: idx(:)
      real(dp), allocatable :: ws(:,:), wy(:,:), sy(:,:), wt(:,:)
      real(dp), allocatable :: x(:), g(:), z(:), r(:), wa(:)

      call get_int(s, 'n', n); call get_int(s, 'm', m)
      call get_int(s, 'col', col); call get_int(s, 'head', head)
      call get_int(s, 'nfree', nfree)
      call get_real(s, 'theta', theta)
      call get_bool(s, 'cnstnd', cnstnd)
      call get_int_vec(s, 'index', idx, nlen)
      call get_real_vec(s, 'x', x, nlen);  call get_real_vec(s, 'z', z, nlen)
      call get_real_vec(s, 'g', g, nlen)
      call get_real_vec(s, 'r_in', r, nlen)
      call get_real_vec(s, 'wa_in', wa, nlen)
      call get_real_mat(s, 'ws', ws, rows, cols)
      call get_real_mat(s, 'wy', wy, rows, cols)
      call get_real_mat(s, 'sy', sy, rows, cols)
      call get_real_mat(s, 'wt', wt, rows, cols)
      ! cmprlb only sets info inside the bmv-calling path; init to 0
      ! so the unconstrained-shortcut and col=0 paths produce a
      ! defined output (matches the F77 unit test which sets info=0
      ! on entry).
      info = 0

      call cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, idx, &
                  theta, col, head, nfree, cnstnd, info)

      call write_real_vec(ou, 'r', r)
      call write_int(ou, 'info', info)
   end subroutine run_cmprlb

   ! ===================================================================
   ! dcsrch
   ! ===================================================================
   subroutine run_dcsrch(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      real(dp) :: f, g, stp, ftol, gtol, xtol, stpmin, stpmax
      character(len=60) :: task
      integer  :: isave(2)
      real(dp) :: dsave(13)

      call get_real(s, 'f', f); call get_real(s, 'g', g)
      call get_real(s, 'stp', stp)
      call get_real(s, 'ftol', ftol); call get_real(s, 'gtol', gtol)
      call get_real(s, 'xtol', xtol)
      call get_real(s, 'stpmin', stpmin); call get_real(s, 'stpmax', stpmax)
      call get_str(s, 'task_in', task)
      isave = 0; dsave = 0.0_dp

      call dcsrch(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task, isave, dsave)

      call write_str(ou, 'task', task)
      call write_real(ou, 'stp', stp)
      call write_real(ou, 'f', f)
      call write_real(ou, 'g', g)
   end subroutine run_dcsrch

   ! ===================================================================
   ! lnsrlb
   ! ===================================================================
   subroutine run_lnsrlb(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, iter, ifun, iback, nfgv, info, rows, cols, nlen
      logical  :: boxed, cnstnd
      real(dp) :: f, fold, gd, gdold, stp, dnorm, dtd, xstep, stpmx
      integer,  allocatable :: nbd(:)
      real(dp), allocatable :: l(:), u(:), x(:), g(:), d(:), r(:), t(:), z(:)
      character(len=60) :: task, csave
      integer  :: isave(2)
      real(dp) :: dsave(13)

      call get_int(s, 'n', n)
      call get_real_vec(s, 'l', l, nlen);  call get_real_vec(s, 'u', u, nlen)
      call get_int_vec(s, 'nbd', nbd, nlen)
      call get_real_vec(s, 'x', x, nlen)
      call get_real(s, 'f_in', f)
      call get_real_vec(s, 'g', g, nlen)
      call get_real_vec(s, 'd', d, nlen)
      call get_real_vec(s, 'r_in', r, nlen)
      call get_real_vec(s, 't_in', t, nlen)
      call get_real_vec(s, 'z', z, nlen)
      call get_real(s, 'stp_in', stp)
      call get_int(s, 'iter', iter)
      call get_int(s, 'ifun_in', ifun); call get_int(s, 'iback_in', iback)
      call get_int(s, 'nfgv_in', nfgv); call get_int(s, 'info_in', info)
      call get_str(s, 'task_in', task)
      call get_bool(s, 'boxed', boxed); call get_bool(s, 'cnstnd', cnstnd)
      fold = 0.0_dp; gd = 0.0_dp; gdold = 0.0_dp
      dnorm = 0.0_dp; dtd = 0.0_dp; xstep = 0.0_dp; stpmx = 0.0_dp
      csave = ''
      isave = 0; dsave = 0.0_dp

      call lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, stp, &
                  dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, info, &
                  task, boxed, cnstnd, csave, isave, dsave)

      call write_real(ou, 'stpmx', stpmx)
      call write_real(ou, 'stp', stp)
      call write_real(ou, 'dnorm', dnorm)
      call write_int(ou, 'info', info)
      call write_str(ou, 'task', task)
      call write_int(ou, 'ifun', ifun)
      call write_real_vec(ou, 'x', x)
   end subroutine run_lnsrlb

   ! ===================================================================
   ! cauchy
   ! ===================================================================
   subroutine run_cauchy(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, m, col, head, iprint, info, nseg, rows, cols, nlen
      real(dp) :: theta, sbgnrm, epsmch
      integer,  allocatable :: nbd(:), iwhere(:), iorder(:)
      real(dp), allocatable :: x(:), l(:), u(:), g(:), tt(:), d(:), xcp(:)
      real(dp), allocatable :: ws(:,:), wy(:,:), sy(:,:), wt(:,:)
      real(dp), allocatable :: p(:), c(:), wbp(:), v(:)
      integer  :: col_max

      call get_int(s, 'n', n);  call get_int(s, 'm', m)
      call get_int(s, 'col', col); call get_int(s, 'head', head)
      call get_int(s, 'iprint', iprint)
      call get_real(s, 'theta', theta)
      call get_real(s, 'sbgnrm', sbgnrm)
      call get_real(s, 'epsmch', epsmch)
      call get_real_vec(s, 'x', x, nlen)
      call get_real_vec(s, 'l', l, nlen);  call get_real_vec(s, 'u', u, nlen)
      call get_int_vec(s, 'nbd', nbd, nlen)
      call get_real_vec(s, 'g', g, nlen)
      call get_int_vec(s, 'iwhere_in', iwhere, nlen)
      call get_real_mat(s, 'ws', ws, rows, cols)
      call get_real_mat(s, 'wy', wy, rows, cols)
      call get_real_mat(s, 'sy', sy, rows, cols)
      call get_real_mat(s, 'wt', wt, rows, cols)

      allocate(iorder(n), tt(n), d(n), xcp(n))
      iorder = 0; tt = 0.0_dp; d = 0.0_dp; xcp = -42.0_dp
      col_max = max(2, 2*col)
      allocate(p(col_max), c(col_max), wbp(col_max), v(col_max))
      p = 0.0_dp; c = 0.0_dp; wbp = 0.0_dp; v = 0.0_dp
      info = 0; nseg = 0

      call cauchy(n, x, l, u, nbd, g, iorder, iwhere, tt, d, xcp, &
                  m, wy, ws, sy, wt, theta, col, head, p, c, wbp, &
                  v, nseg, iprint, sbgnrm, info, epsmch)

      call write_real_vec(ou, 'xcp', xcp)
      call write_int_vec(ou, 'iwhere', iwhere)
      call write_int(ou, 'info', info)
      call write_int(ou, 'nseg', nseg)
   end subroutine run_cauchy

   ! ===================================================================
   ! subsm
   ! ===================================================================
   subroutine run_subsm(s, ou)
      type(kv_store), intent(in) :: s
      integer,        intent(in) :: ou
      integer  :: n, m, nsub, col, head, iword, info, iprint
      integer  :: rows, cols, nlen
      real(dp) :: theta
      integer,  allocatable :: ind(:), nbd(:)
      real(dp), allocatable :: l(:), u(:), x(:), d(:), xp(:), xx(:), gg(:)
      real(dp), allocatable :: ws(:,:), wy(:,:), wn(:,:), wv(:)

      call get_int(s, 'n', n); call get_int(s, 'm', m)
      call get_int(s, 'nsub', nsub)
      call get_int(s, 'col', col); call get_int(s, 'head', head)
      call get_real(s, 'theta', theta)
      call get_int_vec(s, 'ind', ind, nlen)
      call get_real_vec(s, 'l', l, nlen);  call get_real_vec(s, 'u', u, nlen)
      call get_int_vec(s, 'nbd', nbd, nlen)
      call get_real_vec(s, 'x_in', x, nlen)
      call get_real_vec(s, 'd_in', d, nlen)
      call get_real_vec(s, 'xx', xx, nlen); call get_real_vec(s, 'gg', gg, nlen)
      call get_real_mat(s, 'ws', ws, rows, cols)
      call get_real_mat(s, 'wy', wy, rows, cols)
      call get_real_mat(s, 'wn_in', wn, rows, cols)

      allocate(xp(n)); xp = -42.0_dp
      allocate(wv(2*m)); wv = 0.0_dp
      iprint = -1; info = 0
      ! F77 subsm.f returns immediately without setting iword when
      ! nsub<=0; we initialise iword = 0 so the early-return path
      ! produces iword=0 (matching the Python ref and JSON expected).
      iword = 0

      call subsm(n, m, nsub, ind, l, u, nbd, x, d, xp, ws, wy, &
                 theta, xx, gg, col, head, iword, wv, wn, iprint, info)

      call write_real_vec(ou, 'x', x)
      call write_int(ou, 'iword', iword)
      call write_int(ou, 'info', info)
   end subroutine run_subsm

end program conformance_driver
