! Exact serial/global, non-gamma ggen enumeration from QE 7.5
! Modules/recvec_subs.f90:112-200. Link ORIGINAL Modules/kind.f90 and
! Modules/sort.f90; hpsort_eps is deliberately not reimplemented here.
! stdin: nr1 nr2 nr3; tpiba ecutrho_Ry; three bg COLUMN records.
! stdout: ngm; then ngm ordered Miller triples, one per record.
program qe_density_order
  use kinds, only: DP
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  integer :: nr(3), bounds(3), capacity, ng, i, j, k, c, ios
  integer, allocatable :: mill(:,:), order(:)
  real(DP) :: tpiba, ecutrho, bg(3,3), gcutm, tx(3), ty(3), t(3), norm2
  real(DP), allocatable :: squared(:), tt(:)
  real(DP), parameter :: eps8 = 1.0d-8

  read(*, *, iostat=ios) nr
  if (ios /= 0) call fail('invalid FFT dimensions')
  if (any(nr < 2)) call fail('FFT dimensions must exceed one')
  read(*, *, iostat=ios) tpiba, ecutrho
  if (ios /= 0) call fail('invalid tpiba/ecutrho')
  if (.not. (tpiba > 0.0_DP .and. ecutrho > 0.0_DP)) call fail('nonpositive cutoff or tpiba')
  do c = 1, 3
    read(*, *, iostat=ios) bg(:,c)
    if (ios /= 0) call fail('invalid reciprocal column')
  end do
  bounds = (nr - 1) / 2
  capacity = product(2 * bounds + 1)
  allocate(mill(3,capacity), order(capacity), squared(capacity), tt(nr(3)))
  gcutm = ecutrho / tpiba**2
  ng = 0
  do i = -bounds(1), bounds(1)
    tx(1:3) = i * bg(1:3,1)
    do j = -bounds(2), bounds(2)
      ty(1:3) = tx(1:3) + j * bg(1:3,2)
      do k = -bounds(3), bounds(3)
        t(1) = ty(1) + k * bg(1,3)
        t(2) = ty(2) + k * bg(2,3)
        t(3) = ty(3) + k * bg(3,3)
        tt(k+bounds(3)+1) = t(1)**2 + t(2)**2 + t(3)**2
      end do
      do k = -bounds(3), bounds(3)
        norm2 = tt(k+bounds(3)+1)
        if (norm2 <= gcutm) then
          ng = ng + 1
          if (norm2 > eps8) then
            squared(ng) = norm2
          else
            squared(ng) = 0.0_DP
          end if
          mill(:,ng) = [i,j,k]
        end if
      end do
    end do
  end do
  if (ng < 1) call fail('empty reciprocal sphere')
  order(1) = 0
  call hpsort_eps(ng, squared, order, eps8)
  write(*, '(i0)') ng
  do c = 1, ng
    write(*, '(3(i0,1x))') mill(:,order(c))
  end do
contains
  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(error_unit, '(a)') message
    stop 1
  end subroutine fail
end program qe_density_order
