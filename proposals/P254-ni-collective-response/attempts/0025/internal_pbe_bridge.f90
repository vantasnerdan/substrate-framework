! Narrow QE 7.5 internal PBE adapter, ordinary NC / no NLCC / collinear only.
! Algebra follows PW/src/v_of_rho.f90, gradcorr.f90 and LR_Modules/setup_*.f90.
! All functional evaluations and finite differences remain in original xc_lib.a.
module internal_pbe_bridge
  use iso_c_binding, only: c_int, c_double
  use xc_lib, only: xc, xc_gcx, dmxc, dgcxc, xclib_reset_dft, &
                    xclib_set_dft_from_name, xclib_get_ID
  implicit none
  private
  public :: ni_pbe_abi, ni_pbe_evaluate
contains
  integer(c_int) function ni_pbe_abi() bind(C)
    ni_pbe_abi = 1_c_int
  end function

  subroutine ni_pbe_evaluate(nr, rho, grho, local, flux, mu, rr, sr, ss, s, work, status) bind(C)
    integer(c_int), value, intent(in) :: nr
    real(c_double), intent(in) :: rho(nr,2), grho(3,nr,2)
    real(c_double), intent(out) :: local(nr,2), flux(3,nr,2), mu(nr,2,2)
    real(c_double), intent(out) :: rr(nr,2,2), sr(nr,2,2), ss(nr,2,2), s(nr,2,2)
    real(c_double), intent(inout) :: work(nr,19)
    integer(c_int), intent(out) :: status
    status = 1_c_int
    if (nr < 1) return
    ! Private library state; Python serializes calls. Reset also restores the
    ! H0 defaults after dgcxc changed the GGA threshold in the previous call.
    call xclib_reset_dft()
    call xclib_set_dft_from_name('PBE')
    status = 2_c_int
    if (xclib_get_ID('LDA','EXCH') /= 1 .or. xclib_get_ID('LDA','CORR') /= 4) return
    if (xclib_get_ID('GGA','EXCH') /= 3 .or. xclib_get_ID('GGA','CORR') /= 4) return
    call evaluate_fields(nr, rho, grho, local, flux, mu, rr, sr, ss, s, &
         work(1,1), work(1,3), work(1,9), work(1,10), work(1,11), &
         work(1,13), work(1,15), work(1,17), work(1,19))
    status = 0_c_int
  end subroutine

  subroutine evaluate_fields(nr, rho, grho, local, flux, mu, rr, sr, ss, s, &
                             updw, grh, ex, ec, vx, vc, v2x, v2c, v2ud)
    integer, intent(in) :: nr
    real(c_double), intent(in) :: rho(nr,2), grho(3,nr,2)
    real(c_double), intent(out) :: local(nr,2), flux(3,nr,2), mu(nr,2,2)
    real(c_double), intent(out) :: rr(nr,2,2), sr(nr,2,2), ss(nr,2,2), s(nr,2,2)
    real(c_double), intent(out) :: updw(nr,2), grh(nr,3,2), ex(nr), ec(nr)
    real(c_double), intent(out) :: vx(nr,2), vc(nr,2), v2x(nr,2), v2c(nr,2), v2ud(nr)
    real(c_double), parameter :: e2 = 2.0_c_double, epsr = 1.0e-6_c_double
    integer :: ir, is, ipol
    ! xc requires TOTAL/MAGNETIZATION, unlike all the remaining calls.
    call xc(nr, 2, 2, rho, ex, ec, vx, vc)
    do is = 1, 2
      do ir = 1, nr
        local(ir,is) = e2 * (vx(ir,is) + vc(ir,is))
        updw(ir,is) = (rho(ir,1) + real(3-2*is,c_double)*rho(ir,2))*0.5_c_double
      end do
    end do
    ! H0 GGA: default density cutoff 1e-6 and gradient cutoff 1e-10.
    call xc_gcx(nr, 2, updw, grho, ex, ec, vx, v2x, vc, v2c, v2ud)
    do is = 1, 2
      do ir = 1, nr
        local(ir,is) = local(ir,is) + e2*(vx(ir,is) + vc(ir,is))
        do ipol = 1, 3
          flux(ipol,ir,is) = e2*((v2x(ir,is)+v2c(ir,is))*grho(ipol,ir,is) &
                                + v2ud(ir)*grho(ipol,ir,3-is))
          grh(ir,ipol,is) = grho(ipol,ir,is)
        end do
      end do
    end do
    ! Original finite-difference drivers return Ry, not Hartree. No symmetrization.
    call dmxc(nr, 2, updw, mu)
    call dgcxc(nr, 2, updw, grh, rr, sr, ss)
    ! dgcxc sets density/gradient cutoffs to 1e-10. Native setup_dgc calls
    ! xc_gcx AFTER that mutation: reusing the H0 values would be incorrect.
    call xc_gcx(nr, 2, updw, grho, ex, ec, vx, v2x, vc, v2c, v2ud)
    s = 0.0_c_double
    do ir = 1, nr
      if (updw(ir,1)+updw(ir,2) > epsr) then
        s(ir,1,1) = e2*(v2x(ir,1)+v2c(ir,1))
        s(ir,1,2) = e2*v2c(ir,1)
        s(ir,2,1) = e2*v2c(ir,1)
        s(ir,2,2) = e2*(v2x(ir,2)+v2c(ir,1))
      end if
    end do
  end subroutine
end module
