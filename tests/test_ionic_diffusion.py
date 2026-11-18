from substrate_framework.ionic_diffusion import diffusive_ionic_shift_eV


def test_zero_mobile_charge_has_no_static_or_moving_ionic_shift():
    # A dry/empty target must not become0/0 at a finite projectile speed.
    for speed in (0.0, 1e6):
        assert diffusive_ionic_shift_eV(0.0, 1e-10, speed) == 0.0
