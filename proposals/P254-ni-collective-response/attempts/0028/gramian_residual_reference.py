"""Independent two-level residual integral, avoiding oscillatory quadrature loss."""
import mpmath as mp


def two_level_residual_integral(coupling, gap, time):
    """Integral |[exp(-it [[0,d],[d,w]])]_{21}| for real d,w, t>=0."""
    if coupling < 0 or time < 0:
        raise ValueError('nonnegative coupling/time required')
    with mp.workdps(80):
        d,w,t=mp.mpf(coupling),mp.mpf(gap),mp.mpf(time)
        frequency=mp.sqrt(w*w/4+d*d)
        if frequency == 0:
            return 0.
        phase=frequency*t
        periods=mp.floor(phase/mp.pi)
        remainder=phase-periods*mp.pi
        return float(d/frequency*(2*periods+1-mp.cos(remainder))/frequency)
