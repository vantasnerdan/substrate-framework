# Attempt 0004 — N2 block B2 (twist stiffness), axisymmetric-reduction route

## Route

Declared ansatz: vorticity lines wind about the tube axis at twist rate chi'
with the z-flux density omega_z(r) = Gamma/(pi a^2) * (1 - chi'^2 r^2/2) /
(1 - chi'^2 a^2/4) (top-hat base mollifier, exact flux constraint), axisymmetric
reduction of the twisted state; enclosed-flux fraction g_chi(r) computed by
exact integration; Delta E_twist/L = (rho/2) int (u_chi^2 - u_0^2) 2 pi r dr.

## Command

PYTHONPATH=src .venv/bin/python (interactive SymPy session, kernel log;
integral: int_0^a (Gamma^2/(4 pi^2)) (g_chi^2 - g_0^2)/r dr)

## Result

g_chi(r) = r^2/a^2 - chi'^2 r^4/(2 a^2) + O(chi'^4)  (inside core),
Delta E_twist/L = - rho Gamma^2 chi'^2 a^2 / (96 pi) + O(chi'^4 / a^0).

## Verdict

REFUTED as a twist-modulus candidate: the energy change is NEGATIVE. The
axisymmetric reduction of the twisted state does not describe the twist
degree of freedom; it describes an inward redistribution of z-flux at fixed
total circulation, which genuinely lowers the log-energy. The physical twist
variable is the theta-DEPENDENT polarization of the core lines; its energy
enters through the non-axisymmetric perturbation sector (coaxial-helix pair
kernel / Kelvin spinning mode), not through the axisymmetric flux profile.

## Failure scope and next route

Scope: the axisymmetric-reduction ansatz only (a B2 sub-route). The N2
obligation stays active. Next materially different attempt (0005): the
coaxial-helix bundle kernel -- pairwise Biot-Savart interaction of wound
line shells dGamma(r) = Gamma f(r) 2 pi r dr with wavenumber k = chi',
regularized by the declared core mollifier, in the declared long-wave
regime k a << 1; extract the O(k^2) coefficient as the (log-running)
twist modulus C_tw(k-window), with the window cutoff declared, mirroring
the bend route. Expected structure: both moduli log-running in the
declared window; constants exact and profile-pinned.
