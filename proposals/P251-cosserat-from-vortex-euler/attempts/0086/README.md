# 0086 — material-map averaging and the reconstruction equation

Owner `/root`; continuation of the active P251 same-action construction.
Frozen representation change: construct the material map first and use its
actual velocity, rather than transplanting a vorticity-orbit Hessian to a
material displacement. Compare the maps by exact differentiation before
any solver, empirical comparator, or coefficient selection. The direct
material route is being implemented in 0084; this record supplies the
averaging license and an exposing map distinction.

## Primary-source reconciliation

Holm, *Variational Principles for Lagrangian Averaged Fluid Dynamics*,
[arXiv:nlin/0103043v1](https://arxiv.org/html/nlin/0103043), sections 2–4
were read. Theorem 3.3 treats a factorized material map and an averaging
operator with the stated projection property. Equations (3.2), (3.5) and
(3.8) transport velocity, scalars and density through that SAME map. Its
mean equation drops variations of the fluctuation map; section 4.2 explicitly
assigns self-consistent fluctuation dynamics to those omitted variations.
Section 4.3 identifies statistical closure as additional input. Therefore
this theorem supports a material-first averaging representation, not an
automatic independent microrotation equation or a positive EPS locking
coefficient. A retained angle requires its fluctuation-map variation.
This is a declared method import, not an accepted framework claim.

## Exact map distinction, derived here

Let g0(t) be the material flow of stationary u0 and take an actual
volume-preserving displacement map Phi(t). Define g=Phi(t) o g0(t).
The chain rule, with no dynamical approximation, gives

    u = Phi_t o Phi^-1 + Phi_* u0.

At first order Phi=id+epsilon xi the Eulerian velocity variation is

    v = xi_t + (u0.grad)xi - (xi.grad)u0
      = xi_t + curl(xi cross u0).

Every passive material tag then has variation -xi.grad chi0 and obeys
the matching linear transport equation. This is the material repair used
in 0084. It does not say that every prescribed Phi solves Euler.

For a fixed-Kelvin variation with the same material vorticity datum,
the vorticity variation must instead ALSO satisfy

    curl v = curl(xi cross omega0).

On the full-space zero-mean/decaying velocity sector, taking the actual
Leray/Biot–Savart inverse identifies v=P(xi cross omega0). For a Beltrami
background omega0=lambda u0 these two requirements combine to

    xi_t = (lambda P-curl)(xi cross u0).

The pressure/zero-mode or circulation-period conditions of another domain
must be retained. The equation is a reconstruction identity, not an
additional demand for an all-wave-number invariant nonlinear ansatz.
A restricted material action evaluates its own Jacobi form; an orbit
restriction evaluates its own KKS/Hessian form. Equating them requires the
displayed compatibility or an explicit constrained reduction with its
reaction terms. Geometric inertia alone does not supply that equation.

An elementary exact counterexample exposes the distinction without any
numerical approximation. Let u0=(cos z,-sin z,0), so curl u0=u0, and
F=diag(a,1,a^-1), a>0. The static material pushforward is

    v_material=(a cos(a z),-sin(a z),0),
    curl v_material=(a cos(a z),-a^2 sin(a z),0).

The pushforward vorticity is instead

    omega_push=(a cos(a z),-sin(a z),0),
    v_BS=(a^-1 cos(a z),-sin(a z),0).

Both velocities are divergence free, and curl v_BS=omega_push, but they
are different whenever a!=1. Their averaged energies are respectively
rho(a^2+1)/4 and rho(a^-2+1)/4. Thus even a simple volume-preserving map
does not license copying the orbit's energy coefficient into the direct
material kinetic action. The full same-action constructions in 0082/0084
remain the active continuation. No parent no-go or completion follows.
