# P253/0090: comoving Maxwell shell-to-flux and finite-gate leakage

## Ownership and frozen boundary

This is a root-owned analytic continuation of the reviewed Euler--Maxwell
radiation-shell boundary.  Before any derivation body is written it freezes
the exact conversion from a conserved localized comoving current to outgoing
Maxwell power, and then to a finite-gate leakage or one-block Feshbach
self-energy.  It does not assume a nonzero source-specific 0079/0088 response.

Accepted inputs are limited to:

- P253/0068 as independently reviewed by P253/0071: the conserved material
  current, joint Maxwell action, Poynting work identity, typed constraint
  kernel, and Maxwell-subsystem finite propagation;
- corrected P253/0085 as independently reviewed by P253/0086, review SHA-256
  `4c5c535e60a134ba589ac066aa906249952d47d843c96e7f6a3764686f85ea67`
  and verdict SHA-256
  `afd7e563d782d56f8d9f383a8611caf41cf4f2f4bbeb1410dc089262348f9c3d`:
  the full constrained-generator inclusion `i R subset sigma_ess`, the
  comoving shell, positive radial-root/coarea rows, and dark-current
  necessity.

P253/0087 is under active independent review0089 and is not consumed as an
established supplier.  Its Bessel estimate is a downstream consumer only.
P253/0088 is active and excluded.  No electron, neutrino, P4/P5, resonance,
or gate result is imported.

## Frozen source and Fourier convention

Use real physical fields represented by phasors with time average
`<fg>_t=(1/2) Re(f conjugate(g))`, and the unitary spatial Fourier transform

    f_hat(k)=(2 pi)^(-3/2) integral exp(-i k dot x) f(x) dx.       (1)

Let a smooth compact or sufficiently weighted conserved current have the
comoving form

    (rho,J)(t,x)
      =Re exp(-i omega t) (rho_omega,J_omega)(x-c_g t e_z),       (2)

with `omega!=0`, `|c_g|<c_EM`,
the same-convention continuity row

    (omega+c_g k_z)rho_hat=k dot J_hat,                         (2a)

and enough first moment and trace regularity for restriction to the radiation
shell.  At `k=0`, (2a) and `omega!=0` give `rho_hat_omega(0)=0`: the
oscillatory mode is net neutral.  A nonzero charge of the stationary carrier
belongs to the separate `omega=0` background and is not part of this radiation
shell.  Define

    Omega(k)=omega+c_g k_z,
    D_omega(k)=c_EM^2|k|^2-Omega(k)^2,                            (3)

and use the outgoing `exp(-i Omega t)` prescription
`(D_omega-i0 sign Omega)^(-1)`, with the corresponding `c_EM^2` numerator
when converting from the normalized wave operator.  This boundary value is
selected by retarded evolution/Sommerfeld radiation; the algebraic jump alone
does not distinguish it from a time-symmetric field.  On each shell component
`sign Omega` equals `sign omega`.  The longitudinal charge/Gauss/Coulomb row
is solved using (2a) and retained in the exact field, but it is separated from
the transverse radiative numerator `P_T J_hat`.

## Route A: frequency-domain outgoing identity

Construct the weighted outgoing Maxwell resolvent on the conserved-current
subspace, including the gauge quotient, zero mode, and the translating
frequency `Omega(k)`.  Derive rather than postulate the boundary-value
identity

    A_T_hat=mu_EM c_EM^2
      (D_omega-i0 sign omega)^(-1) P_T J_hat,
    E_T_hat=i Omega A_T_hat,                                   (4)

with the chosen signs.  Combine Poynting's theorem, Parseval, and
`Im(D-i0)^(-1)=pi delta(D)` to derive the exact nonnegative mean outward
power.  Under the `exp(-i Omega t)` convention the sign must be fixed directly
from

    P_out=-<integral J dot E dx>_t,                             (4a)

so that inserting (4), including its complex conjugate in the time average,
produces a positive coefficient.  The target form is

    P_out=C_EM integral |Omega(k)| delta(D_omega(k))
                     |P_T(k)J_hat(k)|^2 dk,                    (5)

where `C_EM` and every factor of `2`, `pi`, `epsilon_EM`, and `mu_EM` must be
derived under (1)--(2).  Derive the full shell gradient

    |grad_k D_omega|
      =2|c_EM^2 k-c_g Omega(k)e_z|,                            (5a)

prove it is uniformly nonzero for `|c_g|<c_EM`, and reduce (5) to the exact
star-shaped sphere integral using

    r(n)=|omega|/[c_EM-sign(omega)c_g n_z],
    |partial_r D|=2|omega|c_EM.                                (6)

Prove the weighted trace estimate needed to make (5) finite and continuous;
do not infer it from compact support alone without the declared Sobolev trace
order.  Treat both `+omega` and `-omega` components of the real current and
show explicitly whether the phasor `1/2` already accounts for their common
power or whether a two-shell sum remains.  Exact shell cancellation gives
zero power at this order but remains only a BIC candidate until the coupled
homogeneous problem is solved.

**Route-A verdict target:** established for the prescribed conserved-current
Maxwell problem if the outgoing resolvent, constants, trace domain, and
Poynting limit close; otherwise blocked at the named missing analytic map.

## Route B: time-domain finite-window flux

Independently use retarded Maxwell evolution with a declared real smooth
envelope.  Multiplying both mode densities by `a_T(t)` is not conserved: its
continuity defect is `a_T'(t) rho_omega`.  Because the oscillatory charge has
zero mean, construct on the declared compact/weighted domain a Bogovskii or
Hodge completion `K_omega` satisfying

    div K_omega=-rho_omega,
    i k dot K_hat_omega=-rho_hat_omega,                         (7a)

and use the exactly conserved complex source

    rho_T=a_T(t) exp(-i omega t) rho_omega(x-c_g t e_z),
    J_T=exp(-i omega t)
          [a_T(t)J_omega(x-c_g t e_z)
           +a_T'(t)K_omega(x-c_g t e_z)].                      (7b)

Take real parts only after this identity is verified.  Derive the transverse
and longitudinal pieces of the `a_T' K_omega` source and include their energy
in `E_switch`.  A sharp truncation is a distributional/weak-source limit, not
automatically a member of the smooth-source theorem: type its spatial
completion in the energy space, compute its endpoint field energy, and compare
it with specified smooth regularizations.  Do not assume that this endpoint
energy is infinite.  Compare its sinc tails and endpoint energy with a
smooth/adiabatic switch.  Derive a finite-time energy identity with initial
radiation specified, separate switching transients from the central periodic
interval, and decide whether the strongest supported result is an equality
`E_rad=T P_out+O(1)`, a one-sided bound, or only a spectral density.  The
candidate upper bound is

    E_rad([0,T])
      <=T P_out |a|^2+E_switch+E_near(T),                       (7)

with explicit norms and uniform constants.  Determine when `E_switch` is
subleading, when interference between nearby frequencies survives, and how
the estimate changes for a real `+/-omega` pair.  Maxwell finite propagation
may control the radiative field, while incompressible Euler pressure still
has no strict cone.

Define `E_rad` by one of two operationally closed observables and derive its
balance from the same Poynting identity: either (a) the late-time outgoing
transverse free-field energy after the conserved envelope has switched off,
with the static Coulomb/background energy subtracted, or (b) flux through a
comoving world tube containing the translating source.  In case (b) the
boundary flux is the relative flux

    (S-c_g u_EM e_z) dot n,                                    (7c)

and the balance retains the endpoint field energy stored inside the tube.
A fixed bounded laboratory sphere is not an admissible gate observable for a
source that translates a distance of order `T`.  The chosen observable fixes
the meanings of `E_near` and `E_switch` and must be reconciled with the
frequency-domain source-work identity (4a).

**Route-B verdict target:** an actual envelope-specific finite-window leakage
theorem, not a stationary phasor slogan.  If switching or near-field energy
cannot be controlled at the required gate scale, retain (5) as spectral
density and name the missing time estimate.

## Route C: one-block Feshbach and source-specific consumption

For a normalized positive-Krein mode with physical action `A_mode` and
positive physical frequency `nu_phys`, couple
one amplitude `z` reciprocally to the transverse Maxwell continuum through
the actual current map `J=z j+conjugate(z j)`.  Derive both the amplitude and
field equations from the same joint action, then prove that the one-block
self-energy boundary imaginary part agrees with (5).  A shell delta does not
become damping until this reciprocal amplitude equation and a nonzero
projected current trace exist.  A resonance pole or damping width is earned
only after analytic continuation/Fredholm closure and the fluid/tag feedback
are included.

For a finite gate, the downstream criterion is

    [T_gate P_out]/[nu_phys A_mode] <= tolerance.               (8)

Here `nu_phys A_mode` is the mode energy for the normalization used in the
current profile.  Any equivalent amplitude convention must display its energy
coefficient explicitly; power divided by action alone is not a dimensionless
leakage fraction.

P253/0087 may later supply a shell-trace upper bound, while P253/0088 may
later supply physical `A_mode`, response, and gate histories.  Their inputs
remain typed separately.  An exponentially small nonzero trace produces no
exact `L2` eigenmode, and exponential asymptotics do not prove a useful value
inside a finite upper ceiling without an explicit constant-dependent lower
threshold.

**Route-C verdict target:** exact shell functional to self-energy/finite-gate
conversion at the prescribed-current level, with resonance existence and
source-specific smallness separately verdicted.

## Evidence and completion boundary

The strongest oracle will independently derive the distributional sign and
coarea constants from (1)--(6), then test a smooth transverse Gaussian source
whose shell integral can be evaluated by direct angular quadrature.  This is
an analytic normalization oracle, not production numerics or a substitute for
the weighted resolvent.  No small Hessian eigenvalue, stability edge, force,
or fitted energy splitting is evaluated, so the small-ratio numerical protocol
does not bind this attempt.

The attempt must give each route one verdict and continue from any failed
route by the other representation.  It cannot establish an autonomous Euler
response, action normalization, charged-carrier stability, a detector/reset,
an electron, a neutrino, P4/P5, or the parent campaign.  Its positive purpose
is the exact radiation conversion needed by those later constructions.
