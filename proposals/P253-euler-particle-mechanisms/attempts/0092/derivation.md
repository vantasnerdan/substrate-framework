# Compact gauge character and same-carrier action selection

## 1. The dimensional action is exact but does not select a value

Write the Maxwell normalization as

    c_EM^2=1/(epsilon_EM mu_EM),
    S_EM(Q)=Q^2/(4 pi epsilon_EM c_EM).                    (1)

In charge-explicit dimensions,

    [epsilon_EM]=Q^2 T^2/(M L^3),   [c_EM]=L/T,

so `[S_EM]=M L^2/T`, the dimension of action.  This is a useful invariant
scale already present in the declared Euler--Maxwell extension.  Dimensional
correctness is not a law setting either `Q` or a carrier action.

The normalization distinction is exact.  If the gauge potential is renamed
by `A'=a A`, with the same scaling for the scalar potential, then expressing
the old action in the primed field gives

    epsilon_EM'=epsilon_EM/a^2,
    mu_EM'=mu_EM a^2,
    g_tag'=g_tag/a.                                        (2)

Consequently

    1/(epsilon_EM' mu_EM')=c_EM^2,
    (g_tag')^2/epsilon_EM'=g_tag^2/epsilon_EM.             (3)

The API takes `a>0` as a normalization magnitude and records orientation
separately.  Algebraically (2)--(3) only require `a` nonzero.  The physical
sign of `g_tag` is retained: every transported-tag equation below permits
`g_tag<0` as well as `g_tag>0`.

## 2. Route A: the transported real tag is a continuous charge source

The reviewed current is

    rho_q=g_tag chi,          j_q=g_tag chi u,
    partial_t chi+div(chi u)=0.                            (4)

For every real `lambda`, `chi_lambda=lambda chi` obeys the same transport
equation and has

    Q_lambda=g_tag lambda C_chi,
    C_chi=integral chi dx.                                 (5)

On the prescribed-current Maxwell block, fields scale linearly with
`lambda`, while field energy and (1) scale quadratically.  At fixed `g_tag`,
the different nonzero `lambda` are different material-tag leaves.  Fixing one
leaf conserves its charge; it does not remove the continuum of admissible
leaves or select one nonzero value.  The reviewed 0080 theorem makes this
especially concrete: it freezes `C_chi` as input and obtains
`Q_g=g_tag C_chi`.  It proves a charged branch after the leaf is chosen, not a
law choosing the leaf or coupling.

**Route A verdict.**  Existing transported-tag charge quantization and
nonzero charge selection are refuted by the exact one-parameter family (5).
The conservation law and signed current remain established.

## 3. Route B: a compact material phase still has continuous classical momentum

The smallest direct compact-phase extension has a circle-valued field
`theta`, a dimensionless number density `n`, and an action coefficient `S_0`.
One dimensionally complete convention for the first-order matter term is

    S_phase=integral [S_0 n D_t^u theta
                      -g_0 m n(phi-u dot A)-H(n,...)] dx dt.          (6)

Define the charge-per-action coefficient

    kappa_m=g_0 m/S_0.                                    (6a)

Under

    A -> A+grad Lambda,
    phi -> phi-partial_t Lambda,
    theta -> theta-kappa_m Lambda,                         (7)

the two changes in braces cancel exactly.  The phase canonical one-form and
symplectic form are

    Theta_phase=S_0 integral n delta theta dx,
    Omega_phase=S_0 integral delta n wedge delta theta dx. (8)

For a globally normalized compact gauge parameter `alpha in R/(2 pi Z)`, a
well-defined action on `theta` has character

    theta -> theta-m alpha,      m in Z.                     (7a)

Writing `g_0` for the physical charge unit attached to that normalization,
Gauss charge is, up to orientation, `g_0 m n`, while the global number is
`N_phase=integral n dx`.  The integer `m` is genuine representation data, but
the classical cotangent fiber to `S^1` is `R`.  Hence

    J_phase=S_0 N_phase,      Q=g_0 m N_phase
                             =kappa_m J_phase.               (8a)

This total charge remains classically continuous as `N_phase` varies.

This separates four distinct objects that are often conflated:

1. `w=(2 pi)^-1 integral_gamma d theta` is a spatial winding on a loop where
   `theta` is defined;
2. the classical number `N_phase` and cotangent action
   `J_phase=S_0 N_phase` are real;
3. a chosen global compact `U(1)` character has integer weight `m`, and its
   later quantum representation carries the same discrete label; and
4. prequantization of (8) requires an externally normalized action unit and
   an integral cohomology class.

Neither item 1 nor the integer character in item 3 turns item 2 into an
integer.  Choosing the character is added global representation data.  The
coefficient `S_0` in (6) multiplies (8) and every action period; it is an
imported action normalization, while `g_0` is an imported charge unit.  Under
the field normalization `A'=aA`, `phi'=a phi`, one has `g_0'=g_0/a` and
`kappa_m'=kappa_m/a` with `S_0,m` fixed.  Thus a Wess--Zumino or first-order
coefficient can simply re-import the missing action normalization.

**Route B verdict.**  Classical total-charge quantization from phase
periodicity is refuted.  An integral compact character is established only
after choosing the global representation, and total charge still contains the
continuous classical momentum.  Conditional prequantization or a quantum
charge spectrum additionally requires the action unit and quantization rule.

## 4. Route C: spatial U(1) topology distinguishes magnetic from electric charge

Finite-energy fields that extend smoothly over the one-point compactification
live over `S^3`.  Principal `U(1)` bundles are classified by their first Chern
class in `H^2`; here

    H^2(S^3;Z)=0,             pi_3(U(1))=0.                (9)

Every such whole-space bundle is therefore topologically trivial.  By
contrast,

    H^2(S^2_infinity;Z)=Z.                                 (10)

Its Chern integer is the magnetic flux of the spatial connection.  A nonzero
value requires a removed or singular core, a monopole, or a bundle that does
not extend smoothly across all of space.  It is not the electric Gauss charge

    Q=integral_(S^2_infinity) epsilon_EM *F,               (11)

which is the boundary moment map of Maxwell gauge transformations.  Turning
(10) into electric integrality requires extra duality/monopole data or a
second gauge field, plus a physical normalization.

A compact `3+1` BF candidate makes the cost visible:

    S_BF=(k_BF/(2 pi)) integral B wedge dA.                (12)

Large-gauge invariance can restrict `k_BF` when both gauge fields are compact
and the quantum exponentiated action is specified.  The Euler vorticity
two-form, however, is an ordinary real advected exact form on smooth
whole-space flow.  Its Kelvin periods are continuous data, and it is not a
compact two-form connection with integral periods.  Substituting it for `B`
does not produce a quantized BF level.  Promoting it to a compact two-form
gauge field and choosing `k_BF` adds new degrees of freedom and an action
normalization.

**Route C verdict.**  Smooth whole-space `U(1)` topology does not quantize
electric charge.  Magnetic Chern, dual/monopole, and compact BF routes remain
distinct candidates, but every positive realization currently adds the
singular core, compact higher-form field, quantum rule, or normalization that
must itself be derived and reviewed.

## 5. Route D: fixed charge does not lock the carrier action

Let `J` denote the physical oscillator action from the reviewed
Schwinger--Hopf normalization, `J=sum_a |z_a|^2`.  Before assuming any
dynamics tying the mode to charge, the classical map in the two independent
directions `(lambda,J)` is

    Q=g_tag C_chi lambda,
    A_mode=J,
    S_EM=Q^2/(4 pi epsilon_EM c_EM).                        (13)

Its Jacobian is

           [ g_tag C_chi                0 ]
    D =    [ 0                          1 ] .               (14)
           [ g_tag^2 C_chi^2 lambda/(2 pi epsilon_EM c_EM)  0 ]

For nonzero `g_tag C_chi`, `rank D=2`.  The exact fixed-charge tangent
`(delta lambda,delta J)=(0,delta J)` leaves both `Q` and `S_EM` fixed while
changing `A_mode`.  Thus even an independently earned integral charge would
not select the mode action.

This is consistent with the reviewed Euler similarity ledger: bare Euler
action scales as `A/B^4` while topology is unchanged.  The Maxwell extension
introduces scale-breaking constants, but the reviewed charged Cao theorem
still treats circulation, mean radius, core scale, and tag integral as
continuous prescribed carrier data.  Its bordered inverse establishes a
branch after these data are fixed; it contains no additional row setting
`J=S_EM(Q)` or removing the tangent in (14).

The most economical attempted coupling is to gauge the common oscillator
phase itself, which would make its moment map the charge:

    Q=kappa_m J.                                            (15)

This is a real structural relation, but `J` is still continuous classically.
If one additionally imposes `J=S_EM(Q)`, the equation is

    J[1-kappa_m^2 J/(4 pi epsilon_EM c_EM)]=0,              (16)

with the nonzero algebraic root

    J_*=4 pi epsilon_EM c_EM/kappa_m^2.                     (17)

Both signs of `kappa_m` give the same action magnitude.  Equation (16) is not an
Euler--Maxwell equation; it is precisely the missing selection law written as
a constraint.  With positive oscillator and Coulomb energies, the unforced
classical minimum is the other root `J=0`.  Opposite electric charge at fixed
coupling further requires a signed compact character or conjugate matter
sector; it is not erased by taking a positive action.

**Route D verdict.**  Fixed charge, the Maxwell action scale, and minimal
common-phase gauging do not select a nonzero classical carrier action.  The
rank counterdirection (14) and zero-energy minimum are exact.  A new compact
matter/topological interaction capable of imposing a nonzero level without
an imported action coefficient remains an active candidate.

## 6. Consequence and next constructions

The current `U(1)` extension supplies a causal Maxwell speed, a reciprocal
charged carrier, and the exact dimensionful action (1).  It does not make
that action universal or nonzero.  A successful continuation must derive,
on the same carrier,

- a compact character or topological sector fixing signed `Q/g_0` without an
  inserted action unit;
- a joint constraint or dynamical minimum removing the fixed-`Q` action
  tangent in (14) at a nonzero level;
- the physical doublet response and outgoing two-by-two Feshbach matrix; and
- independent Born, reset, exchange/fermionic, and neutrino mechanisms.

The most exposing next compact candidate is a genuine compact two-form or
material defect sector coupled to `A`.  It must first construct its compact
connection from carrier observables and show that its level and symplectic
coefficient are not free imports.  In parallel, the primary classical P4
route remains the 0088 two-control response followed by the source-specific
matrix-valued shell measure.  No route verdict here narrows either parent
objective.
