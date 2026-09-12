# Nonlinear action locking and compact-defect charge

## 1. Typed variables and the reduced-energy question

Let `J>=0` be the physical canonical action of one selected positive-Krein
carrier mode.  At fixed physical charge `Q` and fixed carrier parameters
`lambda_carrier`, a legitimate same-carrier reduction would have

    E_red(J)=E_0+nu_eff J+a_2 J^2+R(J).                    (1)

Here `[nu_eff]=T^-1`, `[a_2]=energy/action^2`, and `R=O(J^3)`.
The coefficient `nu_eff` is the derivative of the *real* reduced Hamiltonian.
The reviewed outgoing Maxwell theorem determines a nonnegative shell power,
which is an imaginary boundary value after a reciprocal Feshbach map is
constructed.  It does not determine the principal-value real shift in (1).

The reviewed positive-Krein mode fixes the sign of its uncharged constrained
quadratic energy after choosing the KKS orientation.  It does not supply the
absolute physical normalization, the third and fourth constrained
Euler--Maxwell variations, or an invariant finite-dimensional manifold.  The
strong saddle and Maxwell essential continuum make each of those distinctions
load bearing.

## 2. Exact local sign theorem

Write

    F(J,nu)=partial_J E_red
           =nu+2 a_2 J+r(J),       r=R'.                  (2)

Assume on `0<=J<=J_0` that `a_2>0`, `R` is `C^2`,

    |r(J)|<=C J^2,       |r'(J)|<=C_1 J.                  (3)

At `(J,nu)=(0,0)`, `F=0` and `partial_J F=2a_2>0`.  The
one-dimensional implicit-function theorem therefore gives a unique branch

    J(nu)=-nu/(2a_2)+O(nu^2).                              (4)

For sufficiently small `nu<0`, this root lies in `J>0`, and

    partial_J^2 E_red(J(nu))=2a_2+O(|nu|)>0.              (5)

It is a strict local interior minimum.  For `nu>0`, (3) gives
`partial_J E_red>0` on a sufficiently small interval beginning at zero, so
`J=0` is the local boundary minimum and no arbitrarily small positive minimum
exists.  At `nu=0`, `a_2>0` again makes zero the local minimum.

Equations (2)--(5) are an exact implication, not an existence theorem for the
reduction.  A nonzero action requires a mechanism that drives the real
coefficient through `nu_eff<0`, a derived positive `a_2`, the remainder bounds,
and a reduction radius containing (4).  The current positive oscillator and
Coulomb blocks supply no such negative linear or binding term.  Integer
resonance labels can enumerate branches, but cannot set their amplitude.

## 3. What an actual Euler--Maxwell coefficient calculation requires

Let `e` be a KKS-normalized eigenmode and `q(J,phi)` the corresponding real
phase family.  On a true constrained invariant or steady manifold, `a_2`
contains the fourth variation and the slaving of the complement:

    a_2 = direct fourth-order term
          - (1/2)<quadratic source,
                    (complement Hessian)^(-1) quadratic source>.    (6)

The formula is schematic until all objects share one domain and pairing.  In
this campaign the uncharged fluid Hessian has infinite negative index and the
joint generator has Maxwell essential spectrum `i R`.  Consequently (6)
cannot use an ordinary full-space inverse.  It needs either

1. a source-specific center-stable/invariant manifold with a weighted
   outgoing complement and modulation rows; or
2. a different carrier whose complete constrained DA Hessian is coercive
   modulo symmetries before the gauge field is attached.

The real Maxwell term additionally requires the renormalized principal-value
part of the same outgoing resolvent whose shell jump was reviewed in 0091.
No dispersion relation reconstructs that real part from the shell power
without a subtraction prescription and the actual current.  Thus neither
`nu_eff` nor `a_2` is currently available for the Cao carrier.  Route A earns
the exact sign theorem (2)--(5) and identifies the coefficient calculation;
the same-carrier nonzero minimum remains blocked at (6) and the real
principal-value Feshbach map.

## 4. Compact phase character: exact relation and free action scale

Use the dimensionally complete compact-phase convention from corrected 0092:

    S_phase=integral [S_0 n D_t theta
                      -g_0 m n(phi-u dot A)-H] dx dt,
    m in Z.                                                (7)

For `N=integral n dx`,

    J_phase=S_0 N,
    Q=g_0 m N=kappa_m J_phase,
    kappa_m=g_0 m/S_0.                                    (8)

The integer `m` is genuine chosen representation data.  The classical number
`N` is real, and `S_0` is an action coefficient.  At fixed nonzero `Q,m,g_0`,

    J_phase=Q S_0/(g_0 m),                                 (9)

so changing the admissible coefficient `S_0` changes the locked action
continuously.

The Maxwell constants supply the dimensionally and field-normalization
invariant family

    S_0=C_0 g_0^2/(4 pi epsilon_EM c_EM),                  (10)

where `C_0` is any dimensionless constant.  Gauge-field normalization does not
fix `C_0`: under `A'=aA`, `g_0'=g_0/a` and
`epsilon_EM'=epsilon_EM/a^2`, leaving (10) invariant for every `C_0`.
Compactness fixes `m`, not `S_0` or `C_0`.  Route B is therefore refuted as a
classical action selector unless a same-carrier term independently derives
`S_0`.

## 5. Complement topology and the two-label BF candidate

The relevant integral cohomology rows are elementary:

    R^3 minus {point}  retracts to S^2:
        H^1=0, H^2=Z;                                     (11)

    R^3 minus an embedded S^1:
        H^1=Z, H^2=0.                                     (12)

For (12), Alexander duality on `S^3` gives
`H_tilde_i(S^3-S^1)=H_tilde^(2-i)(S^1)`.  Equation (11) may be read as
duality for the two-point compact set consisting of the physical point and
infinity.  Thus a removed point has the Gauss-surface `H^2` row, whereas a
ring has a linking `H^1` row.  The latter does not itself supply an electric
Gauss character.

Now take two smooth circle-valued material labels away from their defects and
form

    B_label=d alpha wedge d beta.                          (13)

Since `d^2=0`,

    d B_label=d(d alpha wedge d beta)=0                    (14)

on every smooth chart.  Variation of

    S_BF=(k_BF/(2 pi)) integral B wedge dA                 (15)

with respect to `A` couples electric current to `dB`.  The smooth pullback
(13) therefore supplies no bulk electric BF source.  Distributional label
defects may create a source, but then their domain, compact periods, core
energy, and boundary terms are new data.

Ordinary Euler vorticity is a real exact advected two-form.  Its Kelvin
periods vary continuously and do not define a compact gerbe connection.
Replacing `B` in (15) by vorticity imports neither large two-form gauge
transformations nor an integral level.  A genuine compact two-form field can
make `k_BF` integral only after specifying the compact gauge groups and the
exponentiated quantum action.  Its dimensional action coefficient and
same-carrier attachment remain additional structure.  Route C refutes the
smooth Euler-vorticity and smooth two-label suppliers, while the singular
compact-defect route remains blocked at those explicit constructions.

## 6. Similarity and carrier-parameter counterfamilies

For bare Euler similarity `u_(A,B)(x)=A u(Bx)`, the physical KKS/action scale
changes by

    J -> (A/B^4) J.                                       (16)

Positive dilation preserves topology.  In the charged extension, fixing `Q`
removes neither the continuous imported `S_0` direction in (9) nor the
reviewed circulation/radius/core carrier inputs unless an additional joint
constraint is proved.  Equations (9), (10), and (16) are explicit continuous
counterfamilies to universal nonzero action selection by integer character or
topology alone.

## 7. Failure-derived global carrier-scale balance

The preceding obstruction is local in modal action. A materially different
same-carrier mechanism can act on the ring radius itself. Fixed circulation
alone does not define a same-Casimir path, so the coefficient must not be
inferred from `dE=c dI` along the Cao parameter family. Instead start from the
filament Biot--Savart energy

    E_f=(rho_m kappa^2/(8 pi))
        contour_integral contour_integral
        [t(s) dot t(s')]/|x(s)-x(s')| ds ds'.              (17)

For a circular filament, the near-diagonal integral is
`4 pi R log(R/a)` at leading order. Thus, at fixed aspect/profile and to
leading logarithmic order,

    E_fluid=(rho_m kappa^2/2) R L+lower order.             (18)

For charge `Q` distributed around a thin circular core, the Newton double
integral has line-charge leading term

    E_C=Q^2 L/(8 pi^2 epsilon_EM R)+lower order.            (19)

Indeed `lambda_Q=Q/(2 pi R)` and
`|x(theta)-x(theta')|=2R|sin((theta-theta')/2)|`; cutting the
angular integral off at core separation `a/R` gives the coefficient in (19).
The additive finite part depends on the transverse charge profile.

The leading fixed-aspect-ratio energy is consequently

    E_scale(R)=L[A R+B/R],
    A=rho_m kappa^2/2,
    B=Q^2/(8 pi^2 epsilon_EM).                             (20)

For `R>0` and nonzero `kappa,Q`, it has the unique strict minimum

    R_*=|Q|/(2 pi |kappa| sqrt(rho_m epsilon_EM)),
    I_z,*=Q^2/(4 pi epsilon_EM kappa),
    E_scale''(R_*)=2 L B/R_*^3>0.                         (21)

This is a genuine positive leading balance, but it is not yet a theorem for
the charged Cao branch. The full proof must differentiate the exact
Euler--Maxwell energy on a fixed-label admissible path, or derive its direct
finite-core energy asymptotic. It must include the comoving magnetic/current
contribution, use the actual tag profile, and track the dependence of `s/R`
and both logarithms on `R`. The subluminal finite window
and core admissibility may exclude (21), especially as `Q` tends to zero.
The 0085 strong saddle also prevents this one radial minimum from becoming a
full P2 coercivity claim.

Even if the balance survives, `I_z` is translational impulse, not the internal
oscillator action `J`. It can fix a carrier scale entering mode frequencies
and KKS normalization, but an additional same-carrier calculation must derive
the map to `J`. Moreover (21) depends continuously on `Q`, `kappa`, and the
Maxwell normalization; it is not a universal quantum action. Route E0
establishes the leading balance. Route E1 is blocked at the exact finite-core
energy/constraint calculation.

## 8. Coupled route and strongest continuation

If a compact defect eventually supplies a signed integer character, it must
enter the same reduced energy (1).  Positive oscillator and Coulomb terms
still select `J=0`.  A nonzero level requires an earned negative detuning,
binding term, chemical potential, anomaly-inflow term, or discrete admissible
set whose coefficient is fixed by the carrier rather than declared.  Each such
term must also preserve the exact charge row and survive (2)--(5).

The most useful next calculation is therefore not another dimensional
identity.  It is the actual constrained third/fourth variation and
principal-value two-mode self-energy on a source-specific invariant or
center-stable construction, run in parallel with a compact-defect supplier
whose symplectic coefficient is derived.  P2 remains a separate active route:
a prepared finite-time doublet and a radiation/readout law are not orbital
persistence.  Born weights, reset, exchange/fermionic character, and the
neutrino sector remain independent parent achievements.
