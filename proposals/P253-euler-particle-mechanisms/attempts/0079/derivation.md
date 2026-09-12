# Distinct-sector Euler doublet and control calculation

## 1. Exact symmetry count

Work in the translating, nonrotating frame of one oriented Cao ring.  The
commuting spatial symmetry is axial rotation,

    (R_alpha q)_n=exp(i n alpha)q_n,     n in Z.          (1)

Meridional reflection is a reverser of the oriented translating flow rather
than a second commuting continuous symmetry.  Real conjugation maps the
`n`-fiber to the `-n`-fiber and reverses the spectral sign as required by the
real Hamiltonian system.  Hence every complex irreducible representation of
the commuting `SO(2)` is a one-dimensional character.  The real cosine/sine plane for
fixed `|n|` is the realification of this one complex line:

    Re(e_n exp(in theta)),  Im(e_n exp(in theta)).        (2)

It is one oscillator, not two complex modes.  Reversibility supplies the
positive/negative-frequency partner, not a second positive oscillator.

The radial Sturm index `J` gives multiple copies of the same axial character;
the `SO(2)` representation alone neither fixes that multiplicity nor forces
two different copies to have equal frequency.
Therefore the actual symmetry representation of a single Cao ring supplies
no symmetry-enforced complex doublet.

**Route-A verdict:** refuted for the one-ring Cao symmetry group.  The named
mechanism is that axial `SO(2)` has one-dimensional characters and supplies no
equality constraint among their radial copies; neither (2) nor the `+n/-n`
reality pair supplies the required second mode.  This does not refute a distinct-
sector tuned crossing or a different carrier with a genuine internal
symmetry.

## 2. Same-`k` Sturm ordering and the massive distinct-sector crossing

For the compact-column `m=0` modes, write `mu=k^2` and normalize

    H_mu u_J=lambda_J(mu) Phi u_J,
    H_mu=H_0+mu,
    integral Phi |u_J|^2 r dr=1.                         (3)

The eigenvalues are simple and strictly ordered.  Hellmann--Feynman gives

    d lambda_J/dmu=integral |u_J|^2 r dr=:a_J(mu)>0.     (4)

Since the positive form of `H_0` is

    h_0[u]=integral (|partial_r^*u|^2)rdr>0,

one has

    lambda_J-mu lambda_J'
       =h_0[u_J]>0,
    d/dmu [sigma_J(mu)^2]
       =d/dmu[mu/lambda_J(mu)]
       =h_0[u_J]/lambda_J^2>0.                           (5)

Thus every branch frequency is strictly increasing in `|k|`, while at the
same `k`

    J_1<J_2 implies lambda_(J_1)<lambda_(J_2)
                  and sigma_(J_1)>sigma_(J_2).           (6)

No scalar Cao-family parameter can create a crossing of two same-`k` radial
labels without first destroying Sturm simplicity/order.  The only viable
one-ring target in the present family is therefore

    sigma_(J_1,n_1)(mu_*)=sigma_(J_2,n_2)(mu_*),
    n_1!=+/-n_2,                                        (7)

or a genuinely different polarization/azimuthal sector.

The transfer theorem supplied by reviewed 0074/0078 is a **massive** theorem:
`k=delta n` remains in one fixed positive compact interval while
`|n|=Theta(delta^(-1))`.  It does not cover taking `delta` to zero with `n`
fixed.  The admissible column candidate must therefore start with two fixed
positive wave numbers `k_1,k_2` in that interval and solve

    F_col(k_1,k_2)
      =sigma_(J_1)(k_1)-sigma_(J_2)(k_2)=0,
    k_1!=k_2.                                            (8)

Because every branch is continuous and strictly increasing by (5), such a
column equality exists whenever the two frequency images over the chosen
massive interval overlap.  This is a checkable interval condition, not a
consequence of symmetry.  For a thin-ring sequence choose integer
approximants

    n_a(delta)=round(k_a/delta),
    k_(a,delta)=delta n_a(delta)=k_a+O(delta).            (9)

An arbitrary pair of integer approximants gives only an approximate equality.
The rational-ray construction below instead arranges an exact column equality
and uses the **uniform `C^0` spectral convergence established by reviewed
0074/0078**. No first curvature coefficient is required for existence.

### 2.1 Rational-ray exact-crossing continuation

There is a parameter-free version of alternative (b).  On any interval where
the frequency images overlap, strict monotonicity makes `sigma_(J_2)`
invertible.  Define the analytic matching ratio

    r_12(k)=sigma_(J_2)^(-1)(sigma_(J_1)(k))/k.          (10a)

If `r_12` is nonconstant, it has a subinterval on which `r_12'!=0`.  Density
of the rationals then gives coprime positive integers `p,q` and `k_*` with

    r_12(k_*)=q/p=:r,
    F_r(k):=sigma_(J_1)(k)-sigma_(J_2)(r k),
    F_r(k_*)=0.                                         (10b)

Differentiating the defining relation for (10a) gives the exposing identity

    F_r'(k_*)
      =k_* r_12'(k_*) sigma_(J_2)'(r k_*)!=0.           (10c)

Thus the rational ray is a simple exact column resonance, not a Diophantine
near miss.  Nonconstancy itself is a source-defined scalar test.  The
small-`k` coefficients below give

    r_12(k)=A_1/A_2
      +[A_1^3(a_(2,0)-a_(1,0))/(2A_2)]k^2+O(k^4),       (10d)

so `a_(1,0)!=a_(2,0)` is one sufficient exact witness.  If this coefficient
vanishes, the analytic Sturm recurrence continues to the first unequal even
coefficient; equality of every coefficient would be the genuine homothety
obstruction.  This reduces selection of `p,q` to enclosing two values of the
one-dimensional analytic function (10a), rather than searching a two-fiber
spectrum numerically.

For each large integer `N`, now set

    n_1=Np,  n_2=Nq,
    delta_N^0=k_*/(Np).                                 (10e)

Then `delta_N^0 n_1=k_*`, `delta_N^0 n_2=r k_*` and the column mismatch
vanishes exactly while `delta_N^0->0`.  Audit of 0074 equations (15)--(20),
(27)--(30), and (37b)--(39) shows that its graph error is uniform when
`k_delta` ranges over a fixed compact neighborhood `K_*` of `k_*` and
`r k_*`.  On the two isolated simple contours, the resolvent identity therefore
gives a uniform eigenvalue error

    e_N=sup_(|k-k_*|<=h_0) max_(a=1,2)
      |sigma_(J_a,N p_a)^Cao(k/(Np))
        -sigma_(J_a)^col((p_a/p)k)| ->0,                (10f)

where `p_1=p`, `p_2=q`.  This uniform error does **not** by itself provide a
connected selection of exact carriers.  The IVT additionally requires a
centered carrier path covering the displayed `delta` interval, continuous in
the coefficient/Hodge graph topology; simple Riesz eigenvalues are then
continuous along that path.  Cao's existence-plus-uniqueness theorem is
pointwise in its small parameter and 0074 constructs parameter-uniform
charts, but neither record proves local coverage/continuity in the geometric
parameter `delta`.  This can be earned either by centered compactness plus
uniqueness and local openness of the source-parameter-to-`delta` map, or by
the uncharged steady Green-map/HSE branch being constructed in P253/0080.
The clean 0080 interface holds circulation `kappa` and mean radius `R` (or the
equivalent rows fixing the same limiting column) fixed, uses external core
scale `epsilon_c` as parameter, and solves the bordered equations for
chemical potential and speed. It must prove

    s_(epsilon_c)/R
      =C(kappa,R)epsilon_c[1+o(1)],
    d(s_(epsilon_c)/R)/d epsilon_c!=0,                  (10f.1)

with a continuous/smooth graph map.  Then geometric `delta` covers (10g)
while `F_r` remains the same column function and residual profile variation
is absorbed in `e_N=o(h_N)`.  If instead `R` or speed is varied, its induced
column-profile derivative enters the transversality; noncancellation with
`F_r'` must be proved rather than assumed.

Take a decreasing envelope `bar e_N>=e_N`, `bar e_N->0`, and choose

    h_N=sqrt(bar e_N)  (or h_N=N^-1 if bar e_N=0),
    delta_N^+/-=(k_*+/-h_N)/(Np).                       (10g)

Then `e_N=o(h_N)`.  By (10c),

    F_r(k_*+/-h_N)=+/-F_r'(k_*)h_N+o(h_N),              (10h)

so, **on such a continuous covering carrier path**, the actual Cao mismatches
at `delta_N^-` and `delta_N^+` have opposite signs for all large `N`.
Continuity and the intermediate value theorem give an **exact Cao crossing**

    delta_N in (delta_N^-,delta_N^+),
    delta_N=delta_N^0+o(N^-1).                          (10i)

Both `delta_N n_a` stay in the same fixed massive neighborhoods.  This uses
`delta` as the steady-family parameter and does not count its static choice as
a gate.  A stronger uniform `C^1` expansion would sharpen the location to
`O(N^-2)`, but it is not an existence premise.

This uncharged sequence is not automatically a charged Euler--Maxwell
sequence.  The Cao translation speed grows like `W log N` along
`delta_N=Theta(N^-1)`.  At fixed electromagnetic light speed `c_EM`, a
subluminal charged continuation can use only

    N<=N_EM approximately exp(c_EM/W)                   (10i.1)

up to the fixed source constants and safety margin.  It must exhibit at least
one integer satisfying

    max(N_graph,N_IVT,N_response)<N<N_EM.               (10i.2)

No `N->infinity` Euler crossing theorem is transferred wholesale to charged
Euler--Maxwell without this hierarchy.

### 2.2 A forced rational ray from a different polarization

The `m=0`/`m=0` ratio test (10d) can be bypassed at column level by pairing
one `m=0` Kelvin mode with the `m=1` translation-bending branch already
isolated in reviewed 0074/0078.  Its exact small-`k` determinant gives

    tau_1(k):=-b_1(k)=C_b k^2 log(1/(ka))+O(k^2),
    C_b=F(a)/2>0,                                       (10j)

whereas (11a) below gives `sigma_J(k)=A_J k+O(k^3)`.  For a fixed rational
`r>0`, set

    F_(J,r)^pol(k)=sigma_J(k)-tau_1(rk).                (10k)

For every fixed `r`, (10j) implies `F_(J,r)^pol(k)>0` for sufficiently small
positive `k`.  Choose one sufficiently small fixed `x_0` in the controlled
range of (10j).  Since `sigma_J(x_0/r)->0` as `r->infinity` while
`tau_1(x_0)>0`, every sufficiently large rational `r=q/p` satisfies

    F_(J,r)^pol(x_0/r)<0.                               (10l)

Choose a sign-changing boundary root `k_r in (0,x_0/r)`, for example the
supremum of the connected positive component issuing from zero.  There are
then fixed values `k_r^-<k_r<k_r^+` with strict opposite signs.  The leading
value balance localizing this root is

    r^2 k_r log(1/(r k_r))=A_J/C_b+o(1),                (10m)

but no derivative is inferred from the value-level `O(k^2)` remainder in
(10j).  Fix this rational `r`, take `n_1=Np,n_2=Nq`, and use the physical
endpoints `delta_N^+/-=k_r^+/-/(Np)`.  Uniform `C^0` graph/eigenvalue
convergence preserves their strict signs for all large `N`; the continuous
exact-carrier path then gives an actual Cao crossing by IVT.  Thus the source
asymptotics force a rational-ray crossing without assuming
`a_(1,0)!=a_(2,0)` or a `C^1` remainder.

This route has one new physical obligation: reviewed 0074/0078 establishes
the nonzero bending branch and the positive Krein sign of the selected `m=0`
mode, but not the actual constrained Hessian/KKS sign of the `m=1` bending
branch.  That sign and its absolute normalization must be derived before the
polarization pair is a positive doublet.  The route nevertheless establishes
an exact sign-changing column crossing rather than merely a candidate ratio.

There is a formally attractive but currently **uncovered** small-`k` route.
For completeness, perturbing (3) at `k=0` gives

    lambda_J(k)=lambda_(J,0)+a_(J,0)k^2+O(k^4),
    sigma_J(k)=A_J k+B_J k^3+O(k^5),                   (11a)

    A_J=lambda_(J,0)^(-1/2),
    B_J=-a_(J,0)/(2lambda_(J,0)^(3/2))<0.

For fixed distinct integers it would give

    F_col(delta)
      =delta(A_1|n_1|-A_2|n_2|)
       +delta^3(B_1|n_1|^3-B_2|n_2|^3)
       +O(delta^5 max|n_a|^5).                          (11b)

If

    (A_1|n_1|-A_2|n_2|)
    (B_1|n_1|^3-B_2|n_2|^3)<0,                         (11c)

the truncated column equation has a small positive simple root; continued-
fraction approximants can make the linear mismatch small.  But this fixed-`n`,
`k->0` construction lies outside the massive interval and is **not supplied by
0074**.  It becomes executable only after a joint `k->0`, `delta->0`
Hodge/graph transfer with the zero-shift logarithmic channel controlled.

**Route-B degeneracy verdict:** conditionally reduced to two exact inputs: one
rational regular value of the source-defined analytic ratio (10a), and one
continuous exact-carrier path covering (10g).  Given both, 0074's uniform
`C^0` graph-Riesz transfer proves exact Cao crossings (10i); no `C^1`
curvature expansion is needed.  The remaining Route-B control construction
is blocked at physical KKS response and autonomous gate histories, not at a
two-fiber curvature coefficient. Equations (11a)--(11c) give an exact
nonconstancy test without claiming a separate fixed-`n` torus transfer.

## 3. Exact Euler perturbation and harmonic selection rule

Let the steady carrier vorticity be `omega_0`, its whole-space velocity be
`u_0=B omega_0`, and let the relative generator be

    A(omega_0)q=-[u_0-c_0 e_z,q]-[Bq,omega_0].           (13)

For a physical constrained perturbation `h`, differentiation on the fixed
translating slice gives the exact first interaction

    V_h q=-[Bh,q]-[Bq,h],                                (14)

and in particular

    V_h e_2=-[Bh,e_2]-[B e_2,h].                        (14a)

The two brackets in (14) respectively contain transport by the induced global
velocity and stretching against the perturbed vorticity.  `B` is the decaying
whole-space Hodge map, so pressure, collar, interface and exterior response are
not omitted.  The allowed perturbation space below already imposes the
physical center, impulse, circulation and phase rows; they are not replaced by
an unproved additive finite-row operator.  A derivative along a separately
parameterized family with varying translation speed requires its own bordered
row and is not silently folded into (14).

If `h` has toroidal harmonic `ell`, then

    V_(h_ell): X_n -> X_(n+ell).                         (15)

For a real perturbation `h=h_ell+conjugate(h_ell)`, projection from mode `e_2`
in `n_2` to `e_1` in `n_1` is therefore possible exactly when

    ell=n_1-n_2  or  ell=n_1+n_2                         (16)

according to whether the target is paired with `e_1` or its real partner.
Equation (16) is the source-bearing selection rule.  A generic Fourier label
is not enough: the physical matrix element is

    M_12(h)=Omega_KKS(e_1^#,V_h e_2),                   (17)

where `e_1^#` is the KKS/Hessian dual on the Riesz range.

Let `Y_ell^s` be the constrained physical interaction space: the closure in
the declared DA/Hodge graph norm of smooth compact perturbations with toroidal
character `ell` and the exact center, impulse, circulation and phase rows.
The bracket/Hodge mapping makes

    G_12 : Y_ell^s -> C,
    <G_12,h> := Omega_KKS(e_1^#,V_h e_2)                 (17a)

a continuous dual functional.  Equation (17a) is its definition at the
present scope.  No expanded `ad/ad^*` formula, pointwise density, or projected
distributional representative is asserted.  A compact smooth DA seed with
`M_12!=0` exists exactly when this continuous functional is nonzero on the
dense smooth sector.  Harmonic selection and unique continuation of the
individual modes do **not** rule out cancellation between the two terms in
(14a), so nonvanishing remains open.

For the polarization pair in Section 2.2, toroidal (16) is only one row.  In
the limiting action-angle column, an `m=0` mode and an `m=1` mode require the
interaction seed to have streamline/meridional harmonic `Delta m=+/-1` (with
the corresponding sum channel for the real partner).  Curvature destroys
exact column `m`-diagonality, so on the actual ring this requirement is imposed
on the **full** constrained functional (17a), including every
curvature-mixed meridional row.  A seed with the correct toroidal `ell` but
zero constrained `Delta m=1` response does not mix the polarization pair.

Axisymmetric compact DA seeds preserve each `n_a` and give a diagonal
compression.  Their differential response is likewise

    <G_3,h_0>
       :=Omega_KKS(e_1^#,V_(h_0)e_1)
         -Omega_KKS(e_2^#,V_(h_0)e_2),                  (18)

where `G_3` is only the corresponding continuous dual functional on `Y_0^s`.
Distinct nodal data make `G_3` a concrete candidate but do not alone prove it
is nonzero on the constrained space.  Evaluating `G_12` and `G_3` is therefore the exact
physical-interaction dependency.  If either is nonzero, choosing a compact
smooth test seed in its defining dense constrained sector constructs the
required instantaneous DA direction.

## 4. Complex-linear and squeezing decomposition

Assume for this section that the missing equality (7) has been constructed and
let `P` be its real four-dimensional positive KKS projector.  With the
compatible complex structure `J`, `J^2=-1`, define

    V_C=(V-JVJ)/2,       V_A=(V+JVJ)/2.                  (19)

Direct multiplication gives

    [V_C,J]=0,           {V_A,J}=0.                     (20)

Thus `V_C` acts complex linearly on the two positive modes, whereas `V_A`
couples them to their negative-frequency partners.  In an energy-action
orthonormal basis,

    H_C=-J P V_C P=h_0 1+h dot sigma.                    (21)

Conditional on `G_3!=0`, the axisymmetric seed in (18) supplies a nonzero
`sigma_3` component. Conditional on `G_12!=0`, a complex seed `h_ell` with
`c=M_12(h_ell)!=0` supplies an off-diagonal row. Its
real cosine and sine quadratures compress, after one basis phase choice, to

    H_cos^0=|c| sigma_1,
    H_sin^0=|c| sigma_2,
    h_cos cross h_sin=|c|^2 e_3!=0.                     (22)

This is an algebraic non-collinearity theorem, not yet two controls.  Axial
rotation conjugates the quadratures:

    h_sin=R_(pi/(2|ell|))*h_cos.                         (23)

If (23) is implemented only by hand-selecting the phase of one seed, by a
rigid rotation of one static perturbation, or by waiting for its free phase,
there is only one autonomous control axis.  Two axes are earned only by two
independently prepared Euler subsystems/histories or by a proved internal
control mechanism that switches the quadratures without leaving the common
invariant neighborhood.

The pair `(h_0,h_ell)` is already non-collinear whenever the predicates
`G_3!=0` and `G_12!=0` hold;
it avoids counting the cosine/sine phase orbit twice.  What remains is to make
both seeds coherent physical interactions for the full gate time.

### 4.1 Normalized high-harmonic response and the useful-gate criterion

For the rational ray, the required mixer order is

    ell_N=n_1-n_2=N(p-q)=Theta(N).                      (18a)

Let `Y_(ell_N)^s` be the physical compact DA interaction graph with the full
Hodge and finite rows, and define the normalized response

    g_N^(s)=sup{|M_12(h)|:
                 h in Y_(ell_N)^s, ||h||_Y=1}
           =||G_(12,N)||_((Y_(ell_N)^s)^*).             (18b)

This is the quantity relevant to a gate, not bare nonvanishing of `G_12,N`.
The Fourier normalization from 0074 gives physical pure harmonics a factor
`(R_N a_c^2)^(-1/2)`, where `a_c` is the core scale.  But 0074 leaves open the
absolute Sturm-to-physical Hessian/KKS conversion, and `M_12` contains the
KKS dual `e_1^#`.  Therefore the following bound is **conditional** on uniform
physical graph norms for `e_1,e_2,e_1^#`, including density and real/complex
factors, and on finite-row representatives satisfying the same uniform
estimate.  Under precisely that hypothesis, the two bilinear brackets in
(14), the order-minus-one Hodge map and one toroidal integration give

    g_N^(s)<=C_s(R_N a_c^2)^(-1/2).                     (18c)

For fixed `a_c` and `R_N=a_c/delta_N=Theta(N)`, this conditional bound is
`O(N^(-1/2))`.
There is no additional `N^s` loss in the **physical covariant** norm because

    delta_N |ell_N| -> k_*|p-q|/p;                     (18d)

all angular derivatives enter through the finite physical wave number.  By
contrast, an unscaled fixed-reference `H^s(dtheta)` norm costs
`|ell_N|^s`; using that topology would weaken the dual response by the
corresponding power.  Equations (18c)--(18d) distinguish geometric volume
dilution from a coordinate-artifact derivative loss.  They are upper bounds;
a useful interaction still needs a lower bound from the actual constrained
dual functional.

If a coherent normalized history has amplitude
`A_N=sup_t||h_N(t)||_Y`, the complex-linear rotation angle obeys

    Theta_N<=A_N g_N^(s) T_g.                           (18e)

A `pi/2` gate therefore requires

    A_N g_N^(s)T_g>=pi/2.                               (18f)

The useful supplier criterion is (18f) together with

    eps_A<=tol_A,
    eps_out S_A/(1-kappa_back)<=tol_Q,
    S_A^2/(1-kappa_back)^2-1<=tol_J,                    (18g)

together with the isolated inverse bound used in (28b) and a nonlinear
existence/remainder estimate through that same `T_g`.  If the
upper scale in (18c) is sharp and `A_N` stays inside a fixed small carrier
neighborhood, then `T_g` is at least `Theta(N^(1/2)/A_N)`.  Standard local
well-posedness near the nonzero steady background supplies only a time
depending on the full background-plus-perturbation norm; no perturbative
`Theta(A_N^(-1))` lifespan or uniform coverage of this growing gate time is
proved here.  Thus an exact doublet plus `G_12,N!=0` is still not a P4-useful
analyzer without a long-time coherent-history or normal-form mechanism.

## 5. Exact finite-time squeezing, action drift and two-sided Q leakage ledger

Let `z in C^2` be the positive-mode amplitude and let `w=Qq`.  In the rotating
frame of the common frequency, the exact linear block has the form

    dot z=C(t)z+D(t)conjugate(z)+P V(t)Qw,
    dot w=A_Q(t)w+Q V(t)Pz,                              (24)

where `C=P V_C P` is skew-Hermitian in the energy-action metric and
`D=P V_A P` is the squeezing block.  Let

    eps_A=integral_0^Tg ||D(t)||dt,
    S_A=exp(eps_A).                                     (25)

Let `U_P(t,s)` be the propagator of the isolated internal block
`dot z=Cz+D conjugate(z)`, so `||U_P(t,s)||<=S_A`.  With `M_Q(t,s)` the actual
nonnormal complement propagator bound, define separately

    eps_out=sup_(t<=Tg) integral_0^t
      M_Q(t,s)||Q V(s)P||ds,                             (26a)

    kappa_back=sup_(t<=Tg) integral_0^t
      ||U_P(t,r)|| ||P V(r)Q||
      [integral_0^r M_Q(r,s)||Q V(s)P||ds]dr.            (26b)

The standard comparison system for `(z,conjugate z)` gives

    ||beta_Bog(Tg)||<=sinh(eps_A),                       (26c)

but (26c) applies to the **isolated** `PVP` propagator.  For the coupled system
with `w(0)=0`, Duhamel substitution of the `Q` equation into the `P` equation
gives

    Z<=S_A||z(0)||+kappa_back Z,
    W<=eps_out Z,                                       (27)

where `Z=sup_(t<=Tg)||z(t)||` and `W=sup_(t<=Tg)||w(t)||`.  Therefore, when
`kappa_back<1`,

    Z<=S_A/(1-kappa_back)||z(0)||,
    W<=eps_out S_A/(1-kappa_back)||z(0)||.               (28a)

Let `z_*(t)=U_P(t,0)z(0)`.  The inverse isolated propagator obeys the same
comparison bound, and (27) gives

    S_A^(-1)||z(0)||<=||z_*(t)||<=S_A||z(0)||,
    ||z(t)-z_*(t)||
      <=[kappa_back S_A/(1-kappa_back)]||z(0)||.         (28b)

Consequently the genuine two-sided absolute action estimate is

    | ||z(Tg)||^2-||z(0)||^2 |/||z(0)||^2
      <=S_A^2/(1-kappa_back)^2-1+R_NL(Tg).               (28c)

where `R_NL` is the separately estimated quadratic Euler remainder.  The
right-hand side follows by adding the isolated absolute action change
`S_A^2-1` to
the return correction
`S_A^2 kappa_back(2-kappa_back)/(1-kappa_back)^2`.
Equations (25)--(28) retain gate time and both directions of leakage.
`eps_out` alone controls only `w` generated from `z`; without `PVQ` and
`kappa_back` it cannot imply (28c).  The lower control in (28b), not a one-way
leakage estimate, is what prevents an unmeasured contraction.  Replacing `V` by
`epsilon V` while taking `Tg=Theta(epsilon^(-1))` does not make these errors
small unless their integrated values shrink.

If `P V_A P=0` by an exact harmonic or reversibility rule, then `eps_A=0`.
Otherwise the internal Bogoliubov estimate (26c) and the two-sided return
estimate (26b) are both required.  A spectral gap alone bounds neither
`eps_out` nor `kappa_back`; they need the graph-resolvent/semigroup bound and
both interaction blocks, including all other modes, pressure, interface and
finite rows.

The algebra in Sections 3--5 is importable and is checked by the exact oracle
in this attempt. Its physical conclusion is limited twice: the full-Hodge
response functionals `G_12,G_3` still require evaluation, and any resulting
compact smooth seeds must then be kept coherent and independently switchable
through the `Tg` required for a finite rotation. Local Euler well-posedness on
a fixed time interval does not automatically cover `Tg=Theta(1/||V||)`.

## 6. Two-carrier route

Two identical copies carry equal uncoupled mode frequencies, but the exact
whole-space Biot--Savart interaction compresses to

    H_pair=omega_*1+t(L)sigma_1+Delta(L)sigma_3+...,     (29)

where exchange symmetry sets the diagonal entries equal only in the exactly
symmetric configuration; it does not force the off-diagonal transfer
`t(L)` to vanish.  The symmetric and antisymmetric combinations therefore
split by `2|t(L)|`.  A far-field multipole calculation is required to decide
whether `t(L)` can vanish in a permitted orientation; generic nonvanishing is
not assumed. Infinite separation is not a localized interacting carrier and
cannot be used as exact degeneracy.

No reviewed construction currently supplies a two-Cao-ring relative
equilibrium with one parameter at which `t(L)=0`, nor two independently
switchable relative motions with a common restoring neighborhood.

**Route-C verdict:** blocked at the actual two-carrier relative equilibrium
and zero-splitting/transversality construction.  Exchange symmetry alone is
not a degeneracy theorem.

## 7. Failure-derived Route D: autonomous resonant Floquet conversion

The static crossing can be replaced, without pretending the frequencies are
equal, by a third physical Euler mode `h_ell(t)` satisfying

    ell=n_1-n_2,
    nu_ell=sigma_1-sigma_2.                              (30)

In the interaction picture, its resonant term is constant and (17) gives a
mode-conversion matrix.  Harmonic sign or the presence of a conjugate pump
does **not** decide whether a term is complex-linear or squeezing.  The
physical derivative must first be resolved with the positive/negative-
frequency projectors,

    P_+ V P_+,  P_+ V P_-,  P_- V P_+,  P_- V P_-,      (31)

and then decomposed by (19).  Only the resulting `P V_A P` block is charged to
`eps_A`; off-resonant complex-linear rows remain in `V_C` and require their own
oscillatory/Floquet estimate.  This is a genuine three-wave/Floquet route, not
a static avoided crossing.

The current Cao inputs provide infinitely many positive Kelvin frequencies
but no exact three-sector resonance theorem and no nonlinear pump-depletion/
back-reaction bound.  The next executable construction is to combine the
strict monotonicity (5), integer harmonic rule (16), and a two-parameter Cao
or companion-carrier family to solve (30), then apply the full graph
semigroup estimate to its resonant normal form.

**Route-D verdict:** blocked with the exact missing objects named: one physical
auxiliary eigenmode satisfying (30), a nonzero trilinear KKS coefficient, and
an autonomous gate-time remainder bound after the projection (31).  It remains
an active failure-derived alternative and is not an obligation-level no-go.

## 8. Failure-derived Route E: fixed mixer harmonic and high radial index

The volume dilution and long gate in (18c)--(18g) motivate a distinct route
that keeps the mixer harmonic `ell` fixed.  Use adjacent large radial indices
and choose

    J_1=J,       J_2=J+1,
    n_1=J ell,   n_2=(J+1)ell,
    delta_J=k_*/(J ell).                                (32)

Then `n_2-n_1=ell` while

    k_1=k_*,     k_2=k_*(1+J^(-1)).                    (33)

The limiting `m=0` Sturm supplier is the 0066 calculation independently
retained by 0073 (and traces to the 0048/0053 reviewed pencil).  Its exact
fixed-`k`, large-`J` scope is

    sigma_J(k)=k L_Phi/(pi J)+O(k/J^2)                 (34)

makes the two leading terms exactly equal:

    k_1/J=k_2/(J+1)=k_*/J.                             (35)

Thus the adjacent radial gap `Theta(k/J^2)` is balanced by the fixed-harmonic
wave-number shift `partial_k sigma * delta_J ell=Theta(J^(-2))`.  The first
unequal phase coefficient in the two-term Sturm expansion decides the signed
crossing; this is not the same rational-ray high-`ell` response problem.

Neither 0066/0073 nor the current 0074 transfer gives the joint regime.  The
latter fixes `J` before `delta->0`; the former is a column-sector Weyl law.
They do not give `J=Theta(delta^(-1))`, and the transferred `o(1)` Riesz error
is much larger
than the `Theta(J^(-2))` contour separation.  Route E therefore requires a
two-term uniform high-`J` Sturm expansion and a thin-ring graph-Riesz error
`o(J^(-2))`, including physical KKS/DA normalization of these oscillatory
radial modes.  A finite-denominator variant is useful only if `ell` remains
bounded while `delta` is below this `J`-dependent transfer threshold.

**Route-E verdict:** blocked at the joint high-`J`/thin-ring Riesz theorem and
the signed next Weyl coefficient.  Its reward is a genuinely bounded mixer
harmonic, so it remains the primary alternative if the normalized high-`ell`
response in (18b) decays too rapidly for (18f)--(18g).

## 9. Strongest exact verdict

0079 establishes:

- the one-ring one-dimensional-character obstruction to symmetry-enforced
  equality among radial copies and the same-`k` Sturm ordering;
- the massive rational-ray bracketing theorem: reviewed uniform `C^0`
  graph-Riesz transfer plus one continuous carrier path gives exact Cao
  crossings at `delta_N=delta_N^0+o(N^-1)`;
- an unconditional sign-changing rational column crossing between an `m=0` Kelvin
  mode and the `m=1` translation-bending magnitude, with the bending Krein
  sign still open;
- the full-Hodge Euler operator derivative and harmonic selection rule;
- the exact continuous constrained response functionals whose nonvanishing is
  equivalent to compact physical DA seeds with diagonal and off-diagonal
  instantaneous compression;
- the complex-linear/anti-complex-linear split, Pauli non-collinearity
  criterion, and a genuinely two-sided gate-time squeezing/action/leakage
  ledger; and
- the normalized high-harmonic response/gate-time criterion and a distinct
  fixed-`ell`, high-radial-index matching route; and
- why exchange symmetry or phase-rotating one seed does not by itself create
  a doublet or two controls.

It does not yet construct an exact equal-frequency **positive** Cao doublet:
the `m=0` pair still needs a rational nonidentity witness and a continuous
geometric-`delta` carrier bracket, while the forced polarization crossing
needs that bracket and the bending Krein sign. It also does not prove the
normalized response functionals nonzero or construct two independently
switchable autonomous Euler controls. This is a route gap, not a physical
impossibility. The strongest next construction joins the fixed-`kappa,R`
0080 carrier path to (10g), evaluates the polarization KKS sign and physical
dual normalization, and then obtains lower bounds for `G_12,G_3` over the
full gate time. Route D and the fixed-`ell` Route E remain distinct
alternatives.
No P2, P4, stable particle, quantum measurement or relativity claim follows.
