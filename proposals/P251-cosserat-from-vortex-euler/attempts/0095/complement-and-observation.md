# Full Kelvin complement, conditional affine closure and physical observations

## 1. The original premise and the exact distinction being repaired

[Issue #198](https://github.com/vantasnerdan/substrate-framework/issues/198)
asks for an exact **conditional** Euler coarse-graining with named,
independently falsifiable ensemble premises. Proposal N3 includes affine
transport and explicitly excludes nonaffine relaxation. Consequently an
exact pullback of Euler's action under a specified Cauchy–Born restriction
is a legitimate conditional constitutive result. Requiring that every
unrestricted microscopic trajectory remain on a finite nonlinear ansatz
would change that objective and is not imposed here.

The observed distinction nevertheless has physical content. Freezing an
excluded nonaffine coordinate, freely relaxing its energy, and retaining
its actual Euler/Kelvin dynamics are three different operations. They can
give different inertia and physical field maps. The formulas below compute
that difference explicitly. They preserve the established conditional
constitutive result; they prevent it being silently upgraded to a claim
about unrestricted Euler reconstruction or unchanged material inertia.

## 2. Perform the actual Kelvin reduction before selecting phase coordinates

On the full Euler coadjoint tangent space modulo relabeling stabilizers,
write the KKS form as `Omega` and the stationary energy Hessian as `H`.
For the present Beltrami field, 0084/0091 give the actual reconstruction
operator

    A eta=(lambda P-curl)(eta cross u),
    Omega(A eta,zeta)=H(eta,zeta).

The convention `Omega(Q,S)=B` agrees with `L=B p qdot-H/2`: in matrix
notation `Omega zdot=-H z`, hence `A=-Omega^-1 H`. Stabilizers satisfy
`curl(eta cross omega)=0`; they change the isovortical generator without
changing the Eulerian vorticity. These, not arbitrary omitted directions,
are the genuine relabeling gauge.

Let `E` embed the finitely retained affine, physical-angle and reaction
directions, including the conjugate phase variables. Assume its finite KKS
matrix `Omega_E=E*Omega E` is invertible, as explicitly constructed for
0085's pair. No inverse of the full weak infinite-dimensional KKS form
is presumed by the finite projection

    Pi=Omega_E^-1 E*Omega,    Pi E=I,
    eta=E z+r,    E*Omega r=0.

The full Euler generator `A` is supplied by its differential/Leray formula
on the compatible Sobolev domain. The phase-space complement is real
nonaffine motion, except where its vorticity variation is a stabilizer.

The restricted action gives

    zdot=A_E z,    A_E=Pi A E=-Omega_E^-1 H_E.

Its exact microscopic reconstruction residual is

    R(z)=E A_E z-A E z=-(I-E Pi)A E z.                     (1)

It obeys `E*Omega R=0`, but that does not imply
`curl(R cross omega)=0`. Equation (1) is precisely the force into the
excluded nonaffine complement. It is not a new all-wavelength condition.
The original no-nonaffine-relaxation premise licenses setting that
coordinate aside in the conditional variational model; it does not
identify the discarded force with a particle-relabeling gauge.

## 3. Constructive repair: retain the complete Kelvin complement

On a compatible linear phase-space chart, keep `r` rather than setting it
to zero. The exact projected equations are

    zdot=A_E z+A_ER r,
    rdot=A_RE z+A_RR r,

with `A_ER=Pi A` on the complement and `A_RE=(I-E Pi)A E`. The exact
finite-time solution is

    r(t)=exp(t A_RR)r(0)
             +integral_0^t exp((t-s)A_RR) A_RE z(s) ds.      (2)

For a time-dependent chart the corresponding propagator replaces the
exponential and its connection terms remain. Formula (2), inserted into
the first equation, is a complete same-Euler memory equation and retains
the actual initial nonaffine state. No arbitrary dissipative or fitted
kernel has been introduced. For an infinite-dimensional complement the
formula uses its actual well-posed evolution operator; it is not an
assertion that an arbitrary formal inverse is bounded.

When stationary harmonic/resolvent analysis is admissible, set
`D(omega)=H-i omega Omega`. The symplectically orthogonal split has
`D_ER=H_ER`, and exact elimination yields

    D_eff=D_EE-H_ER D_RR^-1 H_RE,
    D_RR=H_RR-i omega Omega_RR.                             (3)

An inhomogeneous term from the initial/free complement is retained for the
initial-value problem. Formula (3) is a boundary/resolvent action identity
where the inverse exists; a causal initial-value Green operator is not
silently treated as a self-adjoint time-symmetric action kernel. Equation
(2) supplies the initial-value repair without requiring a spectral gap.

This distinction answers the 0091 mass question: a positive material
metric `M(Z)=A(Z)` is not automatically the reduced Euler inertia. The
full Kelvin/relabeling operation supplies the complement and its actual
symplectic/Hamiltonian coupling. A freely relaxed static Schur inverse is
only the `omega=0` instance of (3), not a universal dynamical replacement.

## 4. Physical angle, spin, tag and mean maps are eliminated with the action

Let an actual linear observation be `o=O_E z+O_R r`: this may be a
smoothed core-vorticity angle, mechanical angular momentum, a finite
material centroid, or a derivative of 0087's specified mean. In the same
resolvent calculation its map is

    O_eff=O_E-O_R D_RR^-1 H_RE.                             (4)

The initial/free-complement contribution remains as in (2). Applying a
Schur correction to energy but leaving a changed physical observation
unchanged is not this reduction. Equations (3) and (4) are one field-map
theorem; their coefficients come from the same `H`, `Omega` and actual
observable rows.

The material observations themselves are constructed, not inferred from
a canonical name. Given the actual Eulerian velocity perturbation
`v=P(eta cross omega)`, the Lin reconstruction is

    xi_t+[u,xi]=v,
    xi(t)=(Phi_t)_*xi(0)
                +integral_0^t (Phi_(t-s))_*v(s) ds.

It yields the transported tag and full moving-boundary spin of 0084.
Then apply the derivative of 0087's actual mean map
`g_mean=C_F composed with E[g]`, including the Moser correction and its
time dependence. This procedure is an explicit composed linear operator,
not the assignment `U_mean=U_canonical`. In a stationary observation chart
it may be included in `O(omega)` in (4); otherwise the time-dependent
version of (2) and the actual observation derivative are used.

For the *restricted compact pair* of 0085, its induced velocity is exactly
zero near the invariant tube boundary. Its reference tag equation has no
omitted normal source, and its actual tube and global spin are `B p`.
Eliminating its retained momentum gives `j=B²/P` and the displayed physical
spin, including the cross term before time-reversal pairing. This is a
genuine repair of the old exterior-velocity boundary defect, not an
unchanged geometric inertia. If additional freely evolving complement
directions are included, their velocity, tag and spin rows enter (2)–(4).
Compactness of the initially selected columns alone does not prove their
whole freely evolving complement remains compact.

## 5. An executed exact counterexample and its repaired action/current

The following positive four-dimensional Hamiltonian is an exposing oracle
for the block algebra, not an asserted Euler field:

    L=p qdot+r ydot
       -(3q²+2p²+7y²+5r²+2pr)/2.

Its Hessian is positive definite. Freezing the nonaffine pair `(y,r)=0`
gives the conditional action `p qdot-(3q²+2p²)/2`, with inertia `1/2`.
But the full equations have `ydot=p+5r`, so the zero-complement surface is
not dynamically invariant. This is the concrete residual (1), not a
warning based only on terminology.

Exact elimination of the actual `(y,r)` dynamics gives

    P_eff(omega)=2-7/(35-omega²),
    j_eff(omega)=1/P_eff=(35-omega²)/(63-2omega²),
    p(omega)=j_eff(omega) qdot(omega).

The SAME physical momentum is retained. The static relaxed inertia is
`5/9`, not `1/2`, and
`j_eff=5/9+omega²/567+O(omega^4)`. The exact full characteristic equation
is `omega^4-41omega²+189=0`; the frozen single-pair frequency squared is
instead `6`. Retaining the complement repairs both the equations and the
physical momentum map without fitting any coefficient.

Even with zero initial complement, for initial retained momentum `p0`
the full `q'''(0)` differs from the frozen model by `-7p0`. Initializing
the omitted state to zero is therefore not a proof it stays irrelevant.
The exact rational calculations are in `verify.py` and its saved first
execution; they are not a solver or a numerical stability assertion.

## 6. Quantitative validity and the minimum conditional closure

There are two honest uses of the original Cauchy–Born premise.

* For an exact conditional **constitutive action**, specify the geometric
  embedding `E`, retain its full KKS/Hessian and physical maps, and exclude
  nonaffine variation/relaxation. This is a concrete, independently testable
  material/profile restriction. The positivity and coefficient results of
  0085/0091 remain meaningful under that restriction. They need not prove
  nonlinear invariance of every microscopic Euler trajectory.
* For an assertion about **unrestricted Euler response**, retain (2)–(4),
  or establish a stated error bound for dropping them. A statistical closure
  can be phrased as the measured omitted contribution
  `E[A_ER r | retained observations]=0`, together with the actual observation
  correction from (4). Such a condition concerns nonaffine correlations;
  its satisfaction is not established just by stationarity of the base field
  or isotropy of its law. It cannot be replaced by assuming the desired
  Cosserat differential equation itself.

A direct finite-time bound makes the comparison falsifiable without
asserting a spectral gap. On the actual smooth-Euler Sobolev domain take
`||exp(t A)||<=exp(C_s t)` and define the finite-column residual norm
`C_R=||A E-E A_E||`. For full and restricted evolutions with the same
compatible initial state,

    ||eta_full(t)-E z_CB(t)||
       <=C_R integral_0^t exp(C_s(t-s)) |z_CB(s)| ds.         (5)

The constants follow from the background derivative bounds and the actual
profile columns; there is no empirical tolerance hidden in (5). Continuous
physical observations multiply this bound by their operator norm. Point
core-jet observations use a sufficiently high Sobolev index so that the
required derivative evaluation is continuous. This is a quantitative
short-time/error statement, not an assertion that its bound is small over
an optical period. Where a uniform complement resolvent bound is available,
the spatial derivatives of (3) can instead be computed by the already
importable noncommuting Schur-jet identities; absence of such a bound is
not repaired by a formal local frequency expansion.

## Result

`route_verdict: established` for the exact full-complement action/physical-
observation repair and the explicit restriction-versus-reduction discrepancy.
`evidence_scope: REPRESENTATION_SCOPED` for its application to the proposed
Euler coarse-graining. The original exclusion of nonaffine relaxation is
the pedigree of the conditional constitutive action; it is not a hidden
new assumption or a reason to delete that positive result. A stronger
free-Euler reconstruction or GLM-identification claim uses the computed
complement and observation formulas, rather than treating them as already
zero. The parent continues its positive compact-route construction with
this exact distinction and its independently measurable validity data.
