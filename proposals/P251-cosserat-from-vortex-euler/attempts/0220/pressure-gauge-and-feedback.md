# Full pressure gauge and the curved fast reaction

All source cross-sections and collars in this argument are fixed
before R increases. The actual toroidal harmonic n is fixed and
nonzero, k=n/R. Estimates concern the compact active-vorticity
orbit; arbitrary independently supplied exterior vorticity is not
included in this spectral subspace.

## 1. Pressure-gauge normalization of a compact force

Let F be a compact covector source in the comparison tube and let
v=S_bar F be its complete solenoidal velocity response. On the source
and a fixed surrounding collar,

    F=g_bar v+grad_k p.

Choose chi=1 on the entire source and zero outside that collar. Let
C be the transverse average of p over a fixed annulus lying OUTSIDE
the source, and replace the source by

    F_sharp=F−grad_k[chi(p−C)].                           (1)

This does not change either complete Leray response: the subtracted
potential is a legitimate compact smooth potential on both the
periodic comparison tube and the actual curved tube. In the source,
F_sharp=g_bar v+ikC e_y. In the annulus, F=0 and
grad_k p=−g_bar v. Transverse Poincare and the axial derivative give

    ||p−C||_annulus <= C_0 ||v||_annulus,
    |k C| <= C_0 ||v||_annulus.

Therefore

    ||F_sharp||_2 <= C_0 ||v||_(source+collar).           (2)

This remains uniform as k tends to zero. Subtracting chi p WITHOUT
removing the transverse mean would instead produce a false 1/k
bound on the cutoff return. The retained ikC term is the actual
axial current, not a reset of a harmonic mean.

Equation (2) removes the large gradient/gauge part of a coadjoint
force before comparing physical pressure operators. In particular
the k^-1 force estimate of 0213 is not the quantity that must be
multiplied by a first-curvature pressure error.

## 2. The actual exterior pressure comparison

In the raw toroidal chart the complete scalar pressure inverse has
the exact Fourier kernel

    G_(R,n)(x,z;x',z')
      =−R/(4pi) integral_0^(2pi)
       exp(−in alpha)/sqrt(d²+4(R+x)(R+x')sin²(alpha/2)) d alpha,
    d²=(x−x')²+(z−z')².                                  (3)

This is the full R³ Green function integrated over the actual ring,
not the streamfunction Green function and not a radial wall. Its
comparison kernel is −K0(kd)/(2pi). Source divergences and output
gradients include their actual axial ik factors.

On two fixed source/collar disks, split the integral at a fixed
small angular cutoff, and in the near part use t=R alpha. Subtract
the common logarithmic singularity and its first affine metric
coefficient. The remaining far part is smooth in both source points;
its transverse derivatives follow by differentiating a denominator
bounded away from zero. The resulting pressure-kernel difference
has: a scalar constant depending on n; an affine (x+x') term of
size R^-1 log R; local first-metric logarithmic singularities of
size R^-1; and a smooth twice-differentiated remainder of size
O(R^-2 log R). The local singular terms are equivalently the first
metric pressure variation P D P from (7) of `actual-operator.md`.

The constant and affine logarithmic pieces do not produce a large
velocity error. The transverse integral of div_N F_sharp is
ik integral N_yy F_sharp,y. A transverse output derivative of the
affine piece therefore costs another k. An axial output derivative
also supplies k. The constant piece contributes only O(k²) to the
velocity. These are source identities, not suppressed exterior modes.

For the singular local part one may extend the local metric
perturbation smoothly to a small positive metric perturbation on a
fixed comparison space. Its full projection difference is bounded
on L² by the exact inverse-metric identity and elliptic energy
coercivity. Subtracting that local model leaves only the smooth
harmonic-return kernel just bounded from (3). Thus no estimate of a
singular kernel by its nonintegrable absolute value is needed.

Pullback to the actual action-angle coordinates composes this
calculation with smooth fixed-collar maps that differ from their
limits by O(R^-1 log R). It yields, for compact covectors,

    ||(S_R−S_bar)F_sharp||_(core+collar)
       <= C R^-1 log^a R ||F_sharp||_2.                  (4)

Here a is a fixed finite power allowing the stated Green-coordinate
remainder; its value is immaterial to convergence. S_bar uses the
actual mean metric in the core, smoothly continued to a radial flat
metric in the exterior. Its exterior covariance G=0 and constant
poloidal circulation give a genuine irrotational continuation.
The regular axis and decaying Bessel return are retained. The metric
continuation is an auxiliary pressure comparison, not a replacement
of the actual background.

Combining (1)-(4), and using annihilation of the same compact gradient
on both sides, gives the stronger response-specific estimate

    ||(V_R−V_bar)Xi||_(core+collar)
       <= C R^-1 log^a R ||V_bar Xi||_(core+collar).      (5)

The full actual exterior pressure is included in (3)-(5). The
comparison has not assumed a globally small toroidal metric.

## 3. Same background and the physical-velocity contraction

The full and averaged core metrics have the SAME contravariant u,
omega=f u and volume form: averaging the actual covariant curl
identity commutes with differentiation. Their Euler/Kelvin brackets
therefore coincide on the active support. It follows there exactly
that

    (A_R−A_bar)Xi=(V_R−V_bar)Xi.                          (6)

For two globally defined positive metrics on the same volume space
the underlying algebra is the useful exact identity

    (S_g−S_bar)F=−S_g(g−g_bar)S_bar F.                   (7)

Equation (5) is its actual-ring, different-exterior counterpart,
proved through the compact pressure gauge and (3). Both formulations
act on the controlled PHYSICAL velocity, not on a large gauge-parallel
generator. This is the representation change that removes the
flat-edge amplification from the first curvature step.

## 4. Reference fast inverse and its domain

The comparison metric and all its background coefficients are
poloidally invariant. In m=0, `mean-metric-response.md` supplies
positive-sector coercivity and the complete response. Its source
dual norm can use L² eta: integrate d(E eta_y)' against Y by parts,
then use control of (dY)' including the tail Hardy estimate. Thus
the physical velocity response is bounded uniformly in small k,c.
The compact force graph has at worst k^-1 times fixed logarithmic
powers, as established by its exact generator reconstruction.

For m!=0, active-support transport denominators obey
|m Omega+kW−omega_mode|>=|m|Omega_min/2 at sufficiently large R.
The actual compact-source pressure coupling is factored AFTER that
transport inverse, as in 0201. The inverse curl gains one derivative;
the transport multiplier cancels its single angular numerator.
This gives a compact Fredholm family, with bounded velocity response
to an L² generator source. At k=0 the new 0211 radial ground-state
identity proves that only the two m=1 translations lie at zero.
Remove that finite symplectic phase, not the full m=1 channel.
The finite low-m complement has a bounded inverse by Fredholm
continuity; high |m| follows from its explicit transport denominator
and elliptic angular term. All constants are uniform at the fixed
smooth profile. Small averaged-metric changes preserve this bound.

Denote this reference translation-complement physical response by
T_bar. A source supported only in a collar where omega=0 represents
zero physical forcing; sources are identified accordingly. Smooth
solenoidal source germs admit a bounded fixed-collar completion.
In m=0 its radial flux is O(k), by the exact divergence identity,
so that completion does not create a 1/k cost.

Equations (5)-(6) give a Neumann factor of size
C R^-1 log^a R on this physical response space. Consequently the
ACTUAL operator, compressed to the reference translation complement,
has a full bounded physical inverse on the bending contour. Its
coadjoint-force bound retains the reference k^-1 polylog bound.
This step sums the actual curved metric/pressure feedback; it is not
a first-Born substitution.

## 5. Exact geometric bending complement

The reference translation projection is not the physical geometric
KKS projection, and n>=2 translations are not global symmetries.
Let b be the exact 0213 bending columns, with nondegenerate inherited
KKS matrix J_b, and M_b their restricted action generator. Their
actual residual is r=A_R b−b M_b=O(R^-2 log R) on the active source.
Let Q_0 be the reference complement, and write a vector as b x+y,
with y in Q_0. Its exact b-KKS-orthogonality condition is

    J_b x+Omega(b,y)=0.

The complementary equation is

    (z−Q_0 A_R Q_0)y=eta+H(z)x,
    H(z)=Q_0 r−Q_0 b(z−M_b).                             (8)

On the actual bending contour, H(z)=O(R^-2 log R). Substituting the
already constructed inverse from section 4 into the orthogonality
row gives a FINITE matrix correction bounded by

    C k^-1 polylog(1/c) R^-2 log R
       =O(R^-1 polylog R),                               (9)

after canceling the common per-length R in the KKS matrix and its
pairings. It is therefore invertible. This changes to the exact
geometric bending complement without pretending it was invariant.
It retains a k^-1 polylog coadjoint-force bound for the complete fast
inverse, sufficient for the action estimate.

The first-reaction calculation of 0213 now applies with that FULL
inverse, so the complete self-energy is O(R^-2 polylog R), lower
order than the positive O(log R/R) bending Hessian. No other fast
channel or exterior pressure has been omitted from this estimate.

## 6. Mode and observation meaning

The finite Schur determinant is a uniformly small analytic
perturbation of the inherited positive bending phase on its
positive-frequency contour. Argument counting preserves its simple
root in each fixed complex toroidal sector of the invariant compact-
active-vorticity space. The Cauchy derivative bound on a comparable
inner contour controls the self-energy derivative by a fixed
polylogarithmic power, whereas the bare KKS row has size R. The
derivative correction is therefore small as well. The exact
Hermitian Schur derivative is the inherited full-mode Krein form,
including the eliminated fast phase. Hence that form and the
quadratic phase energy remain positive; Hamiltonian conservation
keeps the root on the real frequency axis. Its leading scale is

    omega_n²=[Gamma log(R/a_ref)/(4pi R²)]²
                   n²(n²−1)[1+o(1)],   n>=2.            (10)

The old O(1/log R) finite-core action correction and the smaller
reaction correction are both included in o(1). At n=1 the exact
Euclidean Jordan row, not (10), is used.

The induced smooth decaying velocity is an actual solution of the
full Euler eigen-equation, including its exterior harmonic pressure.
The isolated-contour assertion is in the declared active-vorticity
sector, not a claim that arbitrary exterior passive-vorticity
transport has a spectral gap. Smooth interior material reconstruction
uses the noncritical actual torus transport denominators; it does
not require an L² material-label displacement at spatial infinity.

The physical mode has a nonzero complete-domain Euclidean n=2
quadrupole row on the chosen invariant inner torus. On that torus
the vorticity and axial speed are bounded away from zero, so the
force estimate controls the normal displacement. A large displacement
parallel to vorticity is a relabeling tangent and contributes no
uniform material-domain covariance variation. Full tag transport
and boundary terms are retained. The positive phase action and the
actual quadrupole pole are therefore on the SAME global Euler field
with its literal constant-curl inner tube.

This n=2 result has exactly zero linear global material-spin vector
by the unchanged tensor-rank selection. The circular reference
covariance has no pre-existing planar eigenvalue gap; the physical
quadrupole is not silently divided by that zero gap to manufacture
a linear director. A nonzero-spin/director or common macroscopic
continuum join remains a separate parent obligation.
