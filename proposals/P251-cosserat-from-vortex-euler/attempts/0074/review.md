# Independent review: gradient, angular-response and isotropic-law join

Reviewer: `/root/smooth_core_review`, distinct from the construction
authors. Date: 2026-09-05. One frozen transaction covering 0065, 0069,
0070 and 0071. Source applicability is reused from 0060. No canonical
or authored proof file is edited by this review; subsequent 0072/0073/0075
constructions are not independently reviewed here.

## Decisions

- **0065 gradient-dominance theorem: established.** Its new cages preserve
  the complete kinetic jets and dominate the bounded full Schur/potential
  corrections at a finite carrier, conditional on the specified fixed
  base physical/current map.
- **0069 curl plus double-divergence locality theorem: established.**
  Retaining antisymmetric first force moments is compatible with C² full
  stationary action jets; no statistical decorrelation is required.
- **0071 finite Gaussian isotropization: established.** It gives a genuine
  fixed-helicity isotropic translation-ergodic Gaussian law, with the
  required local support and normalization, not a mixture masquerading
  as such a law.
- **0070 local angular-response construction: established.** Exact
  rotational KKS moments, unchanged physical core jets, positive local/full
  compact reaction energy and ambient-inclusive physical spin are supported.
  **Its claimed nonzero leading total-mean centroid transfer does not
  follow and is refuted for its stated mean coordinate by the omitted
  mean-energy cross below.** The local KKS theorem is not refuted by
  that correction.

The single load-bearing correction is to the passage from 0070's local
angular pairing to its physical diagonal mass and nonzero transfer claim.
The same retained first force moment necessarily induces a coherent polar
mean velocity at first slow order. Accounting for this in the already
declared full material mean split removes the added kinetic connection.
No further all-k or unrestricted invariant-manifold requirement is used.

## Evidence

I read the frozen proof/README and verifier for each attempt and their
saved receipts: 0065 19 checks, 0069 7 checks, 0070 20 checks, and the
0071 repaired 7-check receipt plus its original rational-domain failure.
The exact analytic identities and operator bounds are the strongest
oracles. These receipts were inspected and reused without duplicate runs.
The counts do not establish the missing physical mean map in 0070.

## 0069: retaining angular moment without an infrared singularity

The eight canceled moments leave exactly a scalar plus antisymmetric
first force moment. Integrating y_j partial_i f and y_j(curl V)_i
gives the stated signs. Choosing f=-c chi and V=a chi therefore removes
those moments from the residual G while retaining them in an explicit
compact gradient/curl decomposition of the same force. The former
averaged-center Taylor construction supplies its smooth compact double-
divergence remainder without changing the microscopic displacement.

For the complete projected velocity, Leray fixes curl V. Consequently
the mixed kinetic block is C*D, not C*PD with an uncanceled singularity;
it is a degree-three polynomial. The curl block is a degree-two polynomial.
The remaining D*PD block is homogeneous of degree four, smooth on the
unit sphere, and C³ at zero. Its derivatives through second order have
the stated polynomial bounds. Combining those symbols with the uniform
Sobolev-bounded compact primitive maps gives C² quadratic forms and
operator jets on the stationary reaction space. The zero-frequency atom
causes no singular denominator in this representation.

This proves the local analytic license, including cross terms and tails.
It does not assert that the slow Bloch-envelope velocity has zero
coherent polar mean merely because each unmodulated compact force has
zero integral. That separate first-order mean is precisely the issue in
0070 below; it does not invalidate differentiability of the full symbol.

## 0065: complete fixed-base gradient dominance

The hierarchy freezes the good-patch event, geometry, base carrier,
reaction operator and observable map before increasing the new gradient
carrier. This supplies genuinely carrier-independent constants. At each
response support the old and new corrected generators are scalar multiples
of the same single response field; new raw cages have disjoint supports.
Their mutual KKS crosses therefore vanish pointwise. The new eleven
affine moments also vanish, and the physical core jets are unchanged.

It follows that the patchwise reaction P(k) and D(k), and hence their
complete J(k), are unchanged by this lift. A retained angular moment in
the old reaction column does not spoil that support argument. Affine
symplectic crosses with the lift first enter at the omitted third slow
order. The centered bond symbol has the stated first derivative and
zero second derivative, so no false second-order amplitude is introduced.

Moving both Leray and curl onto a fixed base field bounds its energy
cross with the new field. The added force has uniformly bounded L²
norm, and its slow derivative differentiates a bounded local-coordinate
phase, not the high-frequency carrier. The base projected velocity needs
only C¹ here; both the old double-divergence representation and the
0069 extension provide that bound.

The full kinetic Gram remains nonnegative; the growing positive local
helicity diagonal has only bounded compact-return and shared-response
remainders. Three orthonormal bonds give their exact unit directional
frame sum. Rotating the whole geometry then averages the squared angle
axis projection to |q|²/3, yielding the stated common lower slope for
all longitudinal and transverse polarizations.

In the complete Schur difference, P and N0 are unchanged. The seven
displayed products are exactly the changes from L1,L2 and their adjoints;
the two L1-square/cross effects are not lost by assuming commutation.
R2 remains in the base coefficient and cancels only from this difference.
The stated norm polynomial bounds all of these terms. A fixed normalizing
map adds only the bounded stated first-jet cross correction. Thus the
finite-carrier lower bound is valid for the complete normalized gradient
coefficients, provided R_base is the actual corrected physical/current
base coefficient. An unfinished mean-map construction cannot be treated
as a previously evaluated R_base, but it does not refute this conditional
dominance theorem. New gradient cages have all first force moments
removed, so they do not introduce the first-order mean responsible for
0070's failure below.

## 0070: positive local rotational construction and physical spin

Adding B times the three rotational dual responses to S gives exactly
Omega(R_i,S)=B n_i, while translations and STF moments stay zero.
The dual span is isotropic and off the raw supports, so Omega(Q,S)=B
remains exact. Neither the physical Q jet nor the zero S angle jet is
changed. The equality of rotational and angle moments is a declared
reaction-geometry choice made from computed KKS data; it is not a
uniqueness theorem for all possible fluid reactions.

The new fixed response coefficients are uniformly bounded, and the
old positive carrier diagonal still dominates their complete finite
energy corrections. The full stationary reaction operator remains
bounded and coercive above a finite threshold, with no isolated-cell
kinetic factorization. These are valid positive microscopic conclusions.

The force first moment has the correct factor 1/(2 rho) and sign.
For a ball centered on the patch and containing the compact vorticity
change, Stokes gives integral_boundary n cross delta u=integral delta omega=0.
The radius-squared angular boundary term therefore vanishes on that
ball. The force's compact curl plus moment-free remainder also gives
an O(r^-5) velocity tail, so its physical angular-momentum integral
converges absolutely. It equals the rotational impulse B n for the
reaction direction. This is genuine ambient-inclusive fluid spin,
not merely the name of a canonical coordinate. Arbitrary material
faces still require the previously derived surface-current terms.

The local affine pullback correspondingly contains
`B s(qdot+n.betadot)=B s n.Phidot`. That local KKS identity is correct.
Its extension to a diagonal physical mean/spin kinetic matrix is the
step that needs correction.

## 0070 correction: induced mean velocity and its unavoidable cross

Write the coarse physical reaction-spin density as S_spin, to distinguish
it from the microscopic generator S. From the computed patch moments,

```
integral F=0,
M_ij=integral y_j F_i
              =c delta_ij-epsilon_ijm (S_patch)_m/(2 rho).
```

A slowly modulated patch has Fourier-envelope mean
`integral exp(-i k.y)F=-i k_j M_ij+O(k²)`. Its scalar part is longitudinal
and is killed by Leray. The antisymmetric part is not killed. Summing
with the actual stationary/Palm normalization gives

```
v_s = curl(S_spin)/(2 rho) + O(grad²).
```

Thus the reaction carries a coherent polar velocity at exactly the
first slow order relevant to the claimed optical transfer. Zero force
integral of the unmodulated patch does not remove this dipole term.
This conclusion follows from 0070's own local data, not from a new
statistical assumption or an all-wavelength requirement.

The full mean Euler energy already used in 0057 then contains the cross
between the raw Galilean/mean velocity V and this reaction velocity.
With P0=rho V, the relevant mean Hamiltonian terms are

```
H_mean+reaction=P0²/(2 rho)+P0.v_s+H_s.
```

Here H_s already includes the reaction's mean kinetic contribution;
no second copy of it is added. Combining this with the actual local
symplectic term gives

```
L=P0.Udot+S_spin.Phidot
                   -P0²/(2 rho)-P0.v_s-H_s.
```

Varying P0 yields the physical total mean velocity
`Udot=P0/rho+v_s`. Eliminating it gives

```
L=rho |Udot|²/2+S_spin.Phidot-rho v_s.Udot
                                      -[H_s-rho |v_s|²/2].
```

At first slow order, the exact periodic/compact-test curl integration
identity is
`rho <v_s,Udot>=<S_spin,curl Udot/2>=<S_spin,betadot>`.
The leading spin symplectic term is consequently S_spin.qdot, not
S_spin.Phidot in a kinetically independent physical total-mean sector.
The subtraction rho |v_s|²/2 modifies second-gradient reaction energy;
it does not restore the claimed first-order absolute-angle connection.

After the same independent reaction elimination, the resulting leading
physical kinetic cross is again b=-j/2 rather than b=0. With local locking
kappa |Phi-beta|²/2, g=-kappa/2 and therefore
`g-kappa b/j=0` at this order. This refutes the asserted nonzero leading
total-mean transfer for the frozen action. It does not rule out higher
gradient effects, another physical displacement observable, or another
microscopic construction.

The last physical-transfer checks in 0070 start from the already assumed
relative mass matrix and set bphys=j/2-j/2. They do not calculate the
mean response or P0.v_s. Their green result therefore verifies the
conditional matrix algebra but does not expose this omitted Euler term.
The minimum correction is to retain the local moment/KKS/positivity
theorems and retract the physical diagonal-mass/nonzero-total-mean-transfer
inference until the full mean/current map is included. This is the only
load-bearing correction requested in the frozen transaction.

## 0071: actual isotropic ergodic Gaussian law

The independent twelve-field sum is legitimate: all summands solve the
same curl equation, and their sum is again Beltrami and hence a smooth
stationary Euler field. Independence, not merely a random global choice,
makes it Gaussian and averages its covariances. The 1/sqrt(12) factor
retains point covariance I and the stated energy density.

On the unit sphere the source polarization spans the positive-helicity
line. The explicit outer-product identity gives its rank-one projector
times an even scalar weight of degree at most four in the distinguished
axis coordinate. The exact icosahedral second and fourth moments make
the averaged weight constant. Rotation covariance of the helicity
projector then makes the full spectral covariance isotropic; its trace
normalization fixes the coefficient to 3/(4 pi). A centered Gaussian
law is determined by this covariance. The construction therefore
overcomes the specific mixture issue identified in 0060.

Its smooth atomless spectral density gives covariance decay and Gaussian
mixing, hence translation ergodicity. Full local Beltrami support also
follows directly: one independent summand can approximate the scaled
target and the others zero with positive probability. The displayed
tolerances are more than sufficient by the triangle inequality. This
argument establishes applicability of the new law; it does not pretend
the old source named this finite superposition.

Independent homogeneous Poisson candidates with independent Haar marks
give the stated isolated-good-patch intensity. Isotropy makes the good
probability independent of orientation, and the isolation event is
independent of the field event. The Gaussian/Poisson system is mixing;
translation-covariant thinning and reconstruction preserve ergodicity.
The selected structures remain patches of one actual Euler field.
No kinetic diagonalization follows from this spatial ergodicity.

The theorem is for each fixed helicity. Reflection gives the opposite
curl law and still needs its stated componentwise ensemble handling.
Time reversal leaves the centered Gaussian law invariant but reverses
the KKS signs for a fixed physical realization. Neither this isotropic
spatial law nor its ergodicity implies the product-Haar local-state
closure or absent angular current of the separate contrast theorem.

The first verifier failure was a rational polynomial coefficient sent
to an integer Groebner domain. Using QQ repairs that implementation
without changing any theorem or tolerance. The seven repaired exact
checks are appropriate corroboration of the finite moment/projector
identities, not a simulation of ergodicity or knot existence.

## Frozen hashes and disposition

- `0065/full-gradient-dominance.md`: `16b403592ca9326db4551b826cff914ea509df8d2345356cbd9f8b7818207f42`
- `0065/verify.py`: `f874c40db01e430164617e7bb03699cbd85b13838ea2d671ac1f250e11a3db7a`
- `0069/README.md`: `cff8adbce4dba6bfe9b917b2a963958b001c9585ce33068a2b411e2d495b67ed`
- `0069/verify.py`: `0870a21983ce4c966090ff333ca7fa16d2b893ccdf31c9a0b72c36c5f8f3fda4`
- `0070/affine-angular-reaction.md`: `074c2c45cc3b6325fb77ff66a30636e60dc6f488b58967faa99637d4fec69810`
- `0070/verify.py`: `9b7412d8a7263d6deb50ca1ffd7a0df25df7d6ebe68faad7e7cd92d0d5043a7c`
- `0071/README.md`: `cc3fd57a7d9d2ae22f0ac159cc36dc3660ccddd89b29ad0367fa88221f236f96`
- `0071/verify.py`: `49b156a35a740555ee9afc7de937c41c8c5d75b3cb5694f534c136c0ee7f3fc4`

The positive local, analytic and statistical theorems are recommended
for acceptance at their stated scopes. The physical-transfer correction
is recorded once, with its explicit mechanism and minimum repair.
No full parent completion or later replacement verdict is inferred.
