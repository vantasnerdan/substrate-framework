# Independent review: positive correlated column response

Reviewer: `/root/smooth_core_review`. Reviewed 0187's two proof files,
verifier and original output, and 0181's actual mode proof, verifier,
original failures, diagnosis and repaired output. The relevant radial
Euler/KKS inputs in 0135/0137 and common-input conventions in 0172/0174
were checked at their use sites. This is the one substantive pass at
0190's registered boundary, not a review of the parent EPS continuum.

Independence: I did not author or implement 0187 or 0181. Earlier I
authored 0174's common-vector interface and explained its distinction
between averaging actual initial histories and imposing common
instantaneous histories. I supplied neither 0187's signed amplitudes nor
its energy-matching construction. The transverse-frame correction below
is new work derived during this review and reported to the author.

## Outcome and minimum repair

The positive scalar response/energy construction in 0187 is established.
Its actual smooth-column supplier in 0181 supplies the needed nonzero
gap, nonzero unequal slopes after scaling, positive action, and fixed-tag
physical current through carrier two. Signed preparation amplitudes do
not entail signed probabilities or negative Euler phase energy.

There is one concrete error in their joining: the quoted Haar factors
`c_T=1/5,c_L=3/5` describe an AXIAL angle whose marked axis is parallel
to the carrier. The 0181 mode is a TRANSVERSE covariance tilt. Its
scalar angle axis is perpendicular to its carrier. The correct factors
for the explicitly paired scalar preparation are

    c_T=2/5, c_L=1/5.

The positive construction survives: both factors are positive, and the
same amplitude slopes enforce both observed closure and physical energy
equality in both channels. This is a coefficient/physical-frame repair,
not a new condition on a mode or an excuse to discard the result. The
actual conjugate tilt and raw-to-whole-law phase/current factors should
be recorded explicitly as in section 4 below. Parent has accepted this
repair for the new append-only 0191; frozen 0187 remains provenance.

Thus the axial-supplier mapping as printed is refuted by its exact tensor;
the corrected transverse-supplier construction is positively derived
below. One check of the author's correction is the remaining review step.

## 1. The energy check is independent of observed frequency

Put `h_r=j_r A_r^2`, with all derivatives taken at the actual fixed
reference carrier and with each tag fixed. Direct differentiation of
the two measured rows gives

    E[(A h(nu))'']=E[e]h + B h_nu + V h_nunu,
    B=E[a b+2d v], V=E[a v^2].

Here `E` uses the positive family probabilities, not the action masses.
With `E[a]=1,E[e]=0,V=0`, both rows equal those of the oscillator with
frequency `nu+c_s |K|^2 B/2` through spatial order two, for every fixed
time interval. Their Wronskian is one at this order, not just initially.
The nonzero initial phase is instead

    J(K)=J0+c_s |K|^2 E[h'']/2, J0=E[j a^2]>0.

These two averages really are different. The square in the second is
necessary and retained in the proof/code.

For independent scalar initial inputs `(q0,r0)`, the conserved physical
Euler mode energy is `E[h(r0^2+nu(p)^2 q0^2)]/2`. Its excess stiffness
over `nu^2 J(K)` has second coefficient

    c_s E[h(v^2+nu b)+2nu h'v].

The observed stationary action instead has excess coefficient
`c_s J0 nu B`. Their equality is precisely 0187's equation (7). Thus
0187 is not mistaking a frequency average or a moving-embedding
Hamiltonian for physical Euler energy. Equality at the initial time,
conservation of the actual microscopic energy, and the already derived
observed oscillator rows extend equality to the fixed-time second jet.
No exact finite-K invariant microscopic two-plane follows or is needed
for this stated Taylor-jet result.

With the amplitudes fixed by normalization and variance cancellation,
the equations for a chosen finite `B_*>0` and energy equality have
coefficient rows

    (2w1 v1, 2w2 v2),
    (2nu w1 v1(2j1 a1-J0), 2nu w2 v2(2j2 a2-J0)).

Their determinant is `8nu w1w2v1v2(j2a2-j1a1)`, nonzero since the two
amplitudes have opposite signs and the masses are positive. This
independent differentiation agrees with the verifier's determinant and
the rational example. The example is explicitly not an Euler spectrum.
The selected `B_*` is a declared preparation parameter, not a uniquely
predicted intrinsic modulus; the artifact says so.

## 2. The supplier supplies actual Euler jets and actual material rows

The nodal branch argument is not a numerical-root assertion. The two
Sturm Wronskian formulas make `D_l<0`, give one simple root between
consecutive Bessel poles, and give a transverse laboratory zero at
`l=sqrt(3)x`. On its selected side `omega<0`, `nu=-omega>0` and
`nu_x>0`. The squared-frequency curvature is positive sufficiently
near that fixed crossing. A finite choice precedes smoothing.

The radial `(f,P)` system retains exterior pressure and contains no
derivative of the thin vorticity taper. Since `omega-O_e(r)<0`, its
annular parameter coefficients stay bounded; its transfer and the
decaying exterior matching determinant converge with the required
parameter derivatives. The simple-root IFT therefore supplies an
actual smooth Euler mode and its carrier derivatives in the declared
active-vorticity class. It does not require a gap to arbitrary exterior
vorticity states. That stronger spectral statement is not used here.

The full KKS integral has positive core numerator
`2tau integral(r P'^2+P^2/r)+(tau^2+4)P(a)^2/2`; its bounded annular
correction preserves the sign. With the displayed laboratory generator,
the Hessian is `-beta omega I`, positive on the selected branch. The
physical scalar mass `-beta/(omega c_theta^2)` follows by eliminating
the conjugate coordinate; it is not the tag's polar moment.

The Cartesian material calculation also retains both terms in
`S=integral rho chi (r cross xi_t+2xi cross u0)`. In particular, the
axial-displacement integral `R2` enters both shape and spin. Setting its
first three carrier jets to zero makes the actual displacement/current
rows `G=Delta theta,S=Delta theta_t` through order two. These are not
inferred from a canonical momentum name. The two-annulus positive IFT
is valid: at consecutive pressure zeros its first two center columns
contain independent vectors `(1,r1^2),(1,r2^2)`, and the positive-ratio
column controls the remaining derivative row. The tilt row survives
because `1/r1^2-1/r2^2` is nonzero. Smoothing the bumps and the vortex
therefore preserves positivity and the fixed-tag property.

Euler amplitude scaling fixes a common positive gap without specifying a
numerical Bessel root. Geometric scaling with distinct radii then gives
actual nonzero unequal `v_r=Omega a_r tau'(x0)`. The resulting mass
derivatives remain the actual ones in the linear system, not free
constants. Reducing each tag fraction to a common raw `Delta` leaves
its normalized covariance angle, and consequently its phase mass in
that angle coordinate, unchanged; it scales its literal moments only.

## 3. Scope of the common wavevector

Both real carrier bands and the entire stationary column, tag, and phase
must be rotated with ONE common laboratory `K`. The parameter shift is
along the COLUMN axis `t`, namely `p +/- t.K`. This remains distinct
from its transverse scalar tilt axis `n`. Phase compensation is a
representation of the actual shifted full Euler history, not an
omission of the Leray tail. On this supplier its interpretation is the
axial Fourier/Bloch-fiber and phase-averaged action density, with full
radial pressure decay. It is not a finite-total-energy R3 packet or an
already localized periodic-cell construction. 0187 states that boundary.

The second jet must be formed after real-band and whole-law pairing.
Odd columns are not set to zero on an individual circular mode. The
explicit time-reversed pair below cancels the unwanted transverse
conjugate observation; the same real-band conjugation as 0174 then
gives the even carrier jet. No rotating observer clock is introduced.

## 4. Constructive correction: transverse mode and whole-law factors

Take a registered right-handed frame `(n,b,t)`, with `b=t cross n`.
For the equal-probability pair of actual time-reversed column modes
write their full material tilt as

    Theta_s=n theta+s b theta_t/nu,  s=+1,-1,
    theta=A(p)[cos(nu t_time) q0+sin(nu t_time) r0/nu].

Indeed `Theta_s,t=-s nu J_t Theta_s`; this is the actual circular
laboratory evolution supplied by the two reversed backgrounds, not
two independently imposed components. The companion tilt/rate is
retained in each mode, including its parameter derivatives. It cancels
exactly between the two whole realizations. Each scalar phase/energy
has the positive raw mass `M_raw` of 0181; averaging the pair does not
double this mass or negate its action. Their literal rows remain
`G_s=Delta_raw Theta_s,S_s=Delta_raw Theta_s,t` to the supplied order.

For common vector initial data use `q0=n.q_macro,r0=n.r_macro`, and Haar
rotate the COMPLETE frame and preparations. Measure the registered
scalar angle by its actual projection `theta=n.Theta_s`, and reconstruct
`q_macro(t)=3 E[n theta]`. Then

    E[n_i n_j]=delta_ij/3,
    j=M_raw/3, Delta=Delta_raw/3.

Thus the averaged action uses `j`, while literal full physical current
uses `Delta`; neither raw scalar coefficient can simply be substituted
for its averaged density. The rigid-rotation response of the projected
tag angle is one. The signed optical preparation is not asserted to
be the rigid rotation of every tag simultaneously.

To derive the next tensor without guessing, first condition on `t`:
`E[n_i n_j | t]=(delta_ij-t_i t_j)/2`. Then

    E[n_i n_j(t.K)^2]=(2|K|^2 delta_ij-K_iK_j)/15,
    3 E[nn^T(t.K)^2]=(2|K|^2 I-KK^T)/5.

This yields `c_T=2/5,c_L=1/5` for BOTH the observed rows and the action
and physical current jets. It disproves the imported axial tensor but
preserves the full scalar cancellation/energy theorem. In particular
the two families still have a common geometric tensor, so the SAME
two derived slopes work in both channels, leaving positive curvature
`c_s nu B_*`, positive initial mass, and nonzero literal `G=Delta q`,
`S=Delta q_t`. A downstream pointwise curvature-energy representative
must use this corrected ratio, not the optional axial-ratio choice in
0172; no parent constitutive claim is made by this review.

## Evidence and validation boundary

The strongest oracle is the exact differentiation, radial Sturm/IFT
proof and actual Cartesian moment calculation above. Existing targeted
receipts were reused after reading their code and original outputs:
0187 has 10 exact passing checks; 0181 has 19 after its two explicitly
preserved implementation repairs. Neither receipt tests the disputed
tilt-versus-carrier Haar tensor; the tensor calculation above is the
exposing independent check. No numerical spectral design, full suite,
or canonical validation was added for this review-only delta.

Verified SHA256 inputs:

| Artifact | SHA256 |
| --- | --- |
| 0187/correlated-family-action.md | a2376bf497535fe100f4ed68bc9da2ccb8d28b7951bf037d864cbfbb98d538f2 |
| 0187/nodal-family-realization.md | 13627b4eb81d32565a741839dd1e6d051bc648e239d3e4b0457e85fd1d6528e5 |
| 0187/verify.py | 48f0c77a5eaf86461e9df621e1576866819776c159cfee6ae65ae0e3c2c2bba9 |
| 0187/first.stdout | f3fb1fb8272dacc257a0a4f1610b8a8c06fa3553ab535bcee895a356a1c15fdf |
| 0181/nodal-domain-mode.md | 144e3b980ab51e1b24917e0b10f69db3127b95599f81db8a7c225e1cd8151838 |
| 0181/verify.py | accd68d04c0d97594c32437ac44541a8119070c92c004e770add16e10ab7e5a3 |
| 0181/repaired-run.stdout | 9d158c4129216b247796a7a6e68d6f36b71e2c4f2e795d4594a68f0a29791387 |

The parent objective remains active. This positive column-mixture result
does not claim one same-cell EPS field, finite-action packet transfer,
uniform acoustic-time homogenization, or closed translational memory.
Those are independent constructions, not additional rejection criteria
for this correctly scoped response theorem.
