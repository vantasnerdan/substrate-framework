# 0081 — independent finite-parcel and moment-bridge review

Reviewer: `/root/smooth_core_review`, distinct from the 0073/0078 authors.
One frozen evidence transaction: 0073 universal material-moment identity and
0078 exact parcel reduction plus its finite-ball existence example. The
physics-erdos-loop contract governs this claim-level review. The oracle is
direct calculus, symplectic reduction and the explicit field, corroborated
by the saved exact-arithmetic receipts; no numerical soft-eigenvalue claim
or empirical comparator is involved. No author files were changed.

## Verdict and strongest supported result

Established as stated in their declared evidence scope. The material
centroid, spin, mechanical connection and first momentum multipole obey
the displayed exact identities. A smooth Beltrami fluid in a spherical slip
domain supplies a physical reduced angle and a positive **selected quadratic
fixed-Kelvin action** with inertia `I_red=B²/h`. Its inertia is computed from
the same Euler KKS form and energy Hessian, not appended geometric mass.

This is a constructive finite-domain example, not the EPS/ambient parcel
assembly. The sphere's material-shape constraint, its rotation symmetry and
normal pressure traction are actual hypotheses. No load-bearing correction
was found. No parent-completion or claim-promotion verdict follows here.

## Material centroid, connection and Kelvin distinction

Centering the material integral gives both cross cancellations:
`integral r dm=integral pi dm=0`. Substitution into the kinetic energy and
cotangent one-form gives precisely `P²/(2M)+||pi||²/2` and
`P.dX+integral pi.dr dm`. These statements concern the actual selected mass,
not a relabeled Eulerian point mean. Their canonical algebra is the material
cotangent algebra; admissibility and descent after the fluid/boundary
constraints retain the qualifications stated in section 3. In particular,
the fixed-ball example is not evidence that its centroid is freely movable.

The polar chart is locally rotation equivariant when its correlation matrix
has positive determinant and is invertible. Differentiating `r=Qh` yields
the full cross term `Omega.integral h cross hdot`, with its stated sign.
Completing the square or eliminating genuinely cyclic shape momenta gives
the two displayed Schur complements. Positivity requires the stated
independent positive block; neither operation licenses an arbitrary fluid
shape momentum to be treated as a Kelvin invariant. The six-sample example
correctly exposes a nonzero connection and the difference between locked
and reduced inertia; it is not offered as an Euler discretization.

The relabeling momentum is the pulled-back velocity one-form modulo exact
forms. Curl alone omits circulation periods on multiply connected domains.
Commuting spatial and relabeling actions preserve spatial momentum maps
when the domain and traction constraints preserve that spatial symmetry.
The caution about a label-template frame is valid: for example, at zero
Kelvin momentum, relabeling rotations of a ball change the template frame
without creating a physical rotor. It is not necessary to assert that all
rotations stabilize every nonzero fixed Kelvin datum. A physical vorticity
angle avoids that label-frame inference.

## Explicit ball, physical angle and positive Hessian direction

Choose a positive root `z_*` of `j_1` and `R=z_*/lambda`. The profile
`f=j_1(lambda r)/r` is an even analytic function of `r` at the origin, with
`f(0)=lambda/3`. Its radial ODE implies the stated vector Helmholtz identity;
therefore `u=V+curl V/lambda` is a smooth Beltrami field. The exact radial
velocity is `2fz/lambda`, so the root gives zero normal velocity. Normal
pressure on the sphere has zero centroid torque. These facts establish the
stationary Euler/slip example without a singular core or fitted constant.

Rotation about `e_x` changes the actual central vorticity, and differentiation
of `atan2(-omega_y,omega_z)` gives `dq(K)=1`. On the ball the rotated velocity
and `P_D(K cross omega)` have the same curl, divergence and normal boundary
data; uniqueness of the ball's div-curl problem identifies them. In
particular `curl v_K=lambda v_K` and `v_K` is not zero.

The compact off-core direction has the exact pairing

```
Omega(K,eta0) = -rho integral v_K . curl(chi curl v_K)
               = -rho integral chi |curl v_K|² < 0.
```

The integration by parts has no support boundary term. Its vorticity
variation vanishes near the observed core, so `dq(eta0)=0`, even though
its velocity reaction need not vanish there. This distinction is important
and is correctly used by the construction.

The finite-domain projector preserves a compact curl, annihilates gradients
and is an L² contraction. Thus the previously reviewed compact-carrier
parametrix applies here without replacing the ball projector by a free-space
one. The correction coefficient involving `Omega(K,A_k)` stays bounded:
pairing with the fixed smooth `v_K` is bounded by the carrier's L² norm.
The fixed-response Hessian cross is bounded by moving curl onto its fixed
smooth velocity. Consequently the positive order-`k/lambda` term survives
the exact pairing correction at a finite carrier. It supplies `h>0` while
keeping `B=B0` and the physical angle jet unchanged.

The Hessian's curl boundary form really vanishes for these directions. Their
boundary normal vorticity is zero, so their tangential velocity one-forms
are closed. On `S²` they are exact, and the surface integral of the wedge
of two exact one-forms vanishes. This includes the rotation direction;
the same assertion on a higher-genus boundary would require the additional
period information explicitly mentioned in the proof. The finite ball is
therefore a legitimate domain for the displayed symmetric Hessian.

Stationarity and exact rotation invariance give `H(K,eta)=0`. Restricting
the quadratic KKS/Hessian action to the selected two directions gives
`B p qdot-h p²/2`; eliminating `p` gives `B² qdot²/(2h)` and the actual
linearized physical momentum `delta J=(B²/h)qdot`. This is a quadratic
action statement, not a nonlinear invariant two-dimensional Euler
submanifold. Other retained rotations or shapes require their full KKS and
Hessian blocks, exactly as the proof stipulates; they cannot silently be
deleted in a downstream assembly.

## Universal moment bridge

Material differentiation gives `Idot_ij=Q_ij+Q_ji`, while the definition of
spin fixes `Q_ij-Q_ji=-epsilon_ijm S_m`. Thus
`Q=Idot/2-epsilon(S)/2`. Taylor expansion of the actual parcel momentum
distribution gives `-div Q`, hence the signs and factor one-half in
`-div(Idot)/2+curl(S)/2` are correct. Centering cancels the translation
contribution at this derivative order, not its higher multipoles.

An isotropic axial-rate response cannot map linearly to a symmetric
second-rank shape tensor: the only isotropic rank-three candidate is
antisymmetric in the two tensor indices. This removes that zeroth-order
response only. The retained Taylor remainder, actual full spin, ambient
parcel population and shape response are necessary inputs to any complete
physical observation map. 0073 states these conditions and does not equate
an isolated source impulse with the full material spin.

## Receipt audit and provenance

Read both verifiers and saved logs. The five 0073 checks corroborate the
tensor decomposition and Fourier sign. The seventeen 0078 checks corroborate
the connection, reduction, radial identities and physical-angle derivative.
The preserved first-run failure was an algebraic checker defect: a `z²`
substitution missed `z³`; polynomial remainder modulo the sphere constraint
is the correct repair and changes no field or theorem. The existence and
boundary proof above, not the check count, supports the infinite-dimensional
and finite-carrier assertions. No duplicate numerical oracle was needed.

SHA-256 of the frozen reviewed inputs (paths relative to `attempts/`):

```
d94f1f9eb23b1231e00e35ed6e03b4a703470318b5472b786af3da2a24e16d9c  0073/README.md
2a23448f6141f736421404e6bfb1b0ee825b1bb926a4b0f79e29cbe246ba98b2  0073/verify.py
78e35c57d2fd7c00eb88fe6c6fbcf2929d7500a416b20e78cad6d628f74a5d06  0073/first.stdout
d7b718385e5629f469ad92ad75bcfcc7bb9f56709ced1bbe462873ef70c55340  0078/README.md
ece68f5a41a86c9ebd8e1bc3f9f61e8acf2528f688e47c9c18f45a14e803c95b  0078/finite-parcel-reduction.md
ea988e2e8a66d3e2e00f2c45f34550e87de86cb5af8f9535354dbdb71851d863  0078/verify.py
28138fe75fe1cebc78f8ce9dd307687081ca106b3f8128b7affe52fc9c783a41  0078/first-run.txt
3c9086f5eaff695146ab5c8ef3c1c6ede09c99d26a7d89e073361722af028b63  0078/stdout.txt
```
