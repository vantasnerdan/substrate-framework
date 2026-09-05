# Independent review of the localizable-center compatibility theorem

Reviewer: fresh non-author `herdr geometry-review pane w3:p4`.  Review
boundary: `0252/center-compatibility.md`, `verify_center.py`, and the captured
eight-check receipt only.  I did not author the theorem or verifier, did not
edit 0248, and do not re-review the independent 2026 symmetry theorem, 0211,
or the parent periodic join.

## Verdict

**Established as stated.**  A `C2` axisymmetric steady Euler field near a
circle `r=r0>0`, with single-valued axisymmetric pressure, a poloidal zero
`v=(v_r,v_z)=0`, invertible poloidal derivative `Dv`, and the local first
integral `u dot grad p=0`, necessarily has zero swirl at that circle.  Hence
the whole circle is a stagnation circle.  This refutes exactly the subclass
of 0248 Candidate B that places its required nonzero elliptic closed core
inside a pressure-localizable matching region.

This conclusion is independent of Peralta-Salas--Slobodeanu 2026: it is a
local cylindrical Euler identity and does not consume analyticity, boundary
Bernoulli data, compact support, or the global symmetry conclusion of that
theorem.

`route_verdict: established`

`evidence_scope: EXACT_LOCAL_AXISYMMETRIC_CENTER_COMPATIBILITY`

## Independent derivation

Let `x0=(r0,z0)`, put `A=Dv(x0)`, and define

\[
 f=v\mathbin\cdot\nabla_{r,z}p.
\]

Axisymmetry removes the toroidal derivative from `u dot grad p`, so the
localizability hypothesis says `f=0` throughout a neighborhood.  At the
poloidal zero,

\[
 \nabla f(x_0)
 =A^T\nabla p(x_0)+D^2p(x_0)v(x_0)
 =A^T\nabla p(x_0).                                      \tag{1}
\]

Since `f` vanishes locally and `A` is invertible, (1) gives
`p_r(x0)=p_z(x0)=0`.  The physical radial component of steady Euler is

\[
 v_r\partial_rv_r+v_z\partial_zv_r-\frac{w^2}{r}=-p_r.   \tag{2}
\]

At `x0`, (2) reduces to `p_r=w^2/r0`.  The already proved `p_r=0`, the
strictly positive `r0`, and real velocity imply `w(x0)=0`.  No incompressible
streamfunction convention, sign choice, or Bernoulli gauge enters this
implication.

The companion first-integral argument is also correct.  Axisymmetric
toroidal Euler gives

\[
 v\mathbin\cdot\nabla(rw)=0.                              \tag{3}
\]

Steady Euler gives advection of `p+|u|^2/2`; subtracting the assumed pressure
first integral gives `v dot grad |u|^2=0`.  Applying the same differentiated
zero/invertible-`A` argument to both integrals yields

\[
 \nabla(rw)(x_0)=0,\qquad \nabla|u|^2(x_0)=0.             \tag{4}
\]

The first equality gives `w_r=-w/r0`, `w_z=0`.  Since `v=0`, the second then
has radial component `-2w^2/r0`, again forcing `w=0`.  This is a corroborating
derivation of the same local theorem; it is not counted as an independent
physical experiment.

The hypotheses are sharp in the direction tested by the attempt.  A locally
constant pure swirl `u=W e_theta` has `v=0`, pressure `p=W^2 log r`, and
`u dot grad p=0`, but `Dv=0`.  Thus deleting nondegeneracy would make the
claim false.  Conversely, the reviewed 0211 nonzero core does not satisfy
the local pressure-first-integral hypothesis near its nondegenerate poloidal
center; the new theorem does not contradict or narrow 0211.

## Eight-check receipt audit

The captured `first.stdout` contains eight runtime `PASS` lines followed by
`ALL 8 CHECKS PASS`; `first.exit` is exactly zero.  The verifier contains
exactly eight lexical `checks.check` calls.  It uses symbolic first jets,
which are sufficient for this pointwise theorem, and checks:

1. the differentiated localizability row `A.T * grad p`;
2. the algebraic inverse supplied by `det A != 0`;
3. the cylindrical radial Euler term with `-w^2/r` retained;
4. its zero-swirl consequence;
5. the angular-momentum first-integral gradient;
6. the speed first-integral gradient;
7. the false conclusion obtained if the centripetal term is deleted; and
8. the exact degenerate pure-swirl counterexample.

The second check does not by itself encode that the differentiated first
integral is zero; that premise is supplied by the theorem statement and the
first line of the analytic proof.  This is not a load-bearing gap because the
direct derivation above closes the implication independently of the script.
Likewise the first-integral checks corroborate the local jet algebra rather
than replacing the Euler proof.  No numerical tolerance, sampled rank,
solver status, or soft sign occurs.

Reviewed artifact hashes:

| Artifact | SHA-256 |
| --- | --- |
| `README.md` | `c51a7a0d5110968c79dcedc5dfb5c47004259d2635600556c128aaa6a1675c69` |
| `center-compatibility.md` | `3d6583a9037bf94d2564adbecd2223357cb12af4bd608266e1b5cfdc73708d5a` |
| `verify_center.py` | `8daec26befe468fd73b826af50276416e50f9bdbca478fc837f891191803e0f1` |
| `first.stdout` | `21d9a0d6c3252471136a2ac35951cb75dea854b1142f9d21d5ec4c781ba699dd` |
| `first.exit` | `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa` |

## Exact license boundary

`requires`: a neighborhood of a circle at `r0>0`; a `C2` axisymmetric steady
Euler velocity and compatible single-valued axisymmetric pressure; vanishing
poloidal velocity at the circle; invertible `Dv`; and `u dot grad p=0` on the
neighborhood.

`pass_licenses`: within that class, the core swirl is zero and a nonzero
closed elliptic core cannot occupy the pressure-localizable region.  The
core-localizable matching target of 0248B is therefore refuted by a named
mechanism.

`does_not_license`: a no-go for a nonlocalizable inner core with localization
confined to an outer transition; a no-go for degenerate circular streamlines;
a no-go for general compact steady Euler fields; any change to the reviewed
0211 ring; existence of a compact stationary array; or issue200 exhaustion.
It also does not promote the 2026 global symmetry theorem from compatible
guidance into this proof.

`maximum_verdict`: the local center incompatibility above and its route-level
application to the pressure-localizable-core version of Candidate B.

`failure_scope`: only a proposed matching construction whose
pressure-localizable region contains the nondegenerate nonzero poloidal core.

`unlocks`: representation change that preserves a genuine nonzero compact
Euclidean core without requiring pressure localizability at that core.

## Strongest failure-derived continuation

The strongest continuation is the **fixed exact constant-curl compact-core
background followed by the 0250 response construction**.  The 0145/0147
route already gives an exact stationary same-curl field with full pressure
`p=-|u|^2/2`, a contractible compact Euclidean nonzero invariant torus and
periodic realizations of positive fixed density; 0147 additionally supplies
an actual finite-action packet on that same torus.  Using one such fixed
background removes both newly exposed core-localizability conflict and the
need to solve a new stationary matching equation.  Its remaining
load-bearing achievement is precisely 0250's fixed-background full Euler/Lin
source-to-observation lift: pressure feedback, acoustic and optical history
range, physical action/current rows, and finite source costs at one fixed
radius.  Baldi's action-angle chart is useful transport-frequency input for
that analysis but is not itself the field-changing lift.

The common-frame 0248A1 construction remains the strongest route if retaining
the *literal 0211 generalized-Beltrami ring* is required.  It preserves the
nonzero core and has the correct small periodic image-shape source after the
uniform frame is restored, but it still needs the three-dimensional bordered
zero-frequency inverse, periodic Bernoulli periods, and then the curved
response/current transfer.

A nonlocalizable Candidate-B reformulation is logically live but currently
less closed: it must construct the inner-to-outer steady matching without
assuming `u dot grad p=0` at the core, prove whatever outer localization makes
the full velocity and pressure compact, and only then recover the dynamical
transfer.  The present theorem gives no obstruction to that construction,
but neither the 0169 localization identity nor the 2026 symmetry theorem
supplies it.

No parent conclusion or exhaustion follows.  The minimum reconciliation is
to mark only the core-localizable version of 0248B refuted and keep the three
continuations above distinct.
