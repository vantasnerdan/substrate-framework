# P253/0026 independent operator-transfer review

## Frozen transaction and independence

This transaction adjudicates two fixed units separately. Unit A is attempt
`0019`'s smooth compact Gavrilov carrier, whole-space linearized Euler group,
one-meridional-circuit Kelvin/Leray return, and exact linear packet transfer.
Unit B is attempt `0021`'s full independent-label García--Hassainia--Hmidi
contour Jacobi IVP and reduced one-period monodromy. Neither verdict uses the
other unit as evidence. Attempts `0024`, `0025`, and the concurrent `0022`
work are excluded.

Central activation is recorded with exit zero. The frozen `0026` README hash
is `5aab2bae131e9d68ab1a34085d2965ebaf3b17a03529030624f4eeb635cc27f0`.
The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`. The reviewer authored or
implemented neither target.

The substantive-pass boundary was:

| Unit | Artifact | SHA-256 |
| --- | --- | --- |
| A | `0019/README.md` | `aa44d878cac71efbb5864e2d8626fbf9359f01965e379f834c7b42f190c815fe` |
| A | `0019/access-inventory.md` | `5fd8c34a4c2c51018fb16ca6f75a277b75e4cec20e1e0deaf2c9361158e0c7f4` |
| A | `0019/derivation.md` | `b7bc7926e75d1eaae1b681b509b7508063403757673b8d14a14c29914591bfba` |
| A | `0019/result.yaml` | `e6fcc43e552d48d07376cfed47e4cb65bfaef2c39c7cd91bbcaed78d3450d757` |
| A | `0019/validation.md` | `295237f75d31791fd1fdbabe85c14209729543e5cea2d8aed0538269bf91e79e` |
| A | `0019/verify_return.py` | `f47f7554f495956b1382da2754b55ee341482871e92e099aeb6fec51ad73c98d` |
| A | `0019/first-run.stdout.txt` | `f83812fc7fad60ea23b417eaff17c1630bfc8e8c0e0532eccc1789b4559c31b5` |
| B | `0021/README.md` | `e44c5ad0b74b8284db2ecedbbf252492f105ada4bbb5e9a816c9941a54ce3e4c` |
| B | `0021/construction.md` | `d1c27273e34903e76126233076fdcb87cbe41f47ffb37bb8c8700b7087c301aa` |
| B | `0021/source-receipt.yaml` | `d09dd90b2f3c70ca953464bcca583dc17d306d73405c9fc169b825d9da54e143` |
| B | `0021/validation.yaml` | `d0e9cf4cc2a4f50f2aa9e933670ba758a63f6e0dfc18f336021746dcedc3d37a` |
| B | `0021/verdicts.yaml` | `8facbc83c22e2fccb4b4782288ac12a807e0440409452a192d6d58b36f1d61a9` |
| B | `0021/verify_exact_structure.py` | `6686dd53972ff8b20ac289655f8a439fd3ea7a262115adefefab4b6335bec840` |
| B | `0021/exact-check.stdout` | `9f7228b296c6f9c4500e6a746741a945c50f96c7f356dd43e727c90f68f69c66` |

The activated directory inventory resolves the README's provisional filename
metadata: `0019` names and contains no `euler_core_packet` module or tests.
Its validation receipt expressly says that no source module was written. No
API or module theorem is therefore presumed or reviewed.

## Unit A: source carrier and whole-space operator

The pinned Gavrilov PDF hash is
`fcaca85faa77e3876b11d16718037169fc026112dfd3e4248e03e963a0ebc3c9`.
The primary construction supplies a local analytic pressure variable with a
strict nondegenerate minimum, the velocity

    u=(p_z e_r-p_r e_z+b e_theta)/r,
    p=R^4 a/4,  b=R^3 sqrt(H(a))/4,

and the localization rule `u_tilde=omega(p)u`,
`d p_tilde=omega(p)^2 dp`. Choosing the cutoff in a pressure shell avoids the
nonsmooth center and is flat at both interfaces. Extending by zero therefore
gives a nonzero `C_c^infinity(R^3)` divergence-free steady Euler velocity;
the pressure is constant on each zero-velocity component and can be normalized
to have compact support. The field is smooth on the symmetry axis because its
support is a positive distance away, and it has finite kinetic energy. The
source supplies no perturbation-stability theorem.

For divergence-free `v`, direct linearization gives

    v_t=-P[(u_* dot grad)v+(v dot grad)u_*],                (A1)

where on `L^2(R^3)` the whole-space projector has symbol
`I-k tensor k/|k|^2` for `k != 0`; its value at the measure-zero Fourier origin
is immaterial. Equivalently `P=I-grad Delta^{-1} div` with the usual negative
Fourier symbol for `Delta^{-1}`. Thus (A1) retains the Newton/Fourier pressure
and does not preserve compact support in general.

The projected transport part is skew on the divergence-free smooth core. Its
closure generates the isometric projected-transport group: this can be seen
without assuming the conclusion by writing its pressure correction as an
order-zero Riesz operator on `(grad u_*)v`, solving the smooth transport
regularizations, and using the conserved `L^2` energy in both time directions.
The remaining map `-P(v dot grad u_*)` is bounded. Consequently the full
operator generates a strongly continuous group and

    (1/2) d||v||_2^2/dt
       =-integral v dot sym(grad u_*) v,
    ||S_*(t)|| <= exp(|t| ||sym grad u_*||_infinity).        (A2)

The maximal domain need not equal `H^1` in the flat exterior; `0019` correctly
uses the generator domain and only calls `H^1` a core. Equivariance transports
smooth coadjoint-orbit tangents, so their kinetic-norm closure is invariant and
the stated graph-domain part is legitimate. The density factor `rho_0` and
the cylindrical factor `2 pi` occur in the exact kinetic norm and introduce no
missing sign or density term.

## Unit A: annulus, return, and packet transfer

On the cutoff plateau the source Taylor series gives, with
`s=sqrt(X^2+Z^2)`,

    u_0=(R Z,(R/sqrt(2))s,-R X),  T_0=2 pi/R.               (A3)

The Morse lemma makes every sufficiently small displayed pressure level a
regular closed meridional curve. Coarea gives
`T(p)=oint r ds/|grad p|` and
`dnu=T(p)dp d vartheta/(2 pi)`. These are meridional/axisymmetric return
statements; swirl need not close the full three-dimensional particle path.
Axisymmetry is exactly what makes the azimuthal advance irrelevant to the
return used here.

For the returning tangent covector, direct differentiation of the
incompressibility constraint gives the full-pressure amplitude system

    k_dot=-(grad u_0)^T k,
    A_dot=-(grad u_0)A
          +2 k[k dot (grad u_0)A]/|k|^2.                    (A4)

In the parallel meridional frame, `k=e_t` and `A(0)=e_n` yield

    M_0=[[1,0],[-sqrt(2) pi,1]],
    |M_0(1,0)|^2=1+2 pi^2.                                  (A5)

The sign, pressure factor two, determinant, Jordan rank, and gain agree with
the exact verifier. On a fixed small shell, scaled coefficient and period
continuity gives `M_sigma=M_0+O(sqrt(sigma))`. The exact action-angle
cotangent shear is retained: at finite `sigma` this is a map between nearby
polarization fibers, not a fixed Floquet matrix to be iterated.

There is one bounded normalization defect in the written WKB bridge. With a
fixed tube amplitude, the displayed generator

    a_N=N^-1 A exp(iN vartheta),
    q_N=J_*(0,a_N)=({xi,a_N},0)                              (A6)

has vorticity amplitude `O(1)` but, because Biot--Savart is order minus one,
velocity norm `O(N^-1)`, not the claimed unit leading velocity amplitude. An
absolute `O(N^-1)` Duhamel error is therefore not by itself small relative to
that unnormalized packet. The minimum repair is to set

    qhat_N=q_N/||v[q_N]||_2

(equivalently rescale (A6) by a factor comparable to `N`) and state the WKB
residual and Duhamel error as `O(N^-1)||vhat_N(0)||_2`. Scalar rescaling
preserves exact dynamic accessibility, axisymmetry, finite energy, the
translation quotient, and the energy-gain ratio. With this normalization the
standard fixed-time two-term geometric-optics estimate proves the exact
linear solution inequality claimed in `0019`.

This is a bounded proof wording/normalization repair, not contrary evidence
to the theorem. The actual `0019` result is only a linearized-Euler theorem.
It does not assert a finite-amplitude nonlinear perturbation estimate, so no
Kato theorem or nonlinear packet remainder is a missing dependency inside
Unit A. Constructing such a finite-amplitude transfer is a further frontier.

**Unit A verdict: established after the bounded normalization repair.** For
each sufficiently small fixed smooth plateau-shell Gavrilov member, there is
an exact finite-energy dynamically accessible axisymmetric solution of the
global linearized Euler equation, modulo axial translation, whose physical
kinetic energy gains the stated strict factor over one meridional circuit.
This proves finite-time transient amplification, not a growing eigenmode,
repeated-return bound, nonlinear instability, or all-time localization.

## Unit B: source equation, full contour domain, and signs

The García--Hassainia--Hmidi PDF and extracted text hashes are respectively
`7ba22ef57ba9453dd2631c33b1ad094a35ce2be171426cecce80fba8a5cdfac6`
and `f6461cfc68e2aa3b08cadc3409a36f4f0267c5a7c73fdb03e46f9636c121c871`.
Theorem 1.1 supplies the actual pair for fixed `0<a<b`, `kappa>0`, every
sufficiently small `epsilon`, and `lambda` in its Borel Cantor set. It gives
equal-strength `epsilon^-2` patches with sufficiently high Sobolev boundary
regularity and slow period `T`; the physical translating-frame period is
`mathcal T=T/abs(log epsilon)`.

With `J(a,b)=(-b,a)` and physical time, source equation (2.10) is exactly

    gamma_j,t dot J gamma_j,theta
      =partial_theta(Psi-c_epsilon rho)(gamma_j).             (B1)

For counterclockwise orientation `J gamma_theta` points inward, but both
sides of (B1) use that same convention. An outward normal displacement has
the positive domain derivative

    delta q_j=epsilon^-2 k_j delta_Gamma_j,
    delta Psi=C_G epsilon^-2 sum_l integral_Gamma_l G(x,y)k_l ds_y,
    C_G=1/(sqrt(2) pi).                                      (B2)

Direct differentiation of (B1) agrees with the Hamiltonian linearization
`eta_t={eta,h_*}+{q_*,delta Psi}`. Thus the translating sign, patch strength,
self/cross signs, and pressure density are consistent.

Normal graphs modulo tangential parametrization give the individual
fixed-area tangents `integral_Gamma_j k_j ds=0`. The radial-square source
variable satisfies

    k=epsilon^2 (Z dot n) delta f/w,                         (B3)

and strict star-shapedness makes this a Sobolev isomorphism. The exact area
formula shows that zero angular mean of `delta f` is precisely (B3)'s
fixed-area condition. Applying it independently retains both labels and all
first modes; the source's reversibility and half-period label relation select
the base orbit but are not imposed on IVP perturbations. The centroid
differential maps the two translation normal fields to the two coordinate
vectors, so the four-center plus two-shape split is a real local coordinate
split rather than a Fourier deletion.

## Unit B: IVP, monodromy, and reduction

The source's singular self-kernel calculation is local to one smooth contour.
Without using parity, it yields variable first-order transport plus an
order-zero Hilbert/logarithmic part and tame smoothing terms. Applying it to
each label preserves both self blocks. Positive inter-contour separation makes
each cross kernel smoothing. After the periodic trivialization the exact
linear generator therefore has the form

    A_*(t)k=B(t) partial_theta k+C(t) H k+R(t)k,              (B4)

on `X^s=H^s_0(T) direct-sum H^s_0(T)`, with real sufficiently regular
coefficients and `R(t)` bounded on `H^s`. No uniform-in-`epsilon` bound is
claimed. The source periodic Nash--Moser inverse is not used as an IVP
estimate.

The real transport identity, skew-adjoint Hilbert transform, the commutator
`[C,H]`, and bounded remainder give

    (1/2)d||k||_Hs^2/dt <= C_s(t)||k||_Hs^2.                 (B5)

For every index below the available base regularity, Friedrichs
regularization, compactness one derivative lower, uniqueness, and the reverse
equation construct

    U(t,r):X^s -> X^s,
    ||U(t,r)^(+/-1)|| <= exp(integral C_s),                  (B6)

with solutions in `C H^s intersect C^1 H^(s-1)`. Individual area conservation
makes `X^s` invariant. Periodicity of the labeled parametrization identifies
the endpoint tangent fiber with the initial one, so
`M=U(mathcal T,0)` is a bounded invertible endomorphism. A different periodic
tangential gauge conjugates this map and does not change the geometric result.

The normalized KKS form and `H_c=E-c_epsilon P` use the corrected sign
`i_X Omega=dH`. Multiplying `H`, `P`, `Omega`, and the action together by the
physical factor `2 pi rho_m` changes none of the equations, invariant kernels,
or symplectic identities. Hamiltonian linearization gives

    M^* Omega M=Omega,
    M v_phase=v_phase,  M v_z=v_z,
    dH_c M=dH_c,         dP M=dP.                            (B7)

The nonstationary radial exchange makes the phase and common axial-translation
vectors independent. Weak nondegeneracy then makes the two invariant
covectors independent, and smooth normalized companion directions exist.
Hence

    G=span(v_phase,v_z) subset
    K=ker dH_c intersect ker dP

is an invariant flag. The induced operator on `K/G` is bounded and invertible;
the displayed slice projection is only a representative of that quotient.
The companion block is retained, so neither Jordan shear nor physical
conserved-level directions are silently deleted. Symplecticity here is the
continuous weak KKS identity on the patch tangent; it is not a strong
symplectic spectral theorem.

**Unit B verdict: established as stated.** On every fixed sufficiently smooth
source pair, the complete independent-label, non-reversible, individual-area
Jacobi equation has a zero-loss bounded invertible propagator through one
physical period. Its monodromy is symplectic, fixes phase and common axial
translation, preserves energy and impulse covectors, and induces a bounded
invertible map on the stated reduced quotient. Constants may deteriorate as
`epsilon` tends to zero. Nothing locates the remaining spectrum or controls
powers of the monodromy.

## Oracle scope, findings, and separate dependencies

The primary PDFs and direct PDE/operator calculations are the strongest
oracles. The Unit A exact script correctly checks its finite-dimensional
deformation, pressure-amplitude, return, gain, and action-sign identities; it
does not check WKB normalization or the continuum residual. The Unit B exact
script correctly checks the quadratic-action sign, companion algebra, and
invariant-flag block; it does not construct the singular contour operator or
prove the Sobolev energy estimate. Both scripts and their captured zero exits
were inspected and not rerun because their predicates and code are unchanged.
No small numerical ratio occurs.

| Unit | Finding | Disposition |
| --- | --- | --- |
| A | The seed in (A6) has velocity norm `O(N^-1)`, while the prose and absolute error estimate treat it as unit size. | The bounded normalization/relative-error correction is accepted; theorem and route verdict are unchanged. |
| A | No finite-amplitude nonlinear theorem or public module/test is actually stated or inventoried. | Nonlinear transfer is frontier, and filename metadata supplies no evidence; neither is a support defect in the linear claim. |
| B | The source periodic inverse lives on a reversible delayed-label space. | Correctly excluded; direct full-label IVP energy construction supplies the claimed operator. |

Unit A's remaining in-scope dependency is the repeated moving-covector,
nonlocal Leray skew product over arbitrarily many circuits. A separate
finite-amplitude nonlinear transfer would require a normalized `H^s` family,
an amplitude-dependent common Euler existence interval, and a quadratic
remainder estimate, but it is outside the actual `0019` statement.

Unit B's remaining dependency is a reduced Krein/resolvent or modulated normal
form that controls powers of `M_red` while retaining first-mode/center and
companion blocks. Nonlinear orbital stability would additionally require a
nonlinear modulation bootstrap.

Neither verdict licenses all-time packet recurrence, spectral or power-bounded
GHH stability, nonlinear/orbital stability, unrestricted three-dimensional or
swirl perturbations, a mechanical force, particle identity, physical or
quantum spin, electron/neutrino completion, parent completion, or a global
no-go.

## Unit A bounded correction check

The author changed only the four registered `0019` correction artifacts. The
reviewed hashes are:

| Artifact | Before SHA-256 | After SHA-256 |
| --- | --- | --- |
| `0019/derivation.md` | `b7bc7926e75d1eaae1b681b509b7508063403757673b8d14a14c29914591bfba` | `1f56e289e85c0ffd463067d82da2ae05a0b5d9cfca275f0a1956529018f6e3ae` |
| `0019/result.yaml` | `e6fcc43e552d48d07376cfed47e4cb65bfaef2c39c7cd91bbcaed78d3450d757` | `b320c5216095f32f765e99fa1cde4bd61fd71be726b58fe7457926aa1d04acca` |
| `0019/validation.md` | `295237f75d31791fd1fdbabe85c14209729543e5cea2d8aed0538269bf91e79e` | `107423185d6750e7f85482ed7d3961d3f5a136757c335a0e13b928ea1aa433a6` |
| `0019/normalization-correction-receipt.md` | new | `45053dee0cbd61cb08419593900589e5c7d4725c7dd55b0b9861e4517d246d38` |

Equation (36) now defines
`qhat_N=q_N/||v[q_N]||_2`, explicitly notes that the raw velocity norm is
comparable to `N^-1`, and explains that scalar normalization preserves the
coadjoint tangent, symmetry, finite energy, and gain ratio. Equation (37) gives
both the WKB residual and exact-semigroup Duhamel error as
`O(N^-1)||vhat_N(0)||_2`. The result and validation artifacts state the same
normalization and relative-error boundary. This fully closes the scientific
omission found by Unit A.

The receipt's three before/after hashes and the unchanged hashes of
`verify_return.py`, its command, stdout, and zero exit agree with the frozen
0026 boundary. The oracle was correctly not rerun because it does not test the
continuum packet normalization and none of its predicates changed. Active
attempt `0025` was neither inspected nor adjudicated.

One minimum prose cleanup remains in `derivation.md`: “squared-norm relative
squared-norm error” should read “relative squared-norm error.” The duplicated
words do not change the displayed relative estimates, the argument, or the
accepted verdict, and no further scientific correction check is needed.
