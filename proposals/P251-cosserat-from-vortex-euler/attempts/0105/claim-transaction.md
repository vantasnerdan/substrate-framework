# P251 parent claim transaction — proposed, not yet accepted

This binds the previously unused identifiers C-CST-008,009,010 to the
statements below. Original C-CST-001..007 remain bound to their historical
statements; the repaired seven-node scientific objective is not replaced
by reviving the refuted straight-tube alpha formula. Review status is
unaudited. The implementation boundary will be frozen by content hash
after0107's canonical coefficient API is complete.

## C-CST-008 — smooth stationary compact Euler angular cells

Fix rho>0 and a nonzero curl eigenvalue lambda. Within the declared
stationary finite-variance isotropic Gaussian Beltrami law, there is a
positive-density stationary marked selection of actual invariant knotted
vortex-tube neighborhoods with finite volume and a continuous ambient
complement. On a uniformly bounded good-patch event one can construct
smooth compact isovortical core/reaction fields Q,S wholly inside each
tube, whose full-space induced velocities are also compact, with:

    core_angle(Q)=1, core_angle(S)=0,
    Omega(Q,S)=B!=0,
    H=[[Hq,N],[N,P]] positive definite,
    G(Q)=(B^2/P)n, G(S)=0, L(Q)=0, L(S)=B n.

Here G=rho integral r cross xi and L=rho integral r cross v_xi are
different measured moments, Omega is the Euler KKS form, and H is the
complete isovortical second variation of Euler kinetic energy. Every
quantity is constructed from actual field jets, compact profiles and
Euler integrals. No core mass, locking constant or desired frequency is
supplied. The core angle is a normalized local vorticity-direction
observation. It is not asserted to be a global compact symmetry action.

The same event supports additional zero-G/zero-L compact pairs with
positive reduced stiffness/inertia ratio above any prescribed finite
bound. Their supports can be made disjoint, so the full H/KKS locality
used for their reaction blocks is exact. Finite radii, profile norms,
positive margins, laws and number densities are declared geometric
inputs; no universal five-scalar filament formula is claimed.

Exact proof attachments:0085 compact-spin and positive-pair proofs;
0103 six-moment theorem;0102 section1 triangular completion;0098 bounded
assembly and measurable invariant-domain selection;0071 Gaussian-law
construction and its previously audited primary support inputs. The
exact finite-field minor has its explicit characteristic-zero rank bridge,
and its canonical coefficients are independently tested by0106.
Public implementation: euler_compact.py and its tests; euler_fourier.py
provides distinct signed Euler forms and independent normalization tests.

## C-CST-009 — conditional one-action Euler-to-Cosserat continuum

Use C-CST-008's cells in the explicit phase Cauchy--Born embedding of the
quadratic Euler material action,

    eta=U+Ez,
    pi=rho[(u0.grad)eta+V+A Ez],
    Axi=P(xi cross omega0)-curl(P(xi cross omega0))/lambda.

All ambient fluid is retained. U is the declared coherent tube-centroid
plus continuous-ambient coordinate; Phi is the registered local core-angle
field, q=Phi-curl U/2. The ensemble pairs time-reversed and reflected full
states, with one COMMON macro conjugate variable V and independently varied
microscopic reaction momenta. Unresolved nonaffine variation/relaxation is
excluded by the declared Cauchy--Born closure, not silently made into gauge.
The law and marks are fixed before deformation. The claim concerns the
linear continuum and its full second spatial jet, not all wave numbers.

Pulling back that ONE action, taking its stated ensemble average and
eliminating the complete retained momenta yields positive computed

    j=nu E_Palm[B^2/P]/3,
    kappa=nu E_Palm[Hq-N^2/P]/3, alpha=kappa/4.

The six-moment match gives b=0 in the leading physical mixed mass and
ell=g-kappa b/j=-kappa/2!=0. Thus the optical branch has nonzero actual
coarse translation/core-angle transfer
U/Phi=-j sigma |k|/(2rho)+O(|k|^3) in curl helicity sigma. The leading
mass density is the total rho, derived also by the exact common Galilean
translation of the same stationary fluid, not a tube filling fraction.

Zero-moment compact STF and neighbor-angle attachments, with finite
structurally selected amplitudes, give positive shear mu and positive
transverse/longitudinal spin curvatures C_T,C_L after retaining all added
gradient inertia. Their coefficients are full canonical action/Schur
integrals, including the negative bare material covariance stiffness,
potential reaction squares and the surviving STF rate-source norm.
The SAME derivative field normalization transforms both kinetic and
potential forms. In its stated convention the equations are

    rho U_N,tt=(mu+alpha)Delta U_N+2alpha curl Phi_N-grad p,
    div U_N=0,
    j Phi_N,tt=C_T Delta Phi_N+(C_L-C_T)grad div Phi_N
                    +2alpha curl U_N-4alpha Phi_N.

These equations and their angular-momentum/couple-stress balance are exact
coefficient identities for the declared second-gradient conditional action.
The optical gap is4alpha/j; the acoustic transverse speed squared ismu/rho.
There is no incompressible longitudinal displacement wave. The physical
gradient masses, finite-radius Fourier filters and macro affine tube-spin
row accompany the field map. In particular full tube spin is not renamed
j Phi_dot: its additional affine/shape current is explicitly retained.

Removing the core-angle and angle-gradient populations deletes Phi before
any zero-inertia division and leaves a positive incompressible Navier--Cauchy
sector when the independent STF population remains. Removing every attached
population and the vortical covariance instead leaves linearized Euler.
The explicit limiting families and unchanged-density bookkeeping are in
0105/population-limits.md.

This is an exact conditional variational continuum construction from the
Euler material action. It does NOT claim that unrestricted Euler trajectories
remain in the chosen finite phase family. Their reconstruction residual is
R=Udot-V+Ezdot-AEz;0095 supplies its exact complement/memory and the matching
physical observation correction. Its vanishing is not an unstated premise
or a conclusion of the finite KKS calculation. The stronger unrestricted
realization question remains distinct, with0101 as finite-time progress.

Exact proof attachments:0097 canonical phase pullback and signed symplectic
oracle;0102 complete physical coefficient join and its26-check receipt;
0098 observable map;0105 limits; existing micropolar.py balance/field-map
tests and repaired CST004/005 exact operator/dispersion evidence. The executed
DOP853 refinement is a reduced-equation regression, not an Euler-existence
oracle. Public implementation includes the canonical Euler and micropolar
modules and0107's exact paired reaction reduction.

## C-CST-010 — explicit orientation-ergodic conservative contrast

For an integrable finite-range orientation energy of complete unresolved
states with an independent local product-Haar law, independent coherent
frame shifts leave its averaged energy invariant. Consequently its
coherent conservative torque and couple stiffness vanish, while quadratic
fluctuation energy need not vanish. Uniform one-point marginals alone do
not imply the result: correlated phases supply an explicit counterexample.
A no-retained-coherent-angular-current closure is separately stated for
the dynamical no-spin limit; static isotropy alone does not remove memory
or an angular current.

Exact proof:0058's Haar change of variables and canonical phase integration.
Corroboration: the already executed200000-sample signed-response Monte Carlo
with its declared seed, error model and bias mutation, at its narrower
first-moment scope. This does not assert that the coherent law of C-CST-009
is product Haar or that stationary Gaussian Euler isotropy by itself
erases its positive locking. It states the separately declared contrast
ensemble premise and its consequence, without deleting fluctuations.

## Frozen review and impact policy

Review these three proposed claims individually as one transaction. Their
unchanged historical/source receipts are inputs at their recorded scopes,
not invitations to re-adjudicate every prior attempt. Main's accepted
release and all pre-existing accepted claims remain unchanged. Direct
consumers are the Euler/compact/micropolar tests, named construction
verifiers and the proposed continuum claims; no accepted sector imports
these unpromoted APIs as authority. One substantive independent review
and one directly affected correction check are planned after the full
content boundary is frozen. No full-goal or PR-ready verdict is asserted
by this draft.
