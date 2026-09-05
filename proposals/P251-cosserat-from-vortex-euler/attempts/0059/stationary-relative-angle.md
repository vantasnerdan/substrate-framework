# Stationary Euler action density and the affine-relative-angle replacement

## 1. What the source actually supplies

The archived source is 2006.15033, *Beltrami fields exhibit knots and chaos
almost surely*, by Enciso, Peralta-Salas and Romaniega. Proposition 3.2 gives
smooth realizations; Proposition 3.4 gives the translation-invariant spectral
covariance; Corollary 3.5 gives `E[u(x) tensor u(x)]=I`; Proposition 3.7 gives
translation ergodicity for integrable observables; Proposition 3.8 gives
positive probability for every finite C^m neighborhood of every Beltrami
field. Theorem 1.2 gives positive volumetric lower bounds for knotted tori.

The particular Gaussian law in the paper is NOT asserted to be rotationally
invariant. Its spectral density contains the polynomial p of Eq. (2.4).
Equal point covariance I alone does not prove isotropy of the full law.
One may average the law and all marks over rotations; this is a declared
ensemble average, not a claim that the resulting mixture remains translation
ergodic. Likewise reflection partners have opposite lambda and must be
reflected as full fields. Below stationary identities hold for each component;
isotropic/parity-even coefficients refer to the explicitly averaged ensemble.

Select good patches by a translation-covariant hard-core rule with bounded
coherence radius. A concrete positive-intensity rule is an independent
Poisson candidate process: retain an isolated candidate only if its finite
neighborhood is C^m-close to the prescribed good prototype. Isolation has
positive probability and is independent of the field event, which has
positive probability by Proposition 3.8. Hard-core separation makes supports
of different patch generators disjoint. Restrict the event to finite C^m
bounds and a uniform inverse-Gram bound; these are open finite-patch bounds,
not a selected value of a desired modulus. Conditional componentwise
ergodicity or the ensemble expectation gives the corresponding densities.

Marks and patch identities are frozen in the reference configuration and
transported under the declared affine deformation. They are not reselected
after a q perturbation. Reselecting them would change the physical ensemble
and could erase the very relative-angle observable being measured.

## 2. The global common-rotation density: a useful identity, not a rotor proof

For patch a centered at X_a, let eta_a be compact and suppose its three
Euclidean translation moments vanish:

    T_a=rho integral eta_a cross omega=0.

For K_e=e cross x, the exact decomposition is

    Omega(K_e,eta_a)
      =rho integral omega.[(e cross (x-X_a)) cross eta_a]
         +(e cross X_a).T_a.

Thus the pairing is independent of the origin and has a finite density
`nu E_Palm[b_a]` when the local pairing is integrable. This is stronger than
hoping the origin-dependent part cancels only in expectation. It follows
directly from the compact moment constraint.

However, this identity alone does not establish an independent material B.
The actual global K also moves the material centers by e cross X_a. A
time-dependent rigid material rotation has orbital kinetic energy density
proportional to the mean squared distance from its rotation axis, which
grows like the observation radius squared. A finite spin KKS pairing does
not remove that orbital energy. Alternatively, if one quotients by rotations
of the entire realization and all marks, B can become merely a change of
realization label. Invariance of the co-rotated probability law cannot decide
between these interpretations. A physical global rotation is not a gauge
when its spatial/material registration is retained, but neither is it an
independent local microrotation with a finite supplied common-body inertia.

**Global-B density route verdict: blocked at the stated material-identification
step, not a no-go for the parent object.** The exact finite pairing above is
retained. The replacement below does not require the extra independent B.

## 3. Eleven exact affine moments, with a constructive right inverse

Use the original target's physical fields only: U and Phi. In a local affine
material frame put beta=curl U/2, and let q be the physical core angle
RELATIVE to that frame. The needed microscopic action is the already
constructed positive compact angle/cage sector, including its pairing with
the same Euler affine displacement. Those pairings are computed here rather
than declared absent.

For an incompressible affine generator

    X_{a,H}(x)=a+H(x-X_patch),  tr H=0,

there are eleven independent parameters: three translations and eight
trace-free linear fields. Define their local moment functionals

    F_{a,H}(xi)=Omega(X_{a,H},xi)
               =rho integral xi.[omega cross X_{a,H}].

They are compact integrals even though X is affine. Their curl kernels are
minus the affine vorticity variation,
`curl(omega cross X)=(X.grad)omega-H omega`.

Here is an explicit full-rank prototype, extending the six-Euclidean-rank
construction of 0057. Start from a decaying finite-Herglotz EPS field with
the desired robust tube, and add an arbitrarily small same-lambda ABC field

    u_ABC=(sin(lambda z)+cos(lambda y),
           sin(lambda x)+cos(lambda z),
           sin(lambda y)+cos(lambda x)).

It is exactly Beltrami and has nonzero Fourier atoms at all three pairs
±lambda e_i. The perturbation can be small enough to retain the robust
selected EPS structure. The original EPS spectrum is a smooth density on
the sphere, not a point atom.

If an affine X generated a zero vorticity variation, its Fourier action on
the atom at k_i would have a derivative-of-delta term with coefficient
H^T k_i. Smooth sphere density, even after affine differentiation, cannot
cancel a derivative point atom. This can be checked with shrinking test
functions: the point derivative scales as epsilon^-1 while the differentiated
smooth surface term scales at most as epsilon. Thus H^T k_i=0 for all three
independent k_i, giving H=0. The remaining atom coefficients are proportional
to a.k_i, so a=0. All eleven curl kernels are therefore independent. They
are analytic, so they are independent on every open response ball.

Choose eleven pairwise disjoint off-core response balls, separated from the
raw physical core and cage supports. The exact Gram construction in each
ball is

    f_A=rho omega cross X_A,
    G^(j)_{AB}=integral chi_j curl f_A.curl f_B,
    eta^(j)=sum_B curl(chi_j curl f_B) [(G^(j))^-1]_{Bj}.

Each G^(j)>0 and `F_A(eta^(j))=delta_Aj`. Distinct responses have disjoint
supports, so their span is KKS-isotropic. These finite Gram inequalities
persist in an open C^m neighborhood of the prototype; the source full-support
theorem makes that good-patch event have positive probability. All coefficient
and derivative bounds can be chosen uniformly on a smaller such event.

Let R F=sum eta^(j) F_j and Pi_F=I-R F. For the raw physical angle jet plus
negative-helicity cage Q0=Q_R+C1, and its cage partner S0=C2, define

    Q=Pi_F Q0,  S=Pi_F S0.

Then EVERY affine moment vanishes exactly:

    Omega(X_{a,H},Q)=Omega(X_{a,H},S)=0.

The response supports do not touch the physical jet. Moreover the response
span is isotropic and disjoint from Q0,S0, so
`Omega(Q,S)=Omega(Q0,S0)` EXACTLY. Thus the full selected Euler KKS action
has no missing affine/Q or affine/S cross term. This is not a generic
projection identity; the eleven separate response supports are essential.
The fixed corrections have bounded norms, so the finite negative-helicity
construction still gives a positive full H on Q,S at a sufficiently large
finite carrier, uniformly on the selected good event.

There is also an exact force-moment consequence, useful for the full slow
spatial return construction in 0057. For F=xi cross omega, the eleven
constraints imply

    integral F=0,  M_ij=integral y_j F_i=c delta_ij.

The first equation is the three translation moments. Orthogonality of M to
every trace-free H gives the second. Choose a smooth compact f with
`integral f=-c`; then `integral y_j partial_i f=c delta_ij`. Consequently
G=F-grad f has zero zeroth AND first moments, while P0G=P0F and the integrated
helicity is unchanged. This is an explicit moment repair, not an infrared
cutoff. The averaged-center Taylor construction in 0057 represents this G
as a smooth compact double divergence; its full spectral energy symbol then
has the regularity needed for the complete second-gradient action jets.

## 4. Full stationary Leray action density, including all tails

Let Xi be the stationary sum of the selected compact patch generators and
F=Xi cross omega. Uniform local bounds and hard-core supports give finite
second moments for F and curl F. The three translation constraints imply
the mean source E F=0 by the patch integral identity. Use the stochastic
translation spectral representation to define P0: the usual Leray symbol
at nonzero frequency and zero on the constant atom. No mean Galilean velocity
is being silently inserted into the spin sector.

The exact quadratic Euler-orbit energy density is

    H_density=rho E[|P0 F|²-F.curl F/lambda].

This keeps the FULL projection, including the velocity induced outside every
patch and the cross energies of different patches. It is not replaced by
the sum of isolated patch kinetic energies. Stationary integration by parts
gives `E[(P0F).curl(P0F)]=E[F.curl F]`, so this is the same Hessian form as
the finite compact construction, now with the stationary trace E.

For completeness, the high-carrier error estimate can be made before this
global projection, avoiding a locality assumption. In one carrier frame
write a complex local force `(F0+F1/k) exp(ikz)`, Pi=I-ez ez^T, and set

    A_app=(i/k) ez cross (Pi F0) exp(ikz),
    pi_app=F0_z exp(ikz)/(i k).

Then EXACTLY

    F=curl A_app+grad pi_app+R/k,
    R=[F1-i curl_slow(ez cross Pi F0)+i grad F0_z] exp(ikz).

Every potential and remainder is compact in its patch. Summing them over
the disjoint stationary supports preserves this identity. P0 is an L²
contraction, annihilates the gradient, and fixes the curl, hence the global
projection error relative to the principal transverse carrier is O(1/|k|)
in mean-square norm. Its curl error is bounded independently of k, because
curl P0=curl and the fast derivative of F0_z ez cancels. These are the same
D/|k| and E bounds used in 0045, now obtained for the full coherent stationary
sum. All inter-patch kinetic interactions remain in that bounded remainder.

The local negative-helicity principal quadratic form is positive and its
support densities have a positive lower expectation on the good event.
Fixed moment-response corrections are bounded, with their mixed helicity
terms bounded by stationary integration by parts. Therefore a finite uniform
carrier gives a strictly positive full internal H density. Its KKS density
is the corresponding positive patch pairing density. This is a genuine
positive microscopic relative-angle action density, not a positive local
mode followed by an uncontrolled coherent kinetic summation.

## 5. Lie–Poisson structure and the same mean-mass normalization

For stationary smooth fields with the stated integrable products, the
translation derivatives D_i commute and `E[D_i f]=0`. The divergence-free
vector-field commutator is a Lie bracket. Its usual Euler Lie–Poisson bracket
on differentiable density functionals is

    {A,B}=rho E[omega.(v_A cross v_B)]

Here the normalization is `delta A=rho E[v_A.delta u]`. This is
the negative commutator pairing obtained from
`curl(v_A cross v_B)=-[v_A,v_B]`; Jacobi follows from the vector-field Jacobi
identity and the stationary integration-by-parts trace, exactly as for the
ordinary Euler bracket. For E=rho E|u|²/2 it gives the Euler generator
`u_t=P0(u cross omega)` in the zero-mean sector. Equivalently,
`Omega(X,eta)=rho E[omega.(X cross eta)]` has `i_X Omega=dE` for X=u.

The selected finite-dimensional KKS form is the pullback of this closed form
along the actual composed volume-preserving generator flows. The moment
projection describes its tangent at the stationary background; it does not
claim that an arbitrary linear kernel of moment functionals is itself a Lie
subalgebra. No Jacobi identity is inferred just from a constant test matrix.

The mean velocity is a different, explicitly retained component. Writing
`u=V+w`, `E w=0`, the SAME Euler kinetic density decomposes exactly as

    rho E|u|²/2=rho |V|²/2+rho E|w|²/2.

The source normalization gives E u0=0 and E|u0|²=3 in its units. Thus a
uniform Galilean velocity changes the energy DENSITY by rho |V|²/2, not by
an assigned finite whole-space energy. P0 keeps the internal perturbations
mean-free, so the mean/internal kinetic cross is exactly zero. The mass
density is rho for ALL the fluid, not nu times the mass of selected good
patches. The latter intensity weights internal action coefficients only.

This uniform mean is retained in a co-moving material frame:
`u(x,t)=V(t)+w(x-X(t),t)`, `Xdot=V`. In the Euler material derivative the
`-V.grad w` coordinate derivative cancels the `V.grad w` convection term.
The uniform mean acceleration and pressure reaction remain in the mean
equation. Thus the mean/internal split is compatible with the uniform
Galilean/material action, not just an isolated energy identity. Extending
it to the declared slowly varying affine field retains the pressure/current
terms in the parent assembly; this paragraph does not discard those terms.

A material partition has the same total mass/momentum density when every
parcel and exterior return is counted. Its mass-centroid and physical spin
currents are related to the canonical density fields by the explicit
boundary/ambient improvements of 0052/0055. Those observable filters and
pressure fluxes are retained by 0057/0053; this mean normalization does not
replace a parcel centroid by a vortex-center coordinate without explanation.

## 6. Full reaction-space reduction and the physical relative angle

Retain independent fluid reaction amplitudes, including their coherent
global Leray interaction. Let the reaction Hilbert space H have the stationary
energy-density norm (for patch amplitudes, intensity times the Palm L² norm).
The physical relative-angle vector is q in R³. The full first-order action is

    L_int=<s,D qdot>-<s,P s>/2-<s,N q>-q.H_QQ q/2.

D maps q into the actual KKS moment couplings on all patches; P is the FULL
reaction Hessian operator, and N is its mixed angle/reaction block. Every
affine symplectic cross is zero by section 3. The uniform high-carrier bounds
of section 4 apply to arbitrary square-integrable reaction amplitudes, not
only one common scalar: disjoint supports control their squared norm and the
global P0 contraction controls all induced velocity cross terms. They make
the complete (q,s) form positive and P bounded and coercive. D has rank three
after the nondegenerate oriented patch ensemble is included. Its coefficients
come from the local KKS integrals, not a desired mass.

For the two time-reversed realizations, vary the reactions independently:

    L_pair=[L_plus(q,s_plus)+L_minus(q,s_minus)]/2,
    s_plus=P^-1(D qdot-N q),
    s_minus=P^-1(-D qdot-N q).

Elimination gives the exact full-operator result

    L_pair=qdot.J qdot/2-q.K q/2,
    J=D* P^-1 D>0,
    K=H_QQ-N* P^-1 N>0.

Both operators are computed with the complete stationary fluid interaction.
The positive signs follow respectively from coercivity/rank and the positive
full Schur complement. The odd gyroscopic term cancels only after these
independent variations. A single scalar formula B_density²/Hss_density is
valid only for an explicitly constrained single coherent reaction coordinate;
it is not substituted for this inverse operator or for independent patch
reactions. In particular an average of isolated-cell inverses is not assumed.

Simultaneous orientation averaging of the full action and all independent
reaction fields gives `J=j I`, `K=4 alpha I` with j,alpha positive. These
are the actual full density coefficients. The familiar `j=nu I_cell/3`,
`alpha=nu K_cell/12` is recovered only when independent-cell factorization
has separately been established; it is not needed here.

The marked physical core angle is `Phi=beta+q`, `beta=curl U/2`. It is
observable relative to the transported affine material frame. Unlike a
rotation of the entire probability law, changing q changes that relative
core jet and has the strictly positive restoring coefficient just derived.
Under a common rigid frame rotation, Phi and beta change together and q
does not. No independent cyclic common-body coordinate is needed.

With the full density coefficients just derived, the physical kinetic and
locking terms are

    T=rho |Udot|²/2+j |Phi_dot-curl Udot/2|²/2,
    W_lock=alpha |curl U-2 Phi|²/2.

For a transverse curl helicity h, the kinetic matrix is exactly

    M=[[rho+j k²/4, -j h k/2],[-j h k/2,j]],
    det M=rho j>0.

This is the gradient kinetic cross already included in the 0053 normal form;
it is not discarded. In that normal form b=-j/2, m_U=j/4, m_Phi=0, giving
`C_eff=C-alpha j/rho` for the transverse spin-gradient coefficient. Any
additional finite Euler affine/gradient contributions remain in the full
matrix before the same normalization. The existing positive gradient-cage
construction must dominate the resulting finite correction, not a copied
unmixed optical slope.

Integration-by-parts and field-map boundary currents matter for a globally
rigid affine motion: a local bulk normal form is not a license to discard
its boundary angular momentum. The physical-to-canonical maps of 0052/0055
remain part of the reported observable. This is precisely the original
slow-affine/second-gradient scope, not an all-wavelength microscopic claim.

## Route result

The unnecessary independent global-B construction is not used to close the
physical angle obligation. The replacement establishes a positive compact
relative-angle action with an exact eleven-affine-moment KKS projection,
full stationary Leray energy density, and the derived relative-rate kinetic
coupling to the physical affine frame. It supplies the microscopic/action
inputs for the parent material-density and spatial-gradient assembly; it
does not claim those remaining assembly computations have already vanished.
