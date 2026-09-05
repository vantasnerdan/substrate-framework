# One-action compact Euler coefficient join

## 1. The two moments and their constructive completion

Use0097's canonical Jacobi action and0103's six independent compact moment
responses INSIDE the actual invariant EPS tube. For a compact generator xi
whose induced velocity v=xi cross omega is also compact and solenoidal,
define distinct vector rows

    G(xi)=rho integral r cross xi,
    L(xi)=rho integral r cross v.

Both xi and v have zero integral. Since omega=lambda u and v=xi cross
omega WITHOUT a pressure tail,

    integral xi cross u=0,
    G(B0 xi)=0, G(A xi)=L(xi),
    B0=u.grad, A xi=v-curl(v)/lambda.            (1)

The second equality follows by integrating r cross(u.grad xi); the third
uses integral r cross curl(v)=2 integral v=0. These identities are signed
and linear, not consequences of the norm identity H-K=rho||A xi||².

0103 supplies six compact response fields eta_Gi,eta_Li with their G/L
rows the six unit vectors. Choose their six support balls mutually
disjoint, and away from a physical core-angle Q0 and the positive compact
WKB pair Qc,Sc. Every response-response KKS entry vanishes exactly: distinct
supports are disjoint and each self pairing is zero. All raw-response
KKS entries vanish for the same support reason.

Let Bc=Omega(Qc,Sc)!=0 and choose a fixed cage amplitude a>0. FIRST set

    B=a Bc,
    S=Sc+eta_L(B n-L(Sc))-eta_G G(Sc),
    P=H(S,S), j_cell=B²/P.

Then set Qraw=Q0+a Qc and

    Q=Qraw-eta_L L(Qraw)+eta_G(j_cell n-G(Qraw)).            (2)

The notation eta_G c means sum_i c_i eta_Gi, and similarly for L.
Equations(2) give EXACTLY

    G(Q)=j_cell n, G(S)=0,
    L(Q)=0, L(S)=B n, Omega(Q,S)=B.

There is no implicit fitted inertia or circular root: P is computed BEFORE
the Q correction and depends only on S. The Q correction does not change
B or P. The moment matching is therefore an explicit triangular finite
construction, stronger than an existence claim for an unspecified root.

For positivity choose the finite a first so a²h0 dominates the finite
core/response energy at target G=0. Then use0085's analytic finite-carrier
bounds. G(Sc),L(Sc),G(Qc),L(Qc) decay by repeated compact integration;
B=O(a/k), P=h0+O(a²/k²)>0 and j_cell=O(a²/k²). The additional Q moment
correction in(2) is consequently O(a²/k²); its H cost and cross with the
fixed core response tend to zero. The FULL2x2 H remains positive, with

    H=[[Hq,N],[N,P]], kappa_cell=Hq-N²/P>0.                 (3)

No N entry is dropped. All choices are finite once the exact positive
symbol and response-norm margins are fixed. This is the same hierarchy
as0085 with the extra six-moment rows now computed, not a material-Hessian
substitution. Core observation derivatives are unaffected by the off-core
responses. Time reversal uses the same Q,S, because j_cell is even, G is
even and L,B both change sign.

## 2. Exact canonical macro momentum reduction

Write z for the compact coordinates, z=(q,s) in one cell, and let E be
their material-generator columns. The SAME phase embedding is

    eta=U+Ez,
    pi=rho[B0 eta+V+A Ez].

With a coherent macro space G, V belongs to G and P_G is its actual
kinetic Gram projection. Eliminating V from0097's canonical one-form and
Hamiltonian gives, without a second added action,

    L=rho/2 ||P_G(eta_dot-AEz)||²
      +rho <(B0 eta+AEz),eta_dot>
      -K(eta)/2-rho||AEz||²/2.                (4)

The identity H(Ez)=K(Ez)+rho||AEz||² supplies the positive coadjoint
internal Hamiltonian within this ONE canonical material action. It does
not assert K(Ez)=H(Ez). The bare macro K remains the actual negative
covariance stiffness derived in0097, prior to the shear repair below.

For an affine macro test, compactness gives

    rho<U,Ez>=beta.G(Ez),
    rho<U,AEz>=beta.L(Ez), beta=curl U/2.                  (5)

All STF and constant parts vanish in these kinetic pairings. At long wave
number the corresponding coherent profiles are

    P_G Ez = curl(G_density(Ez))/(2rho)+O(k²),
    P_G AEz = curl(L_density(Ez))/(2rho)+O(k²).

The minus cross from the square in(4) cancels the positive
rho<AEz,Udot> in its one-form. Thus the old L-only inference of an
absolute-angle symplectic term is not the complete calculation. The
remaining leading mixed kinetic term is

    beta_dot.G_density(E zdot).                           (6)

For GQ=C, GS=D, this is beta_dot.(C qdot+D sdot). If D is not zero,
integrating by parts produces a beta_ddot forcing on s; dropping it
would change the temporal polynomial pencil. Construction(2) makes D=0
EXACTLY and sets C equal to the inertia computed from the same reaction.

## 3. Shared mean and independent time-reversal reactions

The macro mean and its conjugate V are COMMON across the equally weighted
u and -u realizations. Their microscopic reactions s_+,s_- are independent.
Average the complete canonical action before eliminating that common V;
this is not a separate momentum reset in each realization. Put

    p=(s_+-s_-)/2, r=(s_++s_-)/2.

Then the internal first-order/potential terms are

    B p qdot-[P(p²+r²)+2N r q+Hq q²]/2.

N couples the even reaction r, not the physical momentum p. Eliminating
r gives kappa_cell in(3), while preserving the Bp kinetic reaction.
This is why averaging the KKS coefficient B to zero BEFORE varying the
two reactions would give a false zero inertia.

For a transverse curl helicity sigma, put t=sigma |k|/2. With density
normalized coefficients C,j,B,P, the complete leading mean contribution
and momentum sector have the particularly exposing exact finite-block form

    L= rho Udot²/2+C t Udot qdot+B p qdot
       +t²(C qdot-Bp)²/(2rho)-P p²/2-kappa q²/2.           (7)

This display retains the ENTIRE mean square, not just its O(k) cross.
For general C, p elimination yields

    j_eff = t² C²/rho
       +B²(1-t² C/rho)²/(P-t² B²/rho),
    j_eff=j+t²(C-j)²/rho+O(k⁴), j=B²/P.                  (8)

When C=j, the stronger exact identities are

    p=B qdot/P, j_eff=j,
    C qdot-Bp=0.                                         (9)

Thus the moment completion cancels the mean reconstruction moment itself,
not just its contribution to a norm. The even reaction r=-Nq/P retains
the actual static response; its circulation-sign-odd spin cancels in the
paired physical mean. Eliminating separate realization-wise macro V's
would give a different N-dependent gradient correction. That is not the
shared coherent mean ensemble defined here.

No isolated-cell inverse is used in passing to a population. Let P be the
ACTUAL positive reaction operator (block diagonal only where compact
supports prove it), D the momentum-to-angle map, and

    J=D P^-1 D*, C=J.

The correction square is rho||R(J qdot-Dp)||²/2, where R is the coherent
curl/mean operator. The complete stationarity equation is solved exactly
by p=P^-1 D* qdot, and its residual square vanishes. Uniqueness follows
on the declared slow-wave range from positivity of
P-rho D*R*R D. A sufficient isotropic bound is j|k|²/(4rho)<1.
This finite long-wave neighborhood follows from bounded good-patch
coefficients; it is not an arbitrary finite-k Euler closure requirement.
Local matching GQ=j_cell makes C=J after Palm/isotropic averaging, without
replacing E[B²/P] by a ratio of averaged coefficients.
Equations(7)--(9) display the angle-only source. The additional STF rate
source retained below changes individual cell momenta; its mixed axial
projection vanishes only after the stated isotropic average.

## 4. Physical-angle jet and the nonzero coupling invariant

In the affine physical frame q=Phi-beta. After(9), the kinetic action is

    T=rho|Udot|²/2+j|qdot|²/2+j beta_dot.qdot
     =rho|Udot|²/2+j|Phi_dot|²/2-j|beta_dot|²/2.           (10)

Thus the leading physical mixed mass b is ZERO; it was not zero before
the moment/mean calculation. The remaining U-gradient mass -j k²/4 is
retained, together with any computed macro and strain-cage gradient mass.
It cannot change the positive zeroth rho block or the acoustic k²
coefficient, and is included in the derivative field normalization.

Common rigid rotation has zero material cross K(K_e,E) by direct
integration using stationary Euler; the pressure and transport terms
cancel. The bare macro covariance action is rewritten as symmetric strain
using its explicitly retained boundary null Lagrangian, not by asserting
its unbounded rigid-rotation density vanishes pointwise. The physical
objective locking term is kappa|Phi-beta|²/2. Consequently the transverse
leading matrices, before the remaining finite gradient corrections, are

    M=[[rho-j k²/4,0],[0,j]],
    K=[[(mu+kappa/4)k²,-kappa sigma k/2],
       [-kappa sigma k/2,kappa+C_T k²]].

The physical translation/core-angle coupling invariant is therefore

    ell=g-kappa b/j=-kappa/2 !=0.                         (11)

In particular U/Phi=-j sigma k/(2rho)+O(k³) on the optical branch in
these specified coarse coordinates. For unmatched C the same calculation
gives b=(C-j)/2 and ell=-kappa C/(2j); C=0 reproduces the old separable
zero-transfer case. These are calculated field-response statements, not
a relabeling that manufactures a nonzero transfer from a zero one.

## 5. Full positive shear and spin gradients on this SAME action

Use0103 to project BOTH G and L to zero on0096's auxiliary compact WKB
cage pairs. The extra generator-moment projections, like the velocity-
moment projections, decay by compact oscillatory integration, preserve
the exact B and leave the finite positive H margin intact. Thus their
full reduced ratio still obeys

    K_g/J_g=det(H_g)/B_g² >= c k_carrier².

Their zeroth and first generator AND velocity moments vanish. This is
the moment condition needed for the canonical mean square(4), stronger
than checking only the velocity row. A gradient attachment then has
coherent E and AE mean only at O(k³), so it changes no leading b or(11).

All actual canonical cross blocks are retained. In particular macro
forcing of a new reaction contributes the full potential square
(N a+C_aff E)^T P^-1(N a+C_aff E). Its fixed and linear-in-amplitude
pieces are computed before choosing any shear amplitude. With the new
cages support-disjoint from the physical angle profiles, the leading
STF self coefficient is the same positive FULL Schur stiffness K_g,
not the bare material K or Hq alone:

    mu(t_s)=mu_fixed+b_s t_s
                 +nu_s E[K_g] t_s²/10.

Here mu_fixed includes the bare material covariance and all fixed
reaction contributions. Its corresponding gradient mass is
nu_s E[J_g] t_s² |dev Edot|²/10 and remains in M_UU.

There is ALSO a distinct canonical STF RATE source; G/L normalization
does not make it zero. Integrating the two B0 mixed one-form terms once
in time gives, for an actual symmetric tracefree affine macro field,

    C_rate(xi):Edot = -2rho integral xi.(Edot u0).          (12)

This row is odd under u0 -> -u0. Thus its Qq contribution cancels but
its S contribution couples to the odd reaction p and SURVIVES the paired
action. It is separate from the even POTENTIAL source C_aff E. Before
eliminating p the total kinetic source is

    p^T [D* qdot+C_rate Edot]

and a strain-attached coordinate adds its actual B a_dot to this source.
The complete reduced kinetic contribution is therefore

    [D* qdot+C_rate Edot+B a_dot]^T P^-1
      [D* qdot+C_rate Edot+B a_dot]/2.                     (13)

In particular the full strain-gradient inertia has fixed, linear and
quadratic attachment-amplitude terms, not just the quadratic self term
displayed above. The positive norm C_rate^T P^-1 C_rate remains.
Isotropy annihilates D P^-1 C_rate as a map from STF rank2 to axial rank1,
so its cross qdot:Edot does not change leading b or g; this does NOT
annihilate the diagonal norm in(13). The same identity ensures the
isotropically averaged extra reaction spin is zero, preserving the mean
matching Dp=J qdot, although individual strain-driven cells have a nonzero
additional reaction. Reflection pairing likewise kills the mixed
grad(qdot)/strain-rate pseudotensor, while retaining both diagonal norms.
All these computed terms belong to M_U at order k² and to the physical
strain/spin current. The full norm(13), not a declaration of zero affine
KKS, is used in the final mass normalization.

For the separate neighbor-angle attachments, choose a finite carrier
with K_g/J_g>kappa/j plus a strict margin. Their additions satisfy

    Delta C_T,N=Delta C_L,N
        =t_g²[K_*- (kappa/j)J_*]>0.

This is the full added stiffness MINUS its added gradient-inertia cost.
The general normal-form expression is
C_N=C-2gb/rho-(kappa/j)(m_Phi-b²/rho), with b,g from the complete action;
here the completed leading b=0. Pick the gradient geometry first, retaining
its fixed finite macro corrections, then pick the STF amplitude. Finite
amplitudes give positive shear and both curvature eigenvalues. No old
material cage coefficient or unknown inverse is supplied as an input.

## 6. Actual observables, affine spin and unrestricted residual

U is the hybrid tube-centroid plus continuous-ambient mean defined by
0097/0098, not the Eulerian point mean. The exact phase-space observation is

    delta u_phase=V+curl(U cross u0)+Vop Ez,
    delta S_D=rho integral_D r cross delta u_phase
        +rho integral_boundaryD (eta.n) r cross u0.

Compact internal columns have no boundary term. Their actual material
spin discrepancy is G(E zdot)-L(Ez);(9) and the independent time-reversal
average make that retained internal discrepancy zero. The mean current
therefore retains the nonzero physical spin channel, not merely its KKS
label. The prescribed core-angle observation is the actual phase-space
vorticity-jet observation normalized in0085; every added response/cage
has disjoint core support.

The MACRO affine spin is also retained. For an invariant parcel its
time-reversal-even moving-boundary formula gives a geometric inertia
tensor I_D beta_dot plus the STF shape-spin row. The latter vanishes only
after the declared isotropic average. Hence, in a scalar isotropic notation,

    S_D,paired=S_aff[U]+j_cell n(n.qdot),
    S_aff[U]=I_D beta_dot + the retained nonisotropic shape row.

One does not infer S_D=j_cell Phi_dot pointwise unless that computed
macro row happens to match j_cell. The difference
(I_D-j_cell)beta_dot is an explicit affine spin/orbital current improvement,
and higher moments remain in the hybrid-to-point mean map. It affects
that observable map, not the leading b or ell already derived from the
same canonical action. The canonical spin conjugate to Phi is j Phi_dot;
the displayed formula states its relation to the full physical tube spin.

The complete reconstruction residual remains

    R=Udot-V+E zdot-A Ez.

No statement here declares it identically zero or a relabeling gauge.
Its mean and paired retained internal spin moments close as calculated;
other physical observations use0095's full residual/complement map if
an unrestricted Euler trajectory is asserted. The conditional canonical
Cauchy--Born pullback excludes nonaffine relaxation and retains this
reaction; that is the original premise, not an all-k invariance claim.

## Route verdict

Established conditional coefficient construction: the explicit six-moment
completion, full common-mean momentum reduction and independent fluid
reactions give b=0 and ell=-kappa/2, with actual one-action density and
positive shear/normal curvature supplied by the compact cages. The affine
physical spin, point-mean filter and unrestricted residual are retained
as separate calculated observations, not hidden behind canonical names.
No registry or release promotion is made in this child attempt.
