# Compact Euler-orbit attachments, full reactions, and hybrid mean

## 1. Build zero-spin compact canonical cages on the same EPS field

Use0085's analytic compact-velocity operator and corrected WKB pair C_k,S_k.
Both xi and v=xi cross omega are compact in the same invariant tube D,
div xi=div v=0, and the exact SAME-Euler forms are

    H(xi,eta)=rho integral[v_xi.v_eta-v_xi.curl(v_eta)/lambda],
    Omega(xi,eta)=rho integral omega.(xi cross eta).

These are coadjoint forms, not the material-Jacobi functional of0090--0094.
Choose the open positive symbol ball of0085 and retain all subprincipal
returns. Its finite analytic estimates give

    H_pair=h0 I+O(k^-1), h0>0,
    B_k=Omega(C_k,S_k)=B0/k+O(k^-2), B0!=0.

The physical spin row L(xi)=rho integral r cross v_xi is O(k^-2), with
arbitrarily higher decay available by compact oscillatory integration.
Choose0085's three normalized compact spin responses eta_i in mutually
disjoint off-core balls, also disjoint from the raw pair, and put

    Q_g=C_k-sum eta_i L_i(C_k),
    S_g=S_k-sum eta_i L_i(S_k).

Now L(Q_g)=L(S_g)=0 EXACTLY. The KKS correction is exactly zero, including
the response-response part: distinct eta supports are disjoint, and each
Omega(eta_i,eta_i)=0. Their finite H correction is O(k^-4). Thus the full
2x2 Hessian remains positive and B_k remains the same nonzero number.
Write its entries as Hq,N,P, P>0, and define

    J_g=B_k²/P,
    K_g=Hq-N²/P,
    R_g=K_g/J_g=det(H_g)/B_k².                 (1)

This is the complete momentum reduction; neither Hq alone nor an assigned
inertia appears. For finite k beyond the analytic remainder thresholds,

    H_g >= (h0/2) I, |B_k|<=2|B0|/k,
    R_g >= h0² k²/(16 B0²).                   (2)

Consequently R_g exceeds any fixed positive base optical ratio r0=kappa/j
by a strict finite margin. Unlike0093, this mechanism uses B_k~1/k in the
compact COADJOINT sector. The material Rayleigh quotient is not imported.

## 2. Exact physical mean and current of these reactions

Compact divergence-free velocity gives

    integral_D v=0,
    integral_D (r_i v_j+r_j v_i)=0.

For Q_g,S_g the zero spin row also kills the antisymmetric part. Thus

    integral_D r_i v_j=0                         (3)

for BOTH generators. The ambient induced velocity is identically zero;
no harmonic pressure tail crosses the tube boundary. Their boundary normal
source in the linear material tag equation is exactly zero. Therefore the
hybrid momentum (tube centroid momenta plus continuous ambient point
momentum) of each reaction is exactly zero. This is an actual material
mean statement, not identifying the hybrid mean with the Eulerian point
mean. The latter is obtained by the exact compact Fourier velocity profile.

For a slowly modulated amplitude a, the Palm coherent point-velocity mean
has profile integral exp(-ik.r)v(r)dr=O(k²), by(3) and the zero zeroth
moment. A gradient amplitude a=O(k) therefore produces only O(k³) mean
return. Its full spin density can have higher multipole filters, retained
by that same Fourier integral; zero integrated spin does not imply zero
pointwise spin density. The actual base angle pair of0085 instead has
L(Q)=0,L(S)=B n and supplies the nonzero curl-of-spin current in the
point-mean/hybrid-mean map. These two populations serve different roles.

## 3. Full Routh/Dirac block, including macro forcing

Freeze a finite set of compact cells and all their response profiles.
Distinct cells have disjoint xi AND v supports, so their H and KKS cross
blocks are exactly zero. The inverse below is consequently the inverse
of the ACTUAL complete retained reaction block, not an assumed isolated
cell replacement for a nonlocal Leray operator. In a stationary ensemble
this identity is taken before averaging, with independent reaction
coordinates on different cells and on time-reversed realizations.

Let e denote retained macro affine/gradient data, a(e) the declared attached
coordinate, and s the vector of all its fluid reactions. The full local
first-order action has the form

    L = L_macro + s^T B adot
        -[s^T P s+2s^T(N a+C e)
          +a^T H a+2a^T F e+e^T A e]/2.

B may be diagonal by the disjoint supports, but the formula does not need
that choice. P is the computed positive reaction matrix. Eliminating s gives

    L_red=L_macro
       +adot^T B^T P^-1 B adot/2
       -adot^T B^T P^-1(N a+C e)
       -[a^T H a+2a^T F e+e^T A e
          -(N a+C e)^T P^-1(N a+C e)]/2.     (4)

The C e term is important: freely varying a new reaction can respond to
macro strain even when its attached coordinate vanishes. Its negative
potential correction is not omitted. On u and -u, P,H,N,C,F,A are unchanged
whereas B changes sign; eliminate their INDEPENDENT s variables first.
The middle gyroscopic term cancels, while the full positive inertia and
full potential Schur correction remain. Averaging B to zero before this
elimination would incorrectly erase inertia.

The pairing of a tracefree symmetric affine generator with any compact v
vanishes by its symmetric first moment. The zero-spin cages also have zero
rotational affine pairing. Thus there is no missing affine KKS forcing of
s in(4); the derivative adot comes from the genuinely varied attached
coordinate. Potential macro cross blocks C,F remain as shown.

This is precisely a constrained-coordinate action/Dirac reduction: a(e)
is imposed on the coordinate family, then independent conjugate fluid
momenta are varied. Fixing s as well would be a different constraint and
would not give(4). No unrestricted nonaffine relaxation is added to the
original Cauchy--Born family. The remaining shared material/pressure lift
is0097's input, not replaced by a declaration that xi is the actual
material displacement at every time.

## 4. Positive STF shear from the complete square

Choose an independent cage population with scalar a=t_s T:E, where
E=sym grad U and T is a unit tracefree symmetric tensor rotated with the
complete background/cage configuration. Its SO(3) average is

    average (T:E)²=|dev E|²/5.

The t_s² coefficient in the FULL potential(4) is K_g>0, not Hq:

    mu(t_s)=mu_fixed+b_s t_s+a_s t_s²,
    a_s=nu_s average K_g/10>0.                 (5)

Here mu_fixed includes ALL A-C^T P^-1 C contributions of the newly
retained reactions and the actual macro background; b_s includes
F-N^T P^-1 C. The exact signed quantities are Euler action integrals.
Any finite fixed deficit is overcome at a finite t_s. Their same-action
kinetic self term is

    Delta T=nu_s t_s² average J_g |dev Edot|²/10.

There is no imported material shear or appended mass. This U-gradient
mass is retained even though it first affects acoustic dispersion at k⁴.
Selecting t_s by structural positivity is not fitting an elastic datum.

## 5. Positive spin curvature after its COMPLETE added inertia

Use a separate support-disjoint zero-spin cage population with coordinate

    a=t_g [n.q(X+h/2)-n.q(X-h/2)].

Use three spanning equal-length bonds and isotropic axis marks, rotated
with the complete field; reflection pairing removes polar/axial odd
crosses. Its exact difference multiplier is 2i sin(k.h/2). From(4),

    Delta M_PhiPhi=t_g² J_* |k|² I,
    Delta K_PhiPhi=t_g² K_* |k|² I,
    K_*-r0 J_*>0, r0=kappa/j.                 (6)

The last inequality follows by choosing the finite carriers in(2) with a
uniform strict margin, then averaging the positive quadratic forms. It
does not replace average(B²/P) by (average B)²/(average P). The current
moments(3) show that no new mean correction enters through degree two.
Zero core support leaves the physical core-angle map unchanged. Replacing
q by Phi-curl U/2 produces genuine degree-three mixed and degree-four
translation terms, retained in the exact discrete action.

The C e macro forcing in(4) gives a fixed strain correction and a possible
term linear in t_g coupling E to grad q. The latter is parity odd and
vanishes in the complete reflection-paired isotropic action. The fixed
strain correction is included before choosing the positive shear amplitude
in(5). Common rigid rotation has zero energy cross with compact variations
by Euclidean covariance of the full orbit Hessian. The existing base
angle profiles are disjoint, so the base j,kappa and its leading mixed
kinetic/potential coefficients are unchanged.

For the full computed transverse jet with coefficients b,g,m_Phi,C,

    C_N=C-2gb/rho-(kappa/j)(m_Phi-b²/rho),
    C_L,N=C_L-(kappa/j)m_Phi,L.

Equation(6) increases BOTH by t_g²(K_*-r0 J_*)>0. A finite amplitude therefore
overcomes all fixed normal-form curvature corrections. The added mass
remains in Phi_phys=[1-(m_Phi-b²/rho)k²/(2j)]Phi_N and in its longitudinal
counterpart. Dropping it would change the optical k² coefficient.

## 6. Where the macro mass comes from, and the exact joining interface

For the SAME material Euler configuration on a tube parcel, write
g=X_D+r, integral_D rho r=0. The exact kinetic/cotangent decomposition is

    T_D=M_D |Xdot_D|²/2+integral_D rho |rdot|²/2,
    Theta_D=P_D.dX_D+integral_D pi_r.dr,
    integral_D pi_r=0.

This precedes Kelvin reduction;0078 proves the spatial/relabeling momentum
maps commute. Here both compact orbit directions have zero actual centroid
velocity and boundary normal trace, so their local response belongs to the
centred factor and cannot consume or duplicate M_D. The continuous ambient
is retained with its own original Euler kinetic integral, not replaced by
finite imaginary ambient parcels. Under the shared slow translation, the
tube volume weights plus ambient indicator exhaust the fluid volume, so
the uniform translation kinetic density is exactly rho |Udot|²/2. This
can also be checked by a uniform Galilean variation of the homogeneous
mean-zero Euler ensemble BEFORE subtracting a relative energy; no finite
mass is assigned to an infinite global boost.

The compact reaction moment is exactly zero in that hybrid cotangent row.
Its point-current difference and higher moment expansion were computed in
section2. At nonconstant U, the actual common-affine lift, its pressure
reactions, and its remaining finite-gradient mass/current terms must still
be those of the joint material chart in0097. Equations(4)--(6) retain their
computed coefficients; they do not derive that chart by merely naming
the coadjoint generators material displacements. Thus this construction
supplies positive shear and BOTH positive normal-form curvatures with all
compact reaction and mean terms, conditional on that explicit shared lift.

## Route verdict

Established: compact coadjoint zero-spin cages and finite structural
selection give positive full-Schur STF shear and positive mass-normalized
spin gradients on the same actual EPS Euler background. Compactness makes
the reaction locality exact, and the physical hybrid mean and point-current
moments are explicitly identified. No material-Jacobi coefficient, isolated
unlicensed inverse, rotor mass, empirical fit, or all-k closure is used.
The parent mean/affine Kelvin chart is the separately constructed0097
joining input; this proof does not claim to have completed that input.
