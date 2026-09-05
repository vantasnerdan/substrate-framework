# Complete material parcels: action, spin and finite-coherence estimates

This is a replacement for the point-mean representation refuted in 0072,
not a reinterpretation of that representation. All cells below contain
the same Euler fluid. In particular, cells without a selected EPS patch
still contribute their full mass. The construction proves an exact
reference-time material-coordinate bridge and a nonzero physical
parcel-spin response. Its time-homogeneous material reconstruction
obligation is stated separately below.

## 1. Exact complete partition and its action

At reference time choose the cubes

    C_j = R (j + [-1/2,1/2]^3), j in Z^3.

A uniform lattice shift makes this a stationary spatial ensemble; Haar
rotation of the entire field, lattice and marks makes it isotropic.
Select good EPS patches only at cube centres, with every generator and
response support contained in a ball of radius a < R/4. The existing
positive-probability bounded patch event supplies a positive selected
intensity p/R^3. Selection need not be independent between cells.
Unselected cells are retained, not deleted or assigned zero mass.

Let g_0(t) be the actual smooth Euler material flow and
D_j(t)=g_0(t,C_j). These sets are a complete material partition for every
time for which the smooth flow exists. They have fixed masses rho R^3;
their images need not remain cubes. Define X_j, V_j and r=x-X_j by their
actual material mass moments. For material momentum density pi, write

    g = X_j + h_j, integral_Cj rho h_j = 0,
    pi = rho P_j/M_j + eta_j, integral_Cj eta_j = 0.

Exactly, before any approximation,

    integral pi dot dg = sum_j [P_j dot dX_j + integral eta_j dot dh_j],
    integral |pi|^2/(2 rho)
        = sum_j [|P_j|^2/(2 M_j) + integral |eta_j|^2/(2 rho)].

The equations include the original incompressibility and shared-face
matching constraints. The two summands are an orthogonal material
decomposition, not two independently incompressible velocity fields.
Subtracting a piecewise constant mean need not preserve divergence or
interface traces separately, and is not asserted to do so.

In a rotating frame h_j=R_beta y_j the same cotangent form becomes

    P_j dot dX_j + integral eta_body,j dot dy_j
       + S_j dot (R_beta^T dR_beta)^vee,
    S_j = integral y_j cross eta_body,j.

This is the physical parcel angular-momentum map. The frame is a
coordinate choice, not an extra rigid rotor. Its angular momentum is
constrained to be this integral. Fixing a frame convention or imposing
an affine Cauchy--Born relation does not create a new mass term.

Every compact curl-generated Q or S from the selected patch has support
strictly inside its reference cube. Thus delta X_j=0 for every j on
these configuration tangents. The centroid cotangent two-form
sum dP_j wedge dX_j vanishes on any pair of these tangents. Consequently
subtracting the exact centroid kinetic Gram does not alter their
restricted KKS pairing D. Their induced *velocity* can nevertheless
change the momenta of remote parcels; those changes are retained below.

## 2. Exact pressure and interface bookkeeping

For each actual material parcel, Euler's equation gives

    M_j Vdot_j = -integral_boundary Dj p n,
    Sdot_j = -integral_boundary Dj r cross (p n),
    Tdot_j = -integral_boundary Dj p v dot n,
    (T_internal,j)dot = -integral_boundary Dj p (v-V_j) dot n.

The centroid term is V_j dot the first force. Paired shared faces have
opposite normals and the same physical p,v; the full pressure work
cancels there. The separate centroid and internal works generally do
not cancel independently. This is precisely the transfer that a
cell-level closure must retain. No exterior pressure flux is set to zero.

## 3. Exterior parcel means are small, not assumed zero

Use the 0069/0070 source with compact force F, zero zeroth moment and
first moment c I plus an antisymmetric matrix. Its exact Leray velocity
v has zero integral, finite angular moment and the exterior bound

    |v(x)| <= V_5 |x|^-5  for |x| >= a.

V_5 is an actual Green-kernel moment bound for the fixed smooth source;
it is not a fitted modulus. The compact curl and gradient parts in
0069 leave only the double-divergence potential outside the support.
This supplies the stated power and a finite V_5. On a bounded good-patch
event a common V_5 applies to every reaction column, including all
finite response corrections.

Put m_j(v)=integral_Cj v. Since integral v=0,

    |m_0(v)| <= 8 pi V_5 R^-2.

For j != 0, the elementary cube distance bound
|x| >= R(1+|j|)/16 gives

    |m_j(v)| <= 16^5 V_5 R^-2 (1+|j|)^-5.

Write a_0=8 pi, a_j=16^5(1+|j|)^-5 otherwise and A_1=sum a_j<infinity.
For arbitrary simultaneous reactions at selected lattice sites, Young's
inequality for this summable lattice kernel gives the FULL operator bound

    || C^* M C || <= rho V_5^2 A_1^2 R^-7 = delta_R.

Here C maps reaction amplitudes to all parcel mean velocities and M is
the diagonal parcel mass operator. Either counting measure divided by
total volume, or its equivalent stationary/Palm norm, gives this same
bound. This is not an independent-cell inverse and does not omit
velocity cross terms. Missing selected sites are zero amplitudes in the
same inequality.

If P is the full reaction energy operator, the exact centred operator is

    P_c = P - C^* M C.

Thus P_c >= (p_min-delta_R) I when P >= p_min I. Choosing the geometrical
coherence radius R sufficiently large preserves strict positivity while
counting all ambient fluid mass. The subtraction is the exact material
centroid split; it is not an infrared subtraction from Euler energy.

## 4. The actual parcel-spin row is nonzero

Let D_ang(v)=rho integral x cross v be the source angular moment and

    A(v) = rho sum_j integral_Cj (x-Rj) cross v

the sum of actual reference-parcel intrinsic spins. Their exact difference is

    A(v)-D_ang(v) = -rho sum_j Rj cross m_j(v).

Both sums converge absolutely. With
A_ang=sum_{j!=0}|j|a_j<infinity, the previous bound gives

    |A(v)-D_ang(v)| <= rho V_5 A_ang/R.

The proof is uniform for every source in the selected good-patch class.
For example, bounding lattice shells by 24 n^2+2 proves both A_1 and
A_ang finite using respectively sum n^-3 and sum n^-2. Nothing relies on
decorrelation. In 0070 the reaction source has D_ang=B n with B bounded
away from zero. Its physical parcel-spin row therefore remains nonzero
for an explicit finite R. This differs from the zero-cell-size
point-mean representation: here the moment arm of the spin is the actual
finite parcel, while ALL induced remote centroid momenta are retained.

The variation of the physical parcel spin is indeed this velocity
integral at the reference state: each configuration generator vanishes
near every face, so its moving-boundary and delta X contributions vanish.
For noncompact affine generators the corresponding boundary terms stay
in the rotating-frame action of section 1; they are not silently discarded.

## 5. Physical-spin kinetic normal form, including locked affine inertia

Suppose the centred shape/reaction slice is reconstructed from the same
material cotangent action, so its actual connection is

    <s, D qdot + A betadot>,

where D is the retained KKS row and A the *parcel*, not global, spin row.
This is the remaining reconstruction condition, not a definition of A.
Eliminate the full reaction operator P_c, retaining all affine kinetic
terms. After the prescribed isotropic action average, write

    T_rot = (a qdot^2 + 2 b qdot betadot + c betadot^2)/2,
    a = D^* P_c^-1 D,
    b = D^* P_c^-1 A,
    c = A^* P_c^-1 A + c_locked.

c_locked includes genuine affine within-parcel mass and all other
retained material modes; it is not set to zero. The scalar notation is
the isotropic axial coefficient. The underlying operators are retained
until this stage. With c_locked >= 0, Cauchy--Schwarz gives c-b^2/a >= 0.
The bound A-D=O(R^-1) and coercivity imply

    |b-a| <= ||D|| ||P_c^-1|| ||A-D||.

Since a >= ||D||^2/||P_c|| in each nonzero isotropic direction, a finite
R makes b nonzero. Define the common-rotation-covariant collective angle

    Phi = beta + (a/b) q,  j = b^2/a.

Then, exactly,

    T_rot = j Phidot^2/2 + (c-j) betadot^2/2,
    S_physical = partial T_rot/partial betadot at fixed qdot
               = j Phidot + (c-j) betadot.

The second identity uses the physical rotation momentum map, not a
canonical name for spin. It retains the ambient/affine spin identified
by the parent audit. Since beta=curl U_centres/2, the residual is a real
positive macro gradient inertia. A potential K q^2/2 becomes
K(b/a)^2 |Phi-beta|^2/2, without fitting K. This transformation preserves
independent material core rotation, but Phi is the specified collective
mass-spin angle: the literal core angle beta+q equals
beta+(b/a)(Phi-beta), which must remain in the observation map.

These identities establish how a true finite-parcel embedding can avoid
0072. They do NOT assume that a KKS restriction alone already proves the
displayed material connection. In particular, the non-rotating shape
coordinates, pressure constraints and reconstruction of the orbit
velocity into a material path must be pulled back together.

## 6. Exact scope and next executable continuation

Established: a complete physical reference partition; exact full-density
mass/cotangent/pressure bookkeeping; full-operator centroid correction;
nonzero physical parcel-spin row at finite coherence radius; and the
kinetic normal-form identity conditional on the actual reconstructed
shape connection. The fixed reference partition has a genuine material
continuation D_j(t)=g_0(t,C_j), not fixed Eulerian boxes.

The Gaussian stationary velocity law does not imply that this advected
partition has a stationary bounded-shape law. Its material inertia and
connection may evolve. Thus this attempt has not established a
time-homogeneous stationary material-parcel constitutive action, nor
identified the exact centred shape connection with D qdot+A betadot
for all the constrained orbit paths. These are named construction
obligations, not new all-wavelength or unrestricted-invariance demands.

Next executable route: pull the known compact Q/S orbit path back into
the exact finite-parcel cotangent coordinates (section 1), retaining its
material reconstruction and affine shape modes before eliminating any
momentum. Compute the resulting connection and its pressure/interface
terms. In parallel, test a tube-plus-ambient material decomposition:
the selected invariant EPS tubes have stationary boundaries, while the
ambient material mass remains a continuous phase rather than being
forced into bounded stationary parcels. Its coherent translation and
relative spin are then to be derived from the same action, not assigned.
No exhaustion or parent-completion verdict follows from this record.
