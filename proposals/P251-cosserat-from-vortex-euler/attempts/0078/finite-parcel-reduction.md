# Finite material parcel: exact reduction and a positive fixed-Kelvin spin

## 1. The mass centroid and its cotangent momentum separate exactly

Let D0 be a finite reference parcel, dm=rho0 da, M=integral_D0 dm, and g its
actual incompressible material map. All neighboring maps and pressure tractions
remain in the global action; considering its restriction to D0 does not
introduce an independently movable internal boundary. Define

    X=M^-1 integral g dm,       r=g-X,       integral r dm=0.

For material momentum per unit mass v, put P=integral v dm and
pi=v-P/M. The exact kinetic energy and canonical one-form are

    T=|P|²/(2M)+1/2 integral |pi|² dm,
    Theta=P.dX+integral pi.dr dm,
    integral pi dm=0.

Eliminating P gives M|Xdot|²/2 with no internal kinetic cross. Unlike the
pointwise ensemble mean in 0072, X is the mass centroid of a finite set of
material particles. All within-parcel affine velocities remain in pi.

Intrinsic angular momentum is the actual mechanical quantity

    S=integral r cross pi dm=integral r cross gdot dm.

The centered canonical bracket gives {X_i,P_j}=delta_ij,
{S_i,S_j}=epsilon_ijk S_k, and {S_i,X_j}={S_i,P_j}=0. Its action on centered
positions is the ordinary infinitesimal rotation. These relations follow
from the canonical material bracket, not from assigning a spin algebra to a
vortex center. Boundary pressure gives the usual force and torque flux:
Pdot=-integral_boundary p n dA and Sdot=-integral_boundary r cross p n dA,
with any prescribed external tractions included in the same expression.

## 2. An explicit frame chart and the geometric mechanical connection

Choose a centered reference template c(a). In the open chart where
C=integral r c^T dm is invertible with positive determinant, its right polar
decomposition C=Q U, U symmetric positive, gives an SO(3)-equivariant frame:
Q(Rr)=R Q(r). Put h=Q^T r. The shape satisfies

    integral h dm=0,    integral h c^T dm symmetric positive.

This is an explicit polar/Eckart gauge, not an assumption that a fluid parcel
is rigid. With Q^T Qdot=[Omega cross], the exact centered kinetic identity is

    T_center=Omega.I(h)Omega/2+Omega.l(h,hdot)+||hdot||²/2,
    I(h)=integral (|h|² I-h h^T) dm,
    l(h,hdot)=integral h cross hdot dm.

For a genuine three-dimensional parcel, I(h)>0. Its entries are computed
mass moments. The body angular momentum is

    S_body=Q^T S=I(h)Omega+l(h,hdot).

Equivalently, with mechanical connection A_h(hdot)=I^-1 l,

    T_center=(Omega+A_h hdot).I.(Omega+A_h hdot)/2
                +[||hdot||²-l.I^-1 l]/2.

The last quadratic form is the shape metric orthogonal to physical rigid
rotation. Its nonnegative sign follows by minimizing the original positive
norm over Omega. Connection curvature need not vanish; it describes actual
angular momentum exchange between shape and frame and is retained.

For finitely selected shape rates zdot, write l=B zdot and
||hdot||²=zdot^T G zdot. If their momenta pi_z are genuinely cyclic/conserved
under the stated constraints, their exact Routh elimination gives

    I_Routh=I-B G^-1 B^T,
    L_Routh=Omega.I_Routh.Omega/2
             +Omega.B G^-1 pi_z-pi_z^T G^-1 pi_z/2.

I_Routh is positive when the retained frame and shape velocities form an
independent positive metric block. An infinite-dimensional version uses the
appropriate bounded coercive shape operator and closed subspaces. These
are exact Schur statements about the computed material metric. They do not
say that arbitrary fluid shape momenta are cyclic or that fixing them equals
fixing all Kelvin circulations.

## 3. Kelvin reduction: what is conserved and what can be gauge

The right volume-preserving relabeling action g->g composed with eta has
momentum map represented by the material one-form

    [g^*(v flat)] modulo exact one-forms.

Fixing it retains every material-loop circulation, including harmonic
periods on a multiply connected parcel. A local curl equation alone would
not fix those periods. At a prescribed g and prescribed normal boundary
velocity, incompressibility plus this circulation data constrains v. One
cannot independently append Omega cross r to that v: it adds vorticity
2Omega, the exact failure exposed in 0042.

Spatial translations/rotations act on the LEFT and commute with right
relabeling. Their canonical momentum maps P,S therefore descend through
Kelvin reduction wherever the spatial action preserves the actual domain
and its retained boundary/traction constraints. Fixing Kelvin data does not
by itself erase physical spatial angular momentum.

However the polar frame above uses material labels and need not descend.
For a spherical configuration g=id, a spatial rotation of its positions is
also the right relabeling eta(a)=Ra. A reference-template frame rotates,
while the unlabeled spatial density remains the same ball. Thus a positive
locked metric I_geom does not by itself establish a physical fluid rotor
after relabeling. Adding that metric to an already reduced Euler-orbit
action can count a direction twice or mix circulation ensembles.

The repair is to use a frame observable of the REDUCED Euler state, such as
the orientation of a nondegenerate vorticity/core jet or its anisotropic
tensor. It is relabeling invariant, since it is a function of the actual
Eulerian vorticity and registered parcel geometry. For rotations about e,
choose a local angle q with q(R_theta omega)=q(omega)+theta. Then on the
fixed-Kelvin orbit its physical angular momentum J=e.S satisfies
`{q,J}=1`. This identifies the physical canonical angle; its inertia still
comes from the reduced Hamiltonian, not from the existence of this bracket.

## 4. Positive fixed-Kelvin absolute-spin theorem

Suppose the full finite-parcel action and its boundary data have a genuine
SO(2) rotation symmetry about the centroid, and the stationary Euler state
is not invariant under that rotation. Let K=e cross r be its generator on
the actual Kelvin-reduced state. Assume the physical angle observable q is
equivariant as above. For any admissible same-leaf direction eta,

    dJ(eta)=Omega(K,eta).

Let the stationary reduced Hamiltonian Hessian be H. Its symmetry and
stationarity imply H(K,K)=H(K,eta)=0: differentiate the exact rotation
invariance of H_energy, using its vanishing first variation on the leaf.
Suppose an admissible shape direction eta, with dq(eta)=0, has

    B=dJ(eta) !=0,        h=H(eta,eta)>0.

The selected two-direction SAME Euler action is then

    L=B p qdot-h p²/2,
    delta J=B p,
    L_reduced=I_red qdot²/2,    I_red=B²/h>0.

Thus positive absolute-spin inertia CAN coexist with fixed Kelvin
circulation. Its physical angular momentum is I_red qdot. It is a Hessian/
symplectic Schur result on that fixed leaf and generally differs from
I_geom=e.I(h)e. No geometric inertia is added. Additional retained shape
directions require the full block inverse, I_red=D*P^-1D, and any magnetic
or mixed energy terms must remain before reduction.

This theorem is transferable when its symmetry and direction hypotheses
are satisfied. It does not infer parcel rotation symmetry from isotropy of
a probability law, and it does not remove torques from neighboring fluid.
For a nonsymmetric moving parcel one instead retains its actual boundary
torque and noncyclic angle energy in the same Routh calculation.

## 5. A concrete smooth finite Euler domain satisfying the positive theorem

The existence conditions above are not empty. Consider an incompressible
Euler fluid in a ball of radius R with the material/slip boundary condition
u.n=0. Pressure traction is normal to the sphere and exerts exactly zero
torque. This is a stated finite-domain/shape constraint, not a claim that an
arbitrary internal EPS parcel can move independently of its neighbors.

Let lambda>0 and choose R at a positive zero of the spherical Bessel j_1.
Set f(r)=j_1(lambda r)/r, continued smoothly at r=0, and define

    psi=z f(r),  V=curl(psi r_vector)=f(r) e_z cross r_vector,
    u=V+curl V/lambda.

The radial identity f''+4f'/r+lambda²f=0 implies curl u=lambda u, because
V is divergence free and solves its vector Helmholtz equation. Explicitly,

    u=f e_z cross r
       +[(2f+r f')e_z-(f'/r)z r_vector]/lambda,
    u.r_vector=2f z/lambda.

Hence u.n=0 on the chosen sphere, and u is an exact smooth stationary Euler
field with p=-rho|u|²/2+constant. At the center u=(2/3)e_z; rotating about
e_x changes its actual vorticity direction. A physical local q can be read
as atan2(-omega_y(0),omega_z(0)), so dq(K)=1. This is not a template-label
angle. The geometric ball inertia is positive but is not the inertia used
in the following construction.

For this constant-lambda field, put v_K=P_D(K cross omega), where P_D is the
finite-domain Leray projector with zero normal velocity. Its Eulerian value
is the actual derivative of the rotated field,
v_K=e cross u-(K.grad)u, and curl v_K=lambda v_K. It is nonzero at the core
and, by analyticity, on some off-core open ball. Choose chi supported there,
away from the physical angle observation, and set

    eta0=curl(chi curl v_K),
    B0=Omega(K,eta0)=-rho integral chi |curl v_K|² !=0.

This is a compact same-Kelvin direction with dq(eta0)=0. Add a disjoint
negative-helicity compact high-carrier cage A_k in a flow box with nonzero
omega. Define

    eta_k=(1-Omega(K,A_k)/B0) eta0+A_k.

Then dJ(eta_k)=B0 exactly and its physical core-angle variation remains zero.
The finite-carrier estimates of 0045 apply on the ball directly: the compact
curl/gradient parametrix is formed before P_D; P_D fixes the compact curl,
kills the gradient and is an L² contraction. The principal Hessian is
positive of order k/lambda, while the fixed response and its cross terms
remain bounded. Thus H(eta_k,eta_k)>0 at a finite carrier.

The Hessian formula used here is

    H(xi,eta)=rho integral [v_xi.v_eta-v_xi.curl v_eta/lambda].

Its curl integration boundary form is zero on the specified leaf: compact
vorticity rearrangements have zero boundary normal vorticity, so the
tangential velocity one-form is closed; on the boundary sphere it is exact.
The same holds for v_K. The integral of the wedge of two such exact boundary
one-forms vanishes. On multiply connected boundaries the corresponding
periods would need to remain fixed explicitly rather than using the sphere
argument unchanged.

This supplies all B!=0, h>0 and physical-angle hypotheses, giving a concrete
positive absolute-spin action on one finite smooth Euler Kelvin leaf. It
does not yet supply a knotted EPS ensemble, an independently compatible
ambient extension of all trial states, or the parent parcel assembly. Those
are separate consumers of the license, not part of this example's verdict.

## Route result

Established: exact finite-parcel centroid/cotangent and frame/connection
identities; the precise distinction between a label-dependent locked metric
and a physical reduced inertia; and a constructive positive fixed-Kelvin
absolute-spin theorem with a smooth finite-domain example. The ball shape
constraint and its no-torque pressure boundary are explicit. The parent
0075 assembly must supply its actual partition/ambient compatibility rather
than importing that boundary constraint silently.
