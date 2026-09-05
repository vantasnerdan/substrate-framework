# Actual local Kelvin blocks and the finite joint-form normalizer

This transfers the algebra of0228 by constructing its physical input blocks
on a different stationary Euler background. It does not transfer the C016
two-scalar sector itself. The final periodic-background version is conditional
on the actual compact-ring assembly being supplied by the geometry route.
Its local input fields and forms are actual fields, not prescribed matrices.

## 1. Exact full forms require only local constant curl

Let u be a smooth stationary incompressible Euler field and let
omega=curl u=lambda u on an open patch, lambda nonzero constant. Choose
smooth solenoidal generators xi,eta compactly supported in that patch.
Their actual initial velocities are v_xi=P(xi cross omega), with the full
Helmholtz projector. Pressure tails are retained. On those supported
generators the exact Kelvin/Lin operator is

    L eta=-[u,eta]+v_eta=v_eta-curl(v_eta)/lambda.         (1)

Indeed curl(v_eta)=curl(eta cross omega)=[omega,eta]
=lambda[u,eta] on the source support and both bracket sides vanish outside.
For any divergence-free zeta, integration by parts/projection gives

    Omega(xi,zeta)=rho integral omega dot(xi cross zeta)
                 =-rho integral v_xi dot zeta.           (2)

Consequently the full stationary Kelvin energy Hessian is exactly

    H(xi,eta)=-Omega(xi,L eta)
       =rho integral [v_xi dot v_eta
                        -v_xi dot curl(v_eta)/lambda].   (3)

This is the same physical KKS/Hamiltonian convention as the canonical
coadjoint module and0221. It includes the global kinetic pressure tail.
Equation(1), rather than a global constant-curl assertion, licenses(3).
The stationary Jacobi form restricted to the actual Kelvin/Lin initial data
has this Hamiltonian Hessian; it is not the unconstrained static material
stiffness obtained by setting xi_dot=0. Thus the complete pressure/Jacobi
and coadjoint conventions remain distinct and correctly connected.

The identities apply at zero wave number on the fixed0211 ring with its
actual uniform far velocity: xi is compactly supported, v is square
integrable, and the pairings(2)-(3) are finite. No infinite background energy
is assigned a finite value. For Bloch jets below, take an actual periodic
background and the periodic Helmholtz projector; that global hypothesis is
separate from the local calculation.

## 2. Opposite actual signatures from helical curl potentials

Choose a small contractible patch in an invariant source annulus disjoint
from the observed tag. Pick a constant nonzero covector k with k.u>0
throughout the patch. Let e1,e2,k_hat be a positively oriented orthonormal
frame and p_sigma=e1+i sigma e2, sigma=+1 or-1. For a real compactly
supported smooth profile G, set

    eta_sigma,N = sigma G p_sigma exp(i N k.x)/(N |k|),
    xi_sigma,N = curl eta_sigma,N.                        (4)

These are exact smooth compactly supported divergence-free generators;
their real/imaginary parts supply two real columns. Their leading amplitude
is G p_sigma because i k cross p_sigma=sigma |k|p_sigma.
The actual projected velocity has leading amplitude

    P_k(p_sigma cross omega)G
       =i sigma (omega.k_hat) G p_sigma.                 (5)

Write w=omega.k_hat, with its sign fixed on the patch. For the real pair,
the full forms derived from(2)-(5) have

    Omega_sigma= sigma P_G J_2+O(N^-1),
    H_sigma= -sigma N E_G/lambda I_2+O(1),               (6)
    P_G=rho integral w G^2,
    E_G=rho integral |k| w^2 G^2,
    J_2=[[0,1],[-1,0]].

The remainders are bounds for actual forms, uniform over any fixed
finite-dimensional compact profile family with bounded smooth norms.
They follow by the full-projector symbol expansion applied to(4);
background derivatives, envelope derivatives and all pressure terms are
lower order. In particular both actual energy signatures are definite
for sufficiently large finite N. They are not numerical eigenvalues or
negative material probabilities.

Whiten the actual definite2x2 forms by orientation-preserving real matrices,
so the positive block has H=I and the negative block H=-I. Their phase
scalars have opposite signs and magnitudes

    |kappa_sigma|=|lambda P_G|/(N E_G)[1+O(N^-1)].        (7)

Since lambda w= lambda^2(k.u)/|k|, the positive ratio in(7) lies between
the reciprocals of the extrema of k.u, divided by N. Choose two fixed
phase covectors with distinct scales and sufficiently small source patches
so these reciprocal intervals are disjoint. This supplies

    |kappa_positive+kappa_negative| >= c/N >0            (8)

for all large finite N, with either sign available by swapping the phase
scales assigned to the two signatures. A nonconstant period law is not
needed for these normalization controls; they are not advertised as the
slow physical clock columns.

## 3. Complete finite constraints and actual block construction

Every cross KKS/H row against a previously fixed source is linear in G.
Every required finite spatial derivative of such a row is also linear.
Choose M+1 fixed disjoint smooth subprofiles inside the selected patch for
M real homogeneous constraints, and take a unit vector in their kernel.
The disjoint supports preserve positive upper/lower bounds for the weighted
integrals P_G,E_G and uniform derivative bounds. Thus all signature and
ratio margins survive the constraints. The true nonlocal pressure pairing
is included in the constraints, not set to zero from disjoint supports.

Construct a positive and negative pair, imposing all cross rows between
them and the baseline sources. Whiten the actual forms. Add the two pairs
with equal amplitude a and choose a^2=b/(kappa_positive+kappa_negative),
using the sign assignment making that quotient positive. The resulting
two actual columns have exactly

    H=0,              Omega=b J_2.                       (9)

Each individual energy is retained until its exact cancellation. The source
amplitude is polynomial in N: whitening is O(N^-1/2), a is O(N^1/2),
so the phase-normalized generators remain O(1) at derivative order zero
and grow polynomially at each fixed higher order.

A single generator from either definite block has zero self KKS and either
sign of full energy. Scaling supplies exactly H=+1 or-1. With finitely
many additional cross constraints this constructs the actual P and E_sigma
blocks required by0228, mutually cross-orthogonal through the chosen jet
order and orthogonal to the finite baseline columns. No abstract supplied
Gram matrix has been mistaken for an actual Euler field.

## 4. Full Bloch jet on an actual periodic cell

This section assumes the background is an actual smooth periodic stationary
Euler field containing the above local constant-curl patches and invariant
tag/source separation. The future disjoint compact-velocity assembly would
supply it exactly. The isolated0211 ring does not supply it by itself.

Use the compact potential in(4) to define

    xi_K=curl_K eta_sigma,N,
    v_K=P_K(xi_K cross omega),     curl_K=curl+iK cross.  (10)

The full Bloch divergence vanishes identically. The supports remain inside
the constant-curl patch. Include among the finite homogeneous profile rows
every component of the mean of xi_K cross omega: it is a polynomial of
degree at most one in K, so these are finitely many exact linear constraints.
The initial Kelvin force then has zero mean for all K. This removes the
singular initial mean projector without deleting the actual ray-wise mean
of the evolved Euler field.

Polarization of the real +/-K pair, or the local identity(1)-(3), gives the
full complex forms

    Omega_K=rho integral omega dot(conjugate(xi_i,K) cross xi_j,K),
    H_K=rho integral conjugate(v_i,K) dot
                         [v_j,K-curl_K(v_j,K)/lambda].   (11)

The potential normalization in(4), the full projector and their finite K
derivatives have the usual joint covector-parameter orders. For the
unnormalized O(1) generators, phase order is0 and energy order1 in N;
each K derivative lowers the order by one. All remaining low Fourier rows
are smoothing terms. Thus the actual higher jets and cubic remainders have
finite polynomial costs. Exact jet matching does not require them to vanish.

Apply0228's triangular finite algebra to these actual blocks: correct H0
and Omega0, then energy1 before phase1, then both degree2 forms with
amplitude rows vanishing at K=0. Retain the intrinsic higher jets after
every step. All cross constraints were imposed on fixed profiles for the
whole jet, so no K-dependent uncontrolled kernel choice occurs. The
construction matches the prescribed finite joint physical phase/energy jet
exactly through degree2 on this same periodic field.

## 5. Observation errors and the preparation diagonal

The source patches lie in an invariant annulus with a positive separation
from every observed material tag. The full Euler/Lin propagator has scalar
transport principal symbol plus order-zero pressure/matrix terms. The
arbitrary finite-order transported parametrix in0221 applies to this smooth
background without its special C016 two-scalar sector. Cutting its source
and observation by the invariant separation makes its kernel smoothing.
The compact oscillatory source in(4) has arbitrarily small negative Sobolev
norms after retaining its actual polynomial normalization. Consequently
every required finite remote observation and its time/K derivatives is
bounded by C_q N^-q for any fixed q. Whole-fluid mean/current rows use the
corresponding smooth propagated adjoint tests and have the same bound.
The pressure tail is small and generally nonzero.

First freeze the finite baseline physical columns and their complete forms.
Choose the finite normalizer profiles, constraints and source costs, then
the high-carrier/long-wave diagonal exactly as in0228. Its proof needs
polynomial costs for that finite physical supplier inventory; an arbitrary
dense-range statement without those costs is insufficient. The new
high-carrier fixed-tag suppliers have explicitly retained costs and need
a combined sparse-carrier ordering with these control costs. That final
ordering must be written for the actual selected inventory, rather than
declared automatic from density.

## Result and actual remaining joins

The local opposite-sign physical Kelvin blocks, finite exact cross
constraints, and zero-wave phase/energy normalizer are established on the
fixed ring at the stated patch/gap scope. The complete degree2 normalizer
construction is established on an actual periodic field satisfying those
same local hypotheses, conditional on that field being supplied. The
stationary Jacobi/KKS forms include the full Leray pressure and use the
same physical generators throughout.

No full current conclusion is inferred solely from matrix matching. The
remaining actual joins are the selected fixed-tag supplier/control error
diagonal and the0246 full material Ward/current bridge on those same
sources. The geometry route still owes its compact stationary field and
periodic positive-density assembly. This result identifies the actual
transferable normalizer and its global hypotheses; it does not complete
the parent campaign or promote a new claim by itself.
