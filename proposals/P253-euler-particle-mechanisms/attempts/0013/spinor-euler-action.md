# Exact spinor coordinates for Euler, with their actual Hamiltonian

This is the failure-derived multicomponent continuation registered in0013.
The following calculations are independent derivations. They use local regular
Clebsch coverage, not an unproved global representation of every Euler field.
The full orbit/Hodge action of0005 supplies the physical global formulation.

## Source transfer

Meng--Yang [2302.09741v1](https://arxiv.org/abs/2302.09741v1), sectionII C,
equations24--31, distinguish an Euler encoding from their incompressible
Schrödinger flow: the latter has an additional spin-gradient force. Their
TableI explicitly treats the fluid action parameter as arbitrary. Thus using
their quantum simulation algorithm is not a derivation of physical quantum
statistics for a fluid carrier.

Chern's [2017 thesis](https://cseweb.ucsd.edu/~alchern/projects/PhDThesis/)
likewise describes a Landau--Lifshitz modification and geometric Clebsch
variables. We use this as a representation lead; no theorem asserting
unmodified-Euler stability or universal spinor coverage is imported.

## 1. Local Berry form and the exact energy subtraction

Let z be a unit complex two-vector and kappa>0 a dimensional coordinate scale
with units of circulation (length squared/time). Define

    a=-i z^dagger dz,  u^flat=kappa a,
    s=z^dagger sigma z,  |s|=1.                              (1)

Here sigma are the Pauli matrices. This s is a coordinate texture, not an
electron spin measurement. In a regular chart write

    z=exp(i theta/kappa)
       (sqrt(1-f), sqrt(f) exp(i beta/kappa)), 0<f<1.

Direct differentiation gives

    u^flat=d theta+f d beta,
    curl u=grad f cross grad beta.                          (2)

This covers nonzero vorticity locally. A local Clebsch pair can be rescaled
and shifted so f lies inside (0,1) on a sufficiently small bounded chart;
it does not prove a single global chart for nonzero helicity.

For every spatial derivative, orthogonal decomposition of dz into its
component along z and its horizontal component gives

    |grad z|^2=|u|^2/kappa^2+|grad s|^2/4.

Consequently the PHYSICAL Euler energy is exactly

    H_E=(rho/2) integral |u|^2
       =(rho kappa^2/2) integral |grad z|^2
        -(rho kappa^2/8) integral |grad s|^2.             (3)

The subtraction is essential. The complete expression is nonnegative because
it equals the original kinetic energy; its two pieces cannot be varied as
independent fields. Replacing it by just the positive Dirichlet energy changes
the Hamiltonian. This is the concrete difference between a coordinate
representation and an additional spin-texture restoring law.

## 2. The same Euler action in the chart

The exact local first-order action is

    S=integral dt dx rho [i kappa z^dagger z_t-|u|^2/2]
     =-rho integral [theta_t+f beta_t+|grad theta+f grad beta|^2/2]. (4)

The imaginary scalar z^dagger z_t makes (4) real. Endpoint and spatial boundary
terms have the same vanishing/support conditions as0005. Variations of theta,
f and beta give, respectively,

    div u=0,  D_t beta=0,  D_t f=0.                     (5)

The pressure is the incompressibility multiplier in the equivalent constrained
Euler formulation. Local reconstruction of the momentum equation supplies

    D_t theta=|u|^2/2-p/rho,
    i kappa (partial_t+u dot grad)z=(p/rho-|u|^2/2)z.    (6)

One can derive (6) directly by varying normalized z with a pointwise multiplier:
the phase variation yields div u=0 and the horizontal variations transport the
texture. The unrestricted scalar multiplier determines p. Equation(6) and
incompressibility are therefore an exact local nonlinear Euler encoding, with
the velocity in(1), not a prescribed external transport field.

For an equation-level check, taking the material derivative of the one-form
in(2) yields

    (partial_t+L_u)u^flat=d(D_t theta)
                         =d(|u|^2/2-p/rho),

which is precisely Euler. Also D_t s=0. The pointwise spinor norm remains one;
it is a material-density constraint, not a localized one-particle probability.

General nontrivial-H1 domains retain the harmonic/circulation row from0005;
one cannot reconstruct full momentum by differentiating(2) and discarding all
closed one-forms. Multiple charts and their physical gluing must be supplied
for a selected carrier rather than assumed from this local computation.

## 3. Exact nonzero-vorticity example with no quantum inference

In a local Cartesian chart take theta=0, f=f(y) with 0<f<1, beta=kappa x.
Then u=(kappa f(y),0,0) is an exact stationary Euler shear with constant
pressure. The labels that realize its actual time evolution are

    f(t,x)=f(y),
    beta(t,x)=kappa x-kappa^2 f(y)t,
    theta(t,x)=kappa^2 f(y)^2 t/2.

They satisfy(5)--(6) with p=0; the y terms in d theta+f d beta cancel and the
velocity remains the same steady shear. This executes the nonlinear encoding
on a vortical Euler solution. It is local/periodic-background evidence, not
an isolated stable carrier. The same route must next glue actual carrier
charts, carry their orbit form and pressure, and identify realizable quantum
observations rather than substitute a different Hamiltonian.

## 4. What this positive construction earns

The local two-component Euler representation and exact energy/action map are
established by(1)--(6). They repair the scalar phase's irrotational restriction.
They do not derive a physical action quantum: kappa is a coordinate circulation
scale, and rho kappa integrated over volume has units of action. Neither an
SU(2)-valued chart nor the use of Pauli matrices selects half-integer physical
spin, a Born probability, exchange statistics or relativistic dispersion.
The next positive P4 construction is a controlled global carrier phase space
whose physical observables and interactions implement those structures under
the actual Hamiltonian(3). The parent objective remains unchanged.
