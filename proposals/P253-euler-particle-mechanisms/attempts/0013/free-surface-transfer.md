# Flat-interface distribution test of the free-surface Schrödinger map

The source inventoried before inspection is Zareei,
[arxiv2105.12253v1](https://arxiv.org/abs/2105.12253v1). Its first construction
uses an irrotational free surface, supplied surface tension and a Heaviside
wave amplitude. On page2 the real equation is tested against functions whose
normal derivative vanishes at the interface. That restricted test class
cannot establish the full Schrödinger distribution equation. The following
independent exact counterexample identifies the term it misses.

Set gravity to zero, take a flat stationary interface z=0, fluid below it,
velocity potential zero and pressure zero. With any constant positive surface
tension the curvature is zero, so the fluid is an exact resting free-surface
Euler solution. The source map assigns

    psi=H(-z), V=0, kappa=sqrt(2 sigma/rho)>0.

For the claimed Schrödinger equation the residual is

    R=i kappa psi_t+(kappa^2/2) Delta psi-V psi
     =-(kappa^2/2) delta'(z).                            (1)

It is nonzero as a distribution. Pair it with phi(z)=z exp(-z^2), multiplied
by a smooth compact tangential factor of integral one:

    <R,phi>=(kappa^2/2) phi'(0)=kappa^2/2.               (2)

The Gaussian can be replaced by z times any compact smooth bump equal to one
near zero. A test class with phi'(0)=0 deletes precisely this exposing term.
The stationary example needs no assumptions about a bouncing droplet or a
moving interface, and no products of singular distributions occur in(1).

**Route verdict:** the literal unregularized Heaviside mapping from all such
free-surface Euler solutions to the standard distributional Schrödinger
equation is refuted by the surviving double-layer distribution. This does not
refute free-surface Euler or controlled walking-droplet models. Surface tension
and an imposed material interface are themselves extra physical inputs for
the present whole-space constant-density substrate objective.

**Continuation:** smoothing the amplitude changes its interfacial gradient
energy and introduces a new layer equation; a singular Hamiltonian would need
an explicit operator domain and interface conditions. Either can be studied
as a different model, but neither is silently imported as a solution of the
frozen Euler-particle objective. The in-scope continuation executed in
spinor-euler-action.md retains actual Euler's energy in a regular vortical
coordinate chart, with no surface-tension substitution or Heaviside layer.
