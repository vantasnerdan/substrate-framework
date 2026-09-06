# Exact current and quantum-map tests, with the missing physical construction exposed

This attempt tests particular field identifications. It neither assumes that all electromagnetic observables must be local Euler fields nor claims a no-go for emergence from a structured substrate.

## 1. Source role and access

Marmanis, [1998 Physics of Fluids](https://doi.org/10.1063/1.869762), equations(3)–(6) and the inviscid Lamb construction, uses vorticity and the Lamb vector as hydrodynamic counterparts of electromagnetic fields. His averaged propagation model introduces source/filter assumptions. The open primary PDF is at academicweb.nd.edu/~powers/ame.60635/marmanis.pdf. Its [erratum](https://doi.org/10.1063/1.869825) was identified but the primary text could not be retrieved in this pass. Accordingly no corrected viscous-current formula is imported: the inviscid identities below are independently derived from actual Euler.

Fusca, [arxiv1512.04611v2](https://arxiv.org/abs/1512.04611v2), Theorem3.5, establishes a Madelung momentum map to compressible quantum hydrodynamics. Its Hamiltonian contains a density-gradient term; it is not the constant-density vortical Euler Hamiltonian. The historical [Madelung translation](https://www.neo-classical-physics.info/uploads/3/0/6/5/3065888/madelung_-_hydrodynamical_interp..pdf) was accessible, while the independently derived action split below fixes our conventions.

## 2. Exact Lamb equations and conserved but neutral localized charge

For physical pressure p and density rho define

    B=omega=curl u, E=omega cross u,
    H=p/rho+|u|^2/2.

The vector identity for advection puts Euler in the exact form

    u_t=-E-grad H.
    div B=0, B_t=-curl E,
    q=div E=-Delta H.                                   (1)

For ANY chosen constant c0, define

    J_c=c0^2 curl B-E_t.                                 (2)

Then

    E_t=c0^2 curl B-J_c, q_t+div J_c=0.                   (3)

These are exact identities with a defined source, not a source-free propagation theorem. Changing c0 changes J_c by a divergence-free current and leaves u unchanged. Euler has therefore not selected a universal speed from this rewriting. If c0^2 is replaced by a variable coefficient, its gradient contributes to the divergence; the constant-coefficient continuity calculation cannot be reused silently.

Now take any smooth Euler field whose vorticity is compactly supported at the time considered, with the decaying Hodge velocity. Then E=omega cross u vanishes outside that support even though u has a dipole tail. Gauss gives the exact result

    integral_R3 q dx = integral_S_R E dot n dA = 0        (4)

for every enclosing sphere. For multiple disjoint compact vorticity components, every enclosing surface in the vorticity-free gap likewise has zero flux. Thus the literal Lamb-charge map assigns zero net charge to EACH such isolated component. This conclusion requires no weak-field or far-distance approximation.

More generally(4) follows if R^2 sup_{S_R}|omega cross u| tends to zero and q is integrable. A nonzero electric monopole under THIS map needs a nonvanishing asymptotic Lamb flux, a singular/nonlocalized source, a changed background limit or a different measured-field map. Assigning circulation its sign cannot repair(4). This is not a neutrality mechanism for the neutrino either: spin, weak coupling and flavor/mass dynamics remain absent.

The source energy also needs a map. Euler has rho integral |u|^2/2, while a Maxwell expression would involve an integral of |E|^2+c0^2|B|^2 with a dimensional conversion. Since B differentiates u and E is nonlinear in u, those functionals are not equal by(1). Their variational dynamics and the physical force on a carrier cannot be inferred from matching the homogeneous equations.

An exposing exact field is the steady periodic shear u=(sin y,0,0), p=constant. Then

    B=(0,0,-cos y), E=(0,-sin y cos y,0),
    q=-cos(2y), J_c=(c0^2 sin y,0,0).                    (5)

Its nonzero current is exactly what prevents(3) from predicting an autonomous Maxwell wave on this steady Euler field. Setting J_c=0 is an extra condition, not a consequence of the analogy or spatial averaging.

## 3. Scalar Madelung action: what is equal and what is additional

Let n>0 be a wavefunction NUMBER/probability density, not the fixed material rho. Set psi=sqrt(n) exp(i S/hbar), with positive mass parameter m. Direct differentiation gives

    (hbar^2/(2m)) |grad psi|^2
       =n |grad S|^2/(2m)
        +hbar^2 |grad n|^2/(8m n).                       (6)

The standard Schrödinger action becomes

    integral dt dx [-n S_t - n|grad S|^2/(2m)
                    -V n -hbar^2|grad n|^2/(8m n)].      (7)

Its variations yield n_t+div(n v)=0, v=grad S/m, and

    S_t+|grad S|^2/(2m)+V
       -(hbar^2/(2m)) Delta sqrt(n)/sqrt(n)=0.             (8)

The last term is the functional derivative of the additional gradient energy in(6), not a pressure law derived from the constant-density Euler kinetic action. The symplectic momentum-map property transfers a DECLARED Hamiltonian; it does not remove this difference or set hbar.

A regular single scalar phase also has curl v=0. Setting n=rho/m constant deletes the gradient term, but retains irrotationality and requires Delta S=0 for incompressibility. On all R3 with a regular single-valued phase and finite velocity energy, each component of grad S is an L2 harmonic function and hence zero. On a torus only allowed harmonic/circulation sectors survive. This scalar direct map cannot cover a smooth nonzero-vorticity carrier. Phase singularities, multiple charts, or spinor textures are materially different constructions; they must state their density/core and energy terms.

Conversely a nonlinear amplitude-dependent term can cancel the last term of(8) and encode classical irrotational dynamics in a Schrödinger-shaped equation. That is an exact representation change with nonlinear state dependence, not a derivation of linear quantum superposition or measurement probabilities.

## 4. Executed representation repair and its physical boundary

The vorticity/label action in0005 already repairs scalar irrotational coverage without introducing a new fluid law. Its state carries the actual Euler coadjoint structure and full Hodge Hamiltonian. A complex or spinor coordinate chart can represent that state, but the pulled-back HAMILTONIAN remains the actual nonlocal Euler functional; the canonical linear Schrödinger or Pauli Hamiltonian is not automatically substituted.

There is also an exact Hilbert representation of any measure-preserving classical flow: Koopman pullback acts unitarily on square-integrable functions of the classical state. For physical classical observables f,g represented by multiplication,

    M_f M_g=M_g M_f,                                     (9)
    expectation(M_f)=integral f |Psi|^2 dmu.

This is a positive exact amplitude representation. A local phase Psi -> exp(i theta(state))Psi changes none of these instantaneous expectations. Physical noncommuting measurements, an observable interference phase, a chosen probability measure/polarization and their interaction with carriers still require a construction. Introducing derivative operators into the Hilbert space does not identify them with realizable measurement operations by itself. No Bell-locality assumption about the elliptic Euler pressure is used here.

The failure-derived next candidate is a multicomponent/spinor Euler encoding with its full state-dependent pressure/spin-gradient correction retained, compared to the physical orbit action of0005. Newly discovered primary leads are Meng–Yang arxiv2302.09741 (hydrodynamic Schrödinger equation) and Chern's incompressible Schrödinger flow thesis. They are inventoried as new candidates before source-body reuse; possible extra spin forces must be exposed rather than treated as pure Euler or electron spin.

## 5. Route results and continuation

The Lamb homogeneous equations and defined conserved current are established; the literal compact-vorticity Lamb-charge identification with a nonzero electron charge is refuted by(4). Its extension to emergent electromagnetic fields remains blocked on an actual averaged/defect field, action, nonzero charge-flux and autonomous propagation/current construction. Retain the exact shear-memory source terms when testing a structured background.

The scalar Madelung kinetic/action decomposition is established; its direct regular constant-density identification with a localized vortical carrier is refuted by curl v=0 and the harmonic finite-energy consequence. The vorticity-label/Koopman repairs supply exact representations at their stated scope. Their physical quantum completion remains open, not disproved.

Neither neutral Lamb charge nor a two-component coordinate label supplies a neutrino. The same spin/statistics bridge, a chiral weak current coupled to the electron sector, distinct propagation/flavor states and a derived mass/mixing mechanism remain the positive P6 targets.
