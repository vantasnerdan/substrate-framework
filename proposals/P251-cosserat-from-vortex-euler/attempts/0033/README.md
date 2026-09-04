# Attempt 0033: mutual-energy method repair

P251 N2/N3, original objective unchanged. Frozen route: repair the microscopic
functional used in 0026 before interpreting its parity as a verdict. Oracle:
exact geometry and full angle variation; no numerical remainder.

For decaying divergence-free vorticity with suitable core regularization, put
A=(-Delta)^-1 omega, u=curl A. Integration by parts gives
E=rho/2 integral |u|^2=rho/2 integral A·omega. Two closed filament components
therefore have the mutual energy

    E12 = rho Gamma1 Gamma2/(4pi) double-integral (t1·t2)/|X1-X2| ds1 ds2.

In contrast, integral u1·omega2=Gamma2 integral u1·dX2 is a helicity
contribution. Attempt 0026 used the latter velocity line integral as energy.
Its stated squared distance also differs from its declared two-line geometry.
The parity conclusion therefore did not eliminate the energy interaction route.
For infinite straight lines the mutual kinetic energy additionally needs an
outer-domain prescription; a finite tangent result is not that renormalization.

The code constructs the declared geometry and differentiates the actual
energy kernel. Under reflected tilt theta->-theta it is EVEN. At parallel
geometric tangents the angle Hessian per integration element is

    -(d^2+s1^2+s2^2-s1*s2)/(d^2+(s1-s2)^2)^(3/2).

The numerator is a strictly positive sum of squares. With opposite signed
circulations the finite-segment mutual energy has a strictly positive angular
Hessian. This is a genuine local angle-dependent interaction, independent of
the patch route in 0030. It cannot be annulled by the parity argument in 0026.

`route_verdict: established` for the kernel identity and finite-segment
contribution. `evidence_scope: EXACT_LOCAL_INTERACTION_CONTRIBUTION`.
The complete Euler-frame-locking route stays active: close the vorticity into
loops or a specified network, include the connecting segments in the energy,
find a relative equilibrium, and reduce its Euler symplectic dynamics before
identifying inertia or a continuum modulus. A rigid finite segment by itself
is not a divergence-free Euler vortex. No claim of pair stability, global
positive Hessian, continuum closure or exhaustion is made.

This method repair reopens the inter-tube energy family rather than accepting
0026's route verdict. Closed contour-corrected patch pairs are pursued in
parallel in 0032. The nearest imported construction is the explicit
Biot-Savart line-energy integral in src/substrate_framework/homogenization.py;
that self-energy supplies the normalization convention, not the mutual action.

Independent primary corroboration of the tangent-dot-tangent energy kernel:
Bustamante and Nazarenko, *Derivation of the Biot–Savart equation from the
Nonlinear Schrödinger equation* (2015),
https://arxiv.org/abs/1507.07806. Its model-specific core cutoff is not imported;
the displayed integration-by-parts derivation uses the declared Euler kernel.

Reproduce with PYTHONPATH=src .venv/bin/python
proposals/P251-cosserat-from-vortex-euler/attempts/0033/verify_mutual_kernel.py.
