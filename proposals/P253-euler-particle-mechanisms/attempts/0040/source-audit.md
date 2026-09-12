# Activated primary-source and accepted-API audit

All downloaded bodies remain in `/tmp/primary-source-cache/P253-0040`, outside
the campaign artifacts.  The README receipts for the two Cao families remain
in force.  Activation opened the following additional bodies.

| Source | Version, hash, exact location | Imported result and boundary |
| --- | --- | --- |
| Fischer--Schopohl, *Short wavelength spectrum and Hamiltonian stability of vortex rings* | [arXiv:cond-mat/0008215v3](https://arxiv.org/abs/cond-mat/0008215v3), PDF SHA-256 `0b69b8cae5e758d860dcf8188260bfd79839555d173797d0a74ee487f30e5977`; line action and Biot--Savart cutoff (1)--(8), PDF pp. 1--2; mode equations and coefficients (9)--(14), pp. 2--3; quadratic Hamiltonian (17), p. 4 | Exact dynamics **within a geometrically cut-off filament model**, including its canonical form and fixed-mode thin-ring asymptotics.  The paper explicitly says that the core model changes order-one constants.  It does not construct a smooth finite-core Euler eigenfunction or branch. |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices* | [arXiv:1805.05064v3](https://arxiv.org/abs/1805.05064v3), PDF SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`; exact full column linearization (1.10)--(1.15), PDF pp. 3--4; Theorem 1.3 and Remarks 1.4--1.5, pp. 6--7; Kelvin modes and the `m=1` translation branch, Section 3.2, pp. 14--16; Rankine dispersion (6.10)--(6.17), pp. 38--39 | The exact column operator has an essential interval on the imaginary axis and simple Kelvin eigenvalues outside it.  For cross-sectional mode `m=1`, the zero translation mode at axial wavenumber zero bifurcates outside the essential interval for small nonzero axial wavenumber.  Their global theorem assumes profile conditions H1--H2; applicability to the compact Cao limit profile must be proved, not inferred from the Rankine or Lamb--Oseen examples. |
| Cao--Fan--Li--Qin, *Helical Kelvin waves for the 3D Euler equations* | [arXiv:2411.02055v1](https://arxiv.org/abs/2411.02055v1), PDF SHA-256 `d6fbac9b2aa7f2394684807717224f31e3ae0afef8fed092dc7fa4def2223e4e`; helical reduction (1.3)--(1.7), PDF pp. 2--3; Theorem 1.3, pp. 5--6; contour functional (3.1)--(3.2), pp. 12--13; exact multiplier Lemma 3.3 and Proposition 3.5, pp. 19--25 | A rigorous, exact 3D Euler contour-dynamics/CR construction with a one-dimensional Fourier kernel, closed codimension-one range, and explicit transversality.  Its screw-periodic helical symmetry forces nonzero solutions to repeat indefinitely along the axis, so it is not compactly supported or finite-energy on `R3` and is not a Cao ring. |
| Pocklington, *The complete system of the periods of a hollow vortex ring* | DOI `10.1098/rsta.1895.0017`, volume 186 (1895), 603--619 | Royal Society metadata is available, but the publisher PDF returned HTTP 403 in this environment.  No unseen equation was imported.  The accessible Fischer--Schopohl body advances Route B with an explicit model rather than treating the paywall as a stopping condition. |

## Exact Cao local profile used in Route B

Cao--Lai--Qin--Zhan--Zou define the radial cross-section ground state by

    -Delta U=U^p in B_1,  U>0,  U=0 on partial B_1,           (1)

and extend it harmonically outside; see their (3.5)--(3.7), PDF p. 18.
After the source scaling, the local column vorticity profile is proportional
to

    W_p(s)=U(s)^p for 0<=s<1,  W_p(s)=0 for s>=1.             (2)

It is radially decreasing inside the core.  Since `U` vanishes linearly at
the edge, choosing integer `p>=2` gives at least the `C1` edge regularity
needed to write the Gallay--Smets column operator.  Their H2 Richardson
monotonicity is proved directly from the Lane--Emden ODE in `derivation.md`;
the curved-ring threshold resolvent remains a separate obligation for (2).

## Accepted abstract APIs inspected, not promoted

- `src/substrate_framework/euler_bent_tube.py` exactly pushes a straight
  column vorticity through a volume-preserving bending map and computes its
  relative impulse density.  Its own contract says that it does not make the
  curve an Euler trajectory or persistent finite-core soliton.  It is a
  useful tangent parametrization, not the 0040 branch.
- `src/substrate_framework/euler_column_wave.py` returns conditional column
  linearization/exterior algebra.  Its principal carrier is a pure-swirl
  straight column, not the translating no-swirl Cao torus, and it does not
  supply the threshold spectral projection below.
- `src/substrate_framework/euler_orbit.py` consumes already-computed
  Hessian/KKS pairings.  It cannot manufacture the eigenmode, the coadjoint
  chart, or the nonzero pairing required here.

No accepted algebraic API is treated as a physical realization of the
registered branch.
