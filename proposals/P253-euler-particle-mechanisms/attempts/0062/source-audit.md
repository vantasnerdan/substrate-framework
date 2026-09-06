# Source and dependency audit

All full papers remain in `/tmp/primary-source-cache`. No paper or extracted
full text is copied into the campaign artifacts.

| Source | Receipt and exact use | Boundary retained |
| --- | --- | --- |
| Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation*, arXiv:2206.10165v2 | `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`, SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`. Equations (1.3)--(1.6), PDF pp. 1--2, give the no-swirl potential vorticity `zeta=omega_theta/r`, meridional velocity, and translating steady equation. Equations (3.1)--(3.7), pp. 17--18, give the exact fixed carrier equation and radial limiting cell. The free-boundary and local uniqueness results cited in 0057 are used only at their reviewed scopes. | The paper supplies no nonaxisymmetric Euler graph domain, full-Hodge harmonic resolvent, limiting absorption, adjoint critical trace, KKS mode, or rotating branch. Its auxiliary parameter equation (3.36) is not an exact augmented Banach map. |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices*, arXiv:1805.05064v3 | `/tmp/primary-source-cache/P253-0040/1805.05064.pdf`, SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`. Equations (2.1)--(2.5), PDF pp. 9--10, decompose the straight-column Fourier generator into a multiplication block plus a compact Biot--Savart perturbation in the source enstrophy space. Equations (2.8)--(2.16), pp. 10--11, show explicitly that eliminating velocity produces a coupled singular critical-layer ODE rather than bare scalar division. | Compactness, spectral stability, monotone profile, infinite-column boundary conditions, and the source space `X_(m,k)` are not transferred to the curved compact Cao ring. The paper does not provide the 0062 limiting-absorption theorem. |
| P253/0005 | `derivation.md`, equations (1)--(5), fixes the right-reduced Euler coadjoint convention and physical density: `Omega_KKS(X_xi,X_chi)=-<m,[xi,chi]>`. | It supplies the general orbit form, not a Cao spectral mode or its normalization. |
| P253/0057 joined independent review | `review.md` and `verdicts.yaml`. Unit A preserves the corrected scaling/local cells/unequal-volume representation/KKS normalization but blocks the all-sector graph inverse on `GR`. Unit B preserves pointwise shape and metric-Leray identities plus the limited `p>=6` no-sheet trace, while blocking `HJ2`, the common graph domain, two-index bounds, and exact Riesz transfer. | The conditional implication `HJ2+GR => transfer` is not used as either antecedent. No 0052/0054 IFT, graph inverse, mode, or Riesz projector is assumed. |
| P253/0058 | `derivation.md`, `result.yaml`, and its exact finite algebra receipt. The carrier-map equation, positive-core centralizer, reverser/KKS sign, DA/frozen Hodge symbols, and universal nondivisibility counterexample are retained. | Its generic pendulum equations are conditional candidates only. They did not derive the Cao full-Hodge `alpha`, `V_*`, or symplectic measure. |

Primary URLs and versions:

- <https://arxiv.org/abs/2206.10165>
- <https://arxiv.org/abs/1805.05064>

No production numerics, comparator fitting, or filament substitution is used.
The only executable check is exact SymPy algebra with the repository
interpreter, so the small-ratio numerical protocol is not activated.
