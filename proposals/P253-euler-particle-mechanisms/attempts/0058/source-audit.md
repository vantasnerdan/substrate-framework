# Source audit

All full papers used here remain in `/tmp/primary-source-cache`; no downloaded
paper is copied into the campaign tree.  This attempt derives its new claims
directly and uses the papers only at the scopes in this table.

| Source | Receipt and exact scope used | Scope deliberately not imported |
| --- | --- | --- |
| Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation*, arXiv:2206.10165v2 | Cached `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`; SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`. Weighted meridional operator and Green kernel (2.1)--(2.4), PDF pp. 7--8; exact steady equation (3.1), p. 17; Lane--Emden core (3.5)--(3.7), p. 18; exact rescaling (3.11)--(3.12), pp. 19--20; limiting nondegeneracy Lemma 3.8, p. 24; Propositions 3.10--3.11, pp. 25--29; parameter system (3.36) and Proposition 3.13, pp. 29--30; free-boundary graph Lemma A.2, pp. 49--50; odd inverse Theorem C.2, pp. 57--58. This supplies the actual smooth compact-core axisymmetric translating carrier and its published regularity/axisymmetric stability scope. | No nonaxisymmetric Euler generator, critical-layer range theorem, KKS cluster, rotating branch, or nonlinear same-leaf normal form is attributed to Cao et al. |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices*, arXiv:1805.05064v3 | Cached `/tmp/primary-source-cache/P253-0040/1805.05064.pdf`; SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`. The paper's exact scope is a straight, infinite columnar vortex and its Fourier-mode/critical-layer spectral problem; see Theorem 1.3 and Sections 2--4 of the 52-page v3. It is used only to confirm that full Euler critical layers are coupled nonlocal mode equations, not bare scalar division. | Its columnar profile, spectral-stability theorem, boundary behavior at infinity, and mode exclusions are not transferred to the curved compact Cao ring. |
| Cao--Fan--Li--Qin, *Helical Kelvin waves for the 3D Euler equations*, arXiv:2411.02055v1 | Cached `/tmp/primary-source-cache/P253-0040/2411.02055.pdf`; SHA-256 `d6fbac9b2aa7f2394684807717224f31e3ae0afef8fed092dc7fa4def2223e4e`. Theorem 1.3, PDF pp. 5--6, constructs a curve of simply connected `m`-fold helical Kelvin **patches** from a disk for specified pitch/mode regimes; Theorem 1.8 treats doubly connected patches. This is a positive existence comparator showing that an exact contour variable can close a rotating Euler bifurcation in another carrier class. | The patch contour equation, helical symmetry, discontinuous vorticity, and its Crandall--Rabinowitz transversality do not give the smooth compact Cao-ring range or preserve the selected Cao coadjoint leaf. |
| P253/0048, 0052, and 0054 | Append-only author records. From 0048 only the formal fixed-order calculation and the proposed scale `j_crit=Theta((delta^2 log(1/delta))^-1)` are retained. The 0052 cluster and 0054 common-domain/Hodge claims are labeled conditional author inputs throughout 0058 pending the independent 0057 verdict. | No conditional IFT, graph-Riesz, or modewise bound is used to assert an exact kernel, complement inverse, nonlinear branch, or trace theorem. |

The arXiv source URLs and versions are part of the receipts:

- <https://arxiv.org/abs/2206.10165>
- <https://arxiv.org/abs/1805.05064>
- <https://arxiv.org/abs/2411.02055>

No production numerical datum, empirical comparator, or fitted coefficient is
used.  The verifier is exact finite-dimensional algebra and therefore does not
trigger the small-ratio numerical protocol.
