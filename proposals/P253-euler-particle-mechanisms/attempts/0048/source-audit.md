# Source and applicability audit

Full source bodies remain in the primary-source cache and are not campaign
artifacts.

| Source | Receipt and exact scope used in 0048 |
| --- | --- |
| Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation* | arXiv:2206.10165v2; `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`; SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`. Exact meridional operator and Green kernel (2.1)--(2.4), PDF pp. 7--8; Lane--Emden core (3.5)--(3.7), p. 18; exact rescaled equation (3.11)--(3.12), pp. 19--20; first and odd second approximations, Proposition 3.10 and (3.29), pp. 25--27; Proposition 3.11, pp. 27--29; parameter system (3.36) and Proposition 3.13, pp. 29--30; limiting nondegeneracy Lemma 3.8, p. 24; free-boundary graph Lemma A.2, pp. 49--50; odd inverse Theorem C.2, pp. 57--58. These define the higher core cells and gauges. The paper states an `O(epsilon^2 |log epsilon|)` remainder, not the graph-domain Euler-generator jet derived below, and it supplies no nonaxisymmetric Riesz or branch theorem. |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices* | arXiv:1805.05064v3; `/tmp/primary-source-cache/P253-0040/1805.05064.pdf`; SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`. Fourier Euler equations (1.13)--(1.15), PDF pp. 3--4; Rayleigh function (1.19), p. 4; axisymmetric equation and energy identity (2.17)--(2.18), pp. 11--12; cylindrical Biot--Savart domain and estimate (6.1)--(6.9), pp. 36--38. Only those equations and estimates are used. Its ambient enstrophy realization and passive exterior-vorticity spectrum are not transferred to the compact-core DA domain. |
| Fischer--Schopohl, *Short wavelength spectrum and Hamiltonian stability of vortex rings* | arXiv:cond-mat/0008215v3; `/tmp/primary-source-cache/P253-0040/cond-mat-0008215.pdf`; SHA-256 `0b69b8cae5e758d860dcf8188260bfd79839555d173797d0a74ee487f30e5977`; mode equations (9)--(14), PDF pp. 2--3. Retained only as a regression comparator for the independently derived universal logarithmic compression. It is not used for a finite-core coefficient or Riesz projection. |
| P253/0044 derivation and result | Append-only predecessor artifacts. The exact DA support closure, core translation pencil, toroidal `P_1/P_2` Leray formulas, `m=0,2` audit, and universal logarithmic compression are rechecked inputs. The predecessor verdict is kept narrow: it refutes only a uniform fixed-rank rank-two reduction. |
| `0040/supervisor-core-response.md` | Author/unreviewed input retained untouched. Its sign and scalar response were independently rederived in 0044; 0048 does not promote the note itself as authority. |

No source states the finite-window theorem below. Its Sturm--Liouville
normalization, coupling estimate, resonance count, core-jet recurrence, and
Hamiltonian block estimate are derived in this attempt. No production
numerical values or empirical comparators are used.
