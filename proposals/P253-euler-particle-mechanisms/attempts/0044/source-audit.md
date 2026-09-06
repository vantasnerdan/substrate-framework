# Source and applicability audit

Full papers remain outside the campaign tree.

| Source | Receipt and exact scope used |
| --- | --- |
| Cao--Lai--Qin--Zhan--Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation* | arXiv:2206.10165v2; cached `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`; SHA-256 `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`. Ground-state definition (3.5)--(3.7), PDF p. 18; first error estimate Proposition 3.10, (3.25), pp. 25--27; odd second approximation (3.29) and Proposition 3.11, pp. 27--29; parameter estimate Proposition 3.13, pp. 29--30; odd linear cell theorem C.2, Appendix C. These give the compact Lane--Emden core and an order-one odd correction with an `O(epsilon^2 |log epsilon|)` remainder. They do not give the graph-norm second spectral jet or a 3D ring Riesz theorem. |
| Gallay--Smets, *Spectral stability of inviscid columnar vortices* | arXiv:1805.05064v3; cached `/tmp/primary-source-cache/P253-0040/1805.05064.pdf`; SHA-256 `081e93124ad0a8e04a51842848eb220180f904f6c3eec27e9b985a9a636a5797`. Full Fourier-fiber equations (1.13)--(1.15), PDF pp. 3--4; operator split and ambient essential spectrum (2.1)--(2.6), pp. 9--10; scalar radial equation (2.16), p. 11; axisymmetric equation (2.17)--(2.18), pp. 11--12; planar `m=1` translation statement, Section 3.2, p. 16. The ambient passive-vorticity endpoint is not transferred to the compact coadjoint closure. Equation (2.17) is used to derive the actual internal `m=0` compact generalized eigenproblem. |
| Fischer--Schopohl, *Short wavelength spectrum and Hamiltonian stability of vortex rings* | arXiv:cond-mat/0008215v3; cached `/tmp/primary-source-cache/P253-0040/cond-mat-0008215.pdf`; SHA-256 `0b69b8cae5e758d860dcf8188260bfd79839555d173797d0a74ee487f30e5977`; mode equations (9)--(14), PDF pp. 2--3. It is used only as a regression target for the independently derived leading length/impulse matrix, never for the Cao finite-core constant or graph-domain eigenfunction. |
| `0040/supervisor-core-response.md` | Author/unreviewed input, SHA-256 at consumption `f472b1be7b69f5140fdc1a6c57cbfbf0eeb666eb408a2cfa7805eb773aaedfbf`. Its formulas were independently rederived in `derivation.md` and checked by the attempt verifier. It is not treated as accepted authority and was not edited. |

The existing `euler_orbit`, `euler_column_wave`, and `euler_bent_tube` APIs
were already audited in 0040. They provide conventions or conditional algebra,
not the curved pressure blocks, compact graph domain, or Riesz projection
constructed/tested here.
