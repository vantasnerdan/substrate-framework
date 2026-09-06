# Primary-source access inventory

Inventory frozen on 2026-09-06 before downloading any new full source body for this attempt.
HTTP status was tested through the final redirect. Full papers, when opened, belong only in
`/tmp/primary-source-cache/P253-0010`.

| Candidate/source | Primary endpoint | Access result | Permitted use in 0010 |
| --- | --- | --- | --- |
| Cao--Zhan, *On the steady axisymmetric vortex rings for 3-D incompressible Euler flows*, arXiv:2009.13210v2 | <https://arxiv.org/pdf/2009.13210v2> | HTTP 200, `application/pdf`; arXiv API identifies v2 dated 2020-11-28 | Full theorem/equation audit for the published regular axisymmetric family with and without swirl; retain the theorem's actual Sobolev/vorticity regularity rather than upgrade it to `C-infinity`. |
| Turkington, *Vortex Rings with Swirl: Axisymmetric Solutions of the Euler Equations with Nonzero Helicity*, DOI 10.1137/0520005 | <https://epubs.siam.org/doi/10.1137/0520005> | Landing page and abstract accessible; PDF endpoint returned HTTP 403 / access page | Metadata/abstract scope only unless a lawful full-text copy is located. The DOI paper is authored by Bruce Turkington alone; “Friedman--Turkington” describes the surrounding variational lineage, not this paper's byline. No theorem detail will be invented from the abstract. |
| Constantin--La--Vicol, *Remarks on a paper by Gavrilov...*, arXiv:1903.11699 | <https://arxiv.org/pdf/1903.11699> | HTTP 200, `application/pdf`; author copy also available at <https://cims.nyu.edu/~vicol/CLV1.pdf> | Accessible materially different candidate: smooth compactly supported axisymmetric steady Euler flow with real meridional and swirl components via Grad--Shafranov localization. |
| Gavrilov, *A steady Euler flow with compact support*, arXiv:1810.08020 | <https://arxiv.org/pdf/1810.08020> | HTTP 200, `application/pdf` | Original compactly supported smooth steady construction and exact local field formula; audit helicity, angular momentum, and variational/stability content. |
| Enciso--Peralta-Salas, *Existence of knotted vortex tubes in steady Euler flows*, arXiv:1210.6271 | <https://arxiv.org/pdf/1210.6271> | HTTP 200, `application/pdf` | Accessible helical/knotted Beltrami comparator. Test decay, energy/action, angular momentum, and stability before particle use. |
| Lifschitz--Suters--Beale, *A Numerical and Analytical Study of Vortex Rings With Swirl*, DOI 10.1051/proc:1996017 | <https://doi.org/10.1051/proc:1996017> | The DOI identifies an open ESAIM proceedings article, but its publisher endpoint returned an automated-access HTTP 403 and the tested CiteSeerX mirror returned 404 in this environment | Metadata-level stability diagnostic only in 0010. No result from its body is used. |

The SIAM paywall therefore does not stop the objective. Cao--Zhan supplies the named smooth
ring existence route, while Gavrilov/Constantin--La--Vicol supply an accessible exact compact
steady route with meridional circulation and swirl. Enciso--Peralta-Salas supplies the distinct
helical topology route, subject to finite-action and stability tests.

Post-inventory access resolution: no full Lifschitz--Suters--Beale body was obtained, so it is
excluded from the evidence chain below. The accessible Cao--Zhan, Gavrilov, Constantin--La--Vicol,
and Enciso--Peralta-Salas primary bodies suffice for the registered comparison.
