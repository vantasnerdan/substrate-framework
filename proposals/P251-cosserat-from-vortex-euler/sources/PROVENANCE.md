# P251 source provenance

All primary sources verified in hand 2026-09-03: PDF saved locally, text layer
extracted (sidecar `.txt`, `=== PAGE BREAK ===` per page), md5 recorded in
`MD5SUMS`. Committed here: the four load-bearing PDFs, every extracted text,
and the hash manifest. Exceptions are listed individually.

| arXiv id | Paper | Role | Access |
| --- | --- | --- | --- |
| 1210.6271 | Enciso & Peralta-Salas, "Existence of knotted vortex tubes in steady Euler flows" (2012) | Load-bearing: existence backbone for stationary vortex tubes of arbitrary knotted topology | PDF + txt committed |
| 1003.3122 | Enciso & Peralta-Salas, "Knots and links in steady solutions of the Euler equation" (2010, Ann. Math. 175 (2012)) | Load-bearing: arbitrary knotted/link vortex lines in steady Euler flows | PDF + txt committed |
| 1505.01605 | Enciso, Peralta-Salas & Torres de Lizaur, "Knotted structures in high-energy Beltrami fields on the torus and the sphere" (2015) | Load-bearing: explicit Beltrami eigenfields carrying arbitrary knotted structures (realized-example pool for C-CST-007) | PDF + txt committed |
| 2103.14458 | Rudge, "A micropolar continuum model of diffusion creep" (2021) | Load-bearing: in-hand statement of the target linear micropolar equations (microrotation, couple stress, spin balance); methodology precedent: micropolar continuum coarse-grained from discrete rotational degrees of freedom | PDF + txt committed |
| 2401.09506 | Tekinalp, Bhosale, Cui, Chan & Gazzola, "Soft, slender and active structures in fluids: embedding Cosserat rods in vortex methods" (2024) | Context only (rod-in-fluid coupling methods); no claim imports from it | Open at https://arxiv.org/abs/2401.09506 ; md5 08ba02c27846ec110d9c1a21bf5381cc; PDF (37 MB) NOT committed (repo-bulk disproportionate for a context source); local working copy at ~/sources/cosserat-euler/2401.09506.pdf; txt committed |

Cited, not imported (not load-bearing): A. C. Eringen, "Linear theory of
micropolar elasticity", J. Math. Mech. 15(6) (1966) 909-923 (paywalled; the
target equations are stated by in-hand 2103.14458 and re-derived by this
campaign as its own claims); Cosserat & Cosserat, "Théorie des corps
déformables" (1909, historical).

Classical statements imported under the P242 precedent: MIT 2.25 §9
(Kelvin's circulation theorem, Biot-Savart field, vortex lines move with the
fluid), verified open https://web.mit.edu/2.25/www/225_sect_9.html
(verified 2026-08-22 in P242; re-verified 2026-09-03).
