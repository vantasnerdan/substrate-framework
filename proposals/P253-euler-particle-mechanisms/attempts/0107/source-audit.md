# P253/0107 source and authority audit

## Reviewed framework inputs

- P253/0095 as corrected and independently reviewed in P253/0104 supplies the
  fixed-member coupled BAS/Egorov result, the exact coadjoint/tag/Gauss path,
  and the same-`H^s` packetwise differentiation bridge.  It does **not** supply
  the conserved fixed-row curve; Sections 2--9 of `derivation.md` construct the
  corrected `J_ren` row and its witnesses.
- P253/0068 as reviewed in 0071 supplies the material current, Euler--Lorentz
  and Maxwell equations, and local constraint/energy propagation.
- P253/0077 as reviewed in 0081 supplies the axisymmetric traveling fields
  `B_g=rH_g e_theta`, the strict-band material tag, and the subluminal
  vector-potential/Maxwell inverse.
- P253/0080 as corrected and reviewed in 0084 supplies one fixed nonzero
  charged Cao member with `E_g,B_g=O(g)`.  It supplies no stability or
  fixed-leaf witness theorem.

## Primary sources

### Cao--Lai--Qin--Zhan--Zou

*Uniqueness and stability of steady vortex rings for 3D incompressible Euler
equation*, arXiv:2206.10165v2.

- URL: https://arxiv.org/abs/2206.10165
- cached PDF: `/tmp/primary-source-cache/P253-0040/2206.10165.pdf`
- PDF SHA-256:
  `6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca`
- exact locations: equations (1.3)--(1.6), PDF pp. 1--2, for the no-swirl
  translating vorticity law; Proposition 1.4, PDF pp. 5--6, for the compact
  polynomial ring; Proposition 3.2 and equation (3.8), PDF pp. 17--21, for
  scaled `C1` convergence; Lemma A.2 and (A.5)--(A.8), PDF pp. 49--50, for the
  nondegenerate free-boundary level.

Consumed scope: the compact smooth-enough fixed carrier and its boundary
vanishing.  The source does not prove the coupled momentum identity, compact
H/C/G witnesses, nonaxisymmetric P2, or the charged-column graph theorem.

### Hattori--Fukumoto

Y. Hattori and Y. Fukumoto, *Short-wavelength stability analysis of thin
vortex rings*, Physics of Fluids **15** (2003), 3151--3163,
DOI `10.1063/1.1606446`.

- URL: https://doi.org/10.1063/1.1606446
- cached PDF:
  `/tmp/primary-source-cache/P253-supervisor/hattori-fukumoto-2003.pdf`
- PDF SHA-256:
  `c6a35c44d55d26b8bb3e8fc9e55b29b638d0983b646b17eb522a4192b6444fb9`
- exact locations: BAS equations (2)--(4) and Hill reduction (7)--(9), PDF
  p. 3152; first harmonic (22), resonance (23), and exponent (24), PDF
  pp. 3156--3157; arbitrary first-order flow/covector cancellation in Appendix
  C, PDF pp. 3161--3162.

Consumed scope: the reviewed `b(0)=15|U_theta|/256` local first-curvature
coefficient used in (49a).  It is not a whole-core, free-boundary, charged, or
nonlinear stability theorem.

### Choi

K. Choi, *Stability of Hill's spherical vortex*, arXiv:2011.06808v2.

- URL: https://arxiv.org/abs/2011.06808
- cached PDF: `/tmp/primary-source-cache/P253-0095/2011.06808.pdf`
- PDF SHA-256:
  `0b2fe63374837db46ac054e81440be94629352d49c6e25414885983c51726783`
- exact locations: Theorems 1.1--1.2, PDF pp. 4--5; scope discussion and
  Remarks 1.3--1.6, PDF pp. 3--6.

Consumed scope: all-time stability only in the nonnegative axisymmetric
no-swirl relative-vorticity class and its stated metric.  It does not contain
the full nonaxisymmetric bulk/interface/exterior Hessian (57)--(58).

### Swirling-ring candidate inventory

1. A. Lifschitz and E. Hameiri, *Localized instabilities of vortex rings with
   swirl*, Communications on Pure and Applied Mathematics **46** (1993),
   1379--1408, DOI `10.1002/cpa.3160461005`, publisher page
   https://onlinelibrary.wiley.com/doi/abs/10.1002/cpa.3160461005.
   The abstract states the geometric-optics criterion and exponential or
   algebraic localized growth for axisymmetric rings with swirl.  It is a BAS
   filter, not a coercive nonlinear theorem.
2. A. Lifschitz, W. H. Suters, and J. T. Beale, *The onset of instability in
   exact vortex rings with swirl*, Journal of Computational Physics **129**
   (1996), 8--29, DOI `10.1006/jcph.1996.0230`, metadata/abstract page
   https://scholars.duke.edu/publication/759533.  It compares localized
   short-wave predictions with time-dependent vortex-method simulations on
   computed exact rings.  It is evidence that swirl is not automatically a
   restoring mechanism, not a coverage theorem for all swirl carriers.
3. B. Turkington, *Vortex Rings with Swirl: Axisymmetric Solutions of the
   Euler Equations with Nonzero Helicity*, SIAM Journal on Mathematical
   Analysis **20** (1989), 57--73, DOI `10.1137/0520005`, publisher page
   https://epubs.siam.org/doi/10.1137/0520005.  The abstract establishes
   variational existence of bounded-domain and steadily translating
   axisymmetric swirling rings.  It does not establish full-three-dimensional
   coercivity or P2.

No full paywalled paper is copied into the campaign tree.  These entries record
URL, version/DOI, and exact consumed scope; any later theorem-level use beyond
the abstracts requires a cached primary PDF with its hash and exact location.

## New derivations and scope

- Equations (3)--(10) derive the coupled impulse/Maxwell exchange and the
  conserved `J_ren` row.
- Equations (15)--(25) prove the quantitative cutoff-plus-mollification witness
  lemma; (26)--(38) apply it to the actual Cao `B_g`, `omega_g`, and tag center.
- Equations (39)--(46) close the high-character tangent rows and formulate the
  packetwise, not full-operator, nonlinear contradiction.
- Equations (47)--(52) establish the frozen-column center tongue, edge
  conditioning, and the precise active accessibility/interface transfer.
- Equations (53)--(58) derive the weak full Hill orbit Hessian while leaving its
  nonaxisymmetric VSH radial sign active.

No source or reviewed input is cited as a proof of P2, nonlinear instability,
quantization, or an electron/neutrino mechanism.
