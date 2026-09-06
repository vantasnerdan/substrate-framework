# 0039 primary-source and dependency inventory

## Published microlocal source

Shvydkoy--Vishik, *On Spectrum of the Linearized 3D Euler Equation*,
Dynamics of PDE 1 (2004), 49--63. Primary PDF:
`https://intlpress.com/api/bgcloud-front/resource/pdf/volume/1805802424294064130-1805802424294064130-34419541abcc4a0b17c7e7b4a66a7c98.pdf`.
Cached outside the campaign tree at
`/tmp/primary-source-cache/P253-0039/Shvydkoy-Vishik-2004.pdf`, SHA-256
`9307efc6096565f74c5a099af25c834af09ba93b3e56ff34dd076ba289551592`.

Exact locations: velocity-space linearized Euler operator and energy space,
PDF pp. 2--3 (journal pp. 50--51); bicharacteristic-amplitude system (1.5),
PDF p. 3 (journal p. 51); essential spectral-radius attribution and Theorem
1.1, PDF pp. 3--4 (journal pp. 51--52); weak-null WKB sequence description,
PDF p. 4; Theorem 1.2, PDF p. 4.

Scope audit: the paper proves its stated theorems on `T^3`. It says free space
should present no major difficulty but does not prove that extension in the
displayed theorem. Therefore 0039 does not import its torus essential-spectrum
conclusion as a theorem about the Gavrilov flow on `R^3`. It uses the published
BAS and weak-null construction as corroboration, while deriving the fixed-time
whole-space packet estimate directly from compact coefficients, the Fourier
Leray projector, and nonstationary phase.

The paper cites M. Vishik, *Spectrum of small oscillations of an ideal fluid
and Lyapunov exponents*, J. Math. Pures Appl. 75 (1996), 531--557, for the
essential-radius theorem. No full primary body was available in the activated
inventory, so no theorem from that paper is imported beyond the bibliographic
attribution in Shvydkoy--Vishik.

## Local proposal dependencies

- P253/0032 supplies the corrected physical BAS return on the selected
  Gavrilov shell: `det M_*=1` and
  `Delta_*=88*pi^2 I_*+O(I_*^(3/2))>0`.
- P253/0019 supplies the whole-space skew-transport plus bounded-shear group,
  finite-energy Biot--Savart normalization, and one fixed-time WKB calculation.
- P253/0025 supplies the profile-return, global Leray/Hodge/collar ledger, and
  the permitted quantifier `for every finite circuit count, choose frequency
  afterward`.

All three are active-proposal evidence, not accepted canon. This attempt adds
the missing all-circuit diagonal packet and essential-norm argument; it does
not promote any parent particle claim.
