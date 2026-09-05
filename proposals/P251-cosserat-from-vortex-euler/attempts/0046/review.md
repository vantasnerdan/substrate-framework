# Independent individual review, with a derived finite-domain appendix

Reviewer: Codex `/root/construction_review`, 2026-09-05.
Author of 0046: Codex `/root`; this reviewer did not implement its field,
common-angle mechanism, four-coordinate action, or verifier. The reviewer
authored dependency 0045; that dependency is not independently re-reviewed
here. The estimates explicitly marked below are **review-derived proof
completion**, not evidence of a second independent check of themselves.

Reviewed `intrinsic-rotor.md` and `verify.py` (the recorded 15/15 exact
checks). No test rerun or shared-file edit was needed. The scope is this
finite-cylinder common/relative Euler-orbit construction, not a new review
of prior attempts or an arbitrary-wavelength invariant-subspace requirement.

## Individual verdicts

1. **Stationary finite-core field: established as stated.** With the stated
   cylindrical velocity signs, the Bessel equation gives curl u=lambda u;
   J2(lambda R)=0 makes psi constant on the sidewall, hence u.n=0. The
   noncircular elliptic quadratic jet exists for A>0 and 0<|B|<2A. The solid
   cylinder and nearby core level tubes are material invariant domains.
2. **Euler orbit and boundary Hessian: established as stated.** If xi and
   omega are tangent to the wall, xi cross omega is normal there. Its
   Leray projection therefore has tangential trace -d_boundary pi, an exact
   one-form. The background trace is closed because omega.n=0. The boundary
   curl pairing vanishes by integrating a closed form wedged with an exact
   form over the boundary torus. Perturbation/perturbation pairing vanishes
   likewise. Fixed boundary periods select the required harmonic sector;
   no arbitrary deletion of axial mean momentum is warranted.
3. **Common symmetry and momentum partner: established as stated.** The
   rotation tangent vK is itself Beltrami. In particular, after the boundary
   cancellation,
   `H(K,eta)=rho integral(vK-curl(vK)/lambda).v_eta=0` exactly. This supplies
   a direct equation check in addition to rotational energy invariance.
   With eta0=curl(chi curl vK), compact support gives
   `Omega(K,eta0)=-rho integral chi |curl vK|² != 0`. Its sign follows from
   omega.(K cross eta)=-(K cross omega).eta and Leray self-adjointness.
4. **Positive intrinsic common inertia and separate internal angle:
   established with the finite-bound proof below.** A high nonzero axial
   harmonic raises the common momentum partner's energy without changing
   its angular-momentum derivative. A distinct harmonic provides the full
   positive internal H and KKS pair. Exact axial orthogonality removes every
   cross block, including nonlocal projection terms. Both inertias are then
   obtained from one Euler orbit action, not from an added rigid mass.
5. **Covariant section-angle map: established with one explicit support
   choice.** Choose chi for eta0 in an annulus away from the physical core
   jet, and put both high-harmonic cages away from that jet. This is possible
   because the analytic nonzero vK cannot vanish throughout every such
   annulus. The conjugate generators then do not alter the core's local
   material orientation. The actual section angles are B+q and B-q, and the
   displayed field map and both positive kinetic coefficients follow.

The cylinder is not identified with an arbitrary EPS knotted-tube ambient
domain by this result. That distinction is already in the authored claim;
no additional all-k or free-Euler-invariant-subspace condition is imposed.

## Review-derived appendix A: Neumann Leray estimate on the cylinder

Let D=S cross S1_L, with S the transverse disk. Use a fixed length unit and
write k=2 pi n/L, n a nonzero integer. All amplitudes below are smooth and
compactly supported inside S, independent of z. Complex notation is a
convenience; real parts give the physical fields. Norms include the axial
period, so Plancherel normalization never changes coefficients silently.

For `f=F(x,y) exp(ikz)`, let P_D f=f-grad pi. The pressure amplitude obeys

    (Delta_perp-k²) pi = div_perp F_h + i k F_z,
    partial_n pi = 0.

The homogeneous Neumann data follow from compact support. Put
`pi0=F_z/(i k)` and `r=pi-pi0`. Since pi0 also vanishes near the wall,

    (Delta_perp-k²) r = g,
    g=div_perp F_h-Delta_perp F_z/(i k),
    partial_n r=0.

Multiply by conjugate r and integrate. There is no boundary term, and

    ||r||_2 <= ||g||_2/K²,
    ||grad_k r||_2 <= ||g||_2/K,  K=|k|,
    grad_k=(partial_x,partial_y,i k).

Consequently, with Pi=I-ez ez^T,

    ||P_D f-Pi f||_2
      <= [||grad_perp F_z||_2+||div_perp F_h||_2
          +||Delta_perp F_z||_2/K]/K.

This is the required finite-domain O(1/K) bound; it uses the actual Neumann
operator, rather than importing the R³ Fourier symbol. Since curl P_D=curl,

    curl(P_D f-Pi f)=curl(F_z ez exp(ikz)),
    ||curl(P_D f-Pi f)||_2 <= ||grad_perp F_z||_2.

The fast derivative cancels because F_z ez is axial. These bounds prove
exactly the bounded curl error needed in the physical helicity Hessian.

For an exact compact curl cage, write its force as
`f=(F0+F1/k) exp(ikz)`, with amplitudes independent of k. For K>=1 one may
use the explicit constants

    D_i=||grad_perp F0_z||_2+||div_perp F0_h||_2
        +||Delta_perp F0_z||_2+||F1||_2,
    E_i=||grad_perp F0_z||_2+||F1||_2+||curl_perp F1||_2.

Here curl_perp means the usual curl with the z derivative set to zero.
The F1 term uses ||P_D||<=1 and
`curl(F1 exp(ikz))/k=i ez cross F1 exp(ikz)
                           +curl_perp F1 exp(ikz)/k`.
Thus the full projection error w satisfies
`||w||<=D_i/K`, `||curl w||<=E_i`. Real-part norms are bounded by these
complex norms. All constants are finite integrals of the selected field and
envelope. In dimensional coordinates the same argument uses physical k
and a fixed reference inverse length in place of the shorthand K>=1.

## Review-derived appendix B: internal core jet and a finite positive bound

Use a z-independent compact cage bump phi supported where omega_z has one
strict sign, and the circular generators

    xi_i=-curl(phi p_i)/k,
    p1=(cos(kz),sin(kz),0), p2=(-sin(kz),cos(kz),0).

Take k>0 since lambda>0. Their exact KKS coefficient is

    B_cage=rho(1-lambda/(2k)) integral_D omega_z phi².

There is no axial-cutoff remainder because phi_z=0. Choose k>lambda and
orient the cage so the displayed integral is positive. This makes B_cage
strictly positive at a finite harmonic.

The internal core generator is
`G(x,y) cos(kz)`, `G=chi_core(r) ez cross x`, with radial chi_core equal to
one near the axis and supported disjointly from the cage. It is exactly
divergence free. Its force amplitude is the real vector `F_R=G cross omega`.
The principal Leray field is `Pi F_R cos(kz)`. Its self-helicity is zero
pointwise: its fast curl is perpendicular to Pi F_R, while its slow curl
is axial. Its principal kinetic norm is the finite k-independent number

    A_R=(L/2) integral_S |Pi F_R|².

Every principal core/cage cross term vanishes pointwise by disjoint
transverse supports. In particular there is no O(k) mixed helicity hidden
in the attachment. This is the relevant proof for the k-dependent core;
0045's *fixed* core-generator bound is not substituted for it.

Let `A_c=integral_D (phi omega_z)²>0`. For a real coefficient vector t of
the attached core/cage pair, the principal velocity obeys

    ||v0||²=A_c |t|²+A_R t1²,
    integral v0.curl v0=-k A_c |t|².

Put M=sqrt(A_c+A_R). Choose a finite G0 such that
`||curl v0||<=(k M+G0)|t|`; the L² norms of transverse derivatives of its
fixed amplitudes provide one. Sum the constants D_i,E_i in appendix A for
the two complete amplitudes, including the core in the first, obtaining
D,E. Then `||w||<=D|t|/k`, `||curl w||<=E|t|`. Direct expansion gives

    |H-H0| <= rho C |t|²,
    C=2 M D+D²+(D M+D G0+M E+D E)/lambda.

Since
`H0=rho[(1+k/lambda) A_c |t|²+A_R t1²]`, any admissible finite harmonic

    k>max(1,lambda,lambda C/A_c)

gives a strictly positive full internal Hessian. The attached core changes
the KKS pairing by exactly zero, again by disjoint supports. All terms of
the finite-k H and B remain their exact integral values; this inequality
proves their signs without replacing them by leading coefficients.

For the common momentum partner use a single cage generator at k2. The
same proof with A_R=0 gives
`H(eta2,eta2)>=rho[(1+k2/lambda) A_c-C2]`. Axial orthogonality gives
`H(eta0,eta2)=0` and `Omega(K,eta2)=0`. Therefore the sufficient bound

    k2>max(1,lambda,lambda [C2+|H(eta0,eta0)|/rho]/A_c)

makes `H(eta0+eta2,eta0+eta2)>0` without changing B_body. Choose the internal
k1 above its own finite bound and with |k1| different from |k2|; infinitely
many allowed axial harmonics remain. Orthogonality of these Fourier
subspaces holds for both H and Omega because the background, domain, and
Neumann boundary condition are z-independent. This proves the claimed
block decomposition and both strictly positive intrinsic inertias.

## Final review result

**Established as stated, with the explicit support choice and review-derived
appendices above included as part of the proof.** The finite-cylinder
four-generator Euler action supplies the positive common-angle inertia,
positive physical relative-angle inertia, and restoring coefficient from
the same energy and KKS integrals. Its covariant physical-section field map
is valid. No postulated body mass, frozen external cage reaction, or
unproved high-harmonic boundary transfer remains inside this statement.

The full parent objective still requires its declared continuum spatial
coupling and EPS-ambient compatibility to be joined to this same action;
this unit does not claim that joining by itself.
