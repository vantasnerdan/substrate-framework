# Independent review of the positive correlated displacement preparation

Signed reviewer: `/root/smooth_core_review`.

**Verdict: established as stated at the D-only boundary.** The actual
finite periodic Euler preparation in 0188 has a positive whole-law
physical displacement-response coefficient, matching complete inherited
initial D energy and the stated initial fluid phase. No load-bearing
correction was found. Its common-V and physical optical/current closure
are explicitly still open, rather than silently inferred from this
positive result. This review makes no parent EPS/Cosserat completion claim.

## Independence, inputs and exact scope

I did not author or implement 0188. I authored unchanged 0170 Euler/Lin
observation identities and 0179 velocity/current calculations. Their
role here is disclosed: they are source interfaces, not additional
independent evidence passes supplied by this review. After the 0197
registration and central/schema confirmation, I read the complete
positive proof, both verifiers, the original degree-seven execution,
the phase/energy and stationary-corrector proofs, the nonlinear-return
and velocity-partner discussion, and the recorded velocity-lift repair.
I checked the actual 0154/0180 action definitions and the canonical
Fourier operations at their use sites.

The proposition reviewed is existence of this explicit displacement
preparation and its fixed-time second spatial jet, energy, and initial
phase—not existence of a closed scalar acoustic phase dynamics. The
whole law rotates one common laboratory `(kappa,D)` with the actual
field and preparation. Its probabilities remain positive. Its signed
stationary velocity coefficients are microscopic initial correlations,
not negative probabilities or a renamed density.

## 1. Full Euler pressure and the selected lift

Let `T_D=-(D.grad)u`, `lambda=-1` and
`q_D=kappa cross T_D/lambda`. The field has `curl u=lambda u`, and
the actual generator is

    L_K w=P_K[u cross (curl_K w-lambda w)].

The divergence-compatible lift is exactly
`T_K=curl_K(T_D/lambda)=T_D+ik q_D`. Its first-order defect is

    F_D=P[u cross curl q_D]
       =-P[u cross (kappa.grad T_D)]/lambda.

For the actual unequal two-wave field, differentiating each circular
wave twice gives `kappa.grad T_D=kappa_Y D_Y u_Y+kappa_Z D_Z u_Z`.
It follows directly that

    F_D=d P(u_Z cross u_Y), d=kappa_Y D_Y-kappa_Z D_Z.

This independently identifies the forcing used by the degree-seven
verifier. It is the full Euler forcing, not merely a planar curl row.
The exact constructed `z_*` solves `L z_*+F0=0`; its axial component
`-alpha/2+Hphi` and the mean-preserving Leray projection are present.

On a unit Fourier shell, `Pi_-=(P-curl)/2` is the genuine orthogonal
negative-helicity projection. Consequently every returned field in
0188 equation (3) is solenoidal, mean zero and killed by the FULL
linear Euler generator: `curl z_return=lambda z_return` makes the
displayed generator identically zero. Arbitrary correlations with the
whole `(kappa,D)` do not alter that exact property. Additional first-shell
wave directions absent from both source and observation are orthogonal
and unnecessary; the projection used here does not drop a driven row.

Thus `z_D=d z_*+z_return` is an actual stationary first-cell velocity
corrector. No stationary material displacement is assumed. In particular
the earlier Kelvin-versus-curl-lift discrepancy in the V equation does
not enter this D construction; the supplied D lift above is the actual
one used in its code and energy.

## 2. Finite positive range construction

At the endpoint A=0, the Chebyshev moments have the stated values:
ordinary cosine orthogonality fixes `c1=-1`, while
`<T_n'(cos Y)>=n` for odd n gives `sum n c_n=2` from the other
desingularized first-shell constraint. For an individual harmonic,
`phi=c_n cos(nY)/(n^2-1)` and its actual range velocity norm is

    c_n^2 n^2(n^2+1)/(2(n^2-1)^2)=w_n c_n^2.

The coefficients `c_n=2n/(w_n S)` therefore satisfy the constraint
and give the stated finite cost. The two moment columns for c1,c3
have endpoint matrix `[[1/2,0],[1/2,3/2]]`, so the continuation
does not assume a nonexistent inverse at equal wave amplitudes.

More importantly, the positive construction actually solves the finite
moment equations at `A=1/100`, rather than relying on an unspecified
smallness interval from that endpoint. The exact code evaluates the
Chebyshev polynomials by full Fourier convolution, removes the actual
first shell, inverts each remaining nonzero `H-1` coefficient, and
checks all three components of `Lz_*+F0`. Odd simultaneous phase parity
also excludes an unrecorded constant term. Its archived first execution
gives the displayed rational R and the exact strict inequality
`0<R<13E/1280`. No rounded eigenvalue, Fourier truncation of a generated
product, or empirical comparison is involved. That existing 12-check
receipt is a valid exact finite calculation and was reused.

## 3. Inherited energy and physical response are independently computed

For clarity let `A_K=u.grad_K`, `C_K=P_K A_K`,
`N_K=(I-P_K)A_K`. The constrained momentum and Hamiltonian are

    pi=rho(eta_t+C_K eta),
    H=||pi||^2/(2rho)-<pi,C_K eta>
                     +rho<eta,(Hess p-N_K^*N_K)eta>/2.

Substituting the ACTUAL initial material rate, rather than discarding
the normal-pressure term, gives

    H=rho[||eta_t||^2-||A_K eta||^2+<eta,Hess p eta>]/2.

At `eta_D(0)=D`, write
`b_final=q_D-(kappa.u)D+z_D`. Its zero-order divergence vanishes,
and the exact finite-K rate is `ik P_K b_final`. Also
`A_K D=ik(kappa.u)D` and `<D.Hess(p)D>=0` by periodicity. Hence the
complete second coefficient of the initial D energy is indeed

    rho k^2[||b_final||^2-||(kappa.u)D||^2]/2.

The subtraction is a derived Euler/Jacobi Hamiltonian contribution,
not a positive norm relabeled as an elastic energy.

Separately, the physical Euler velocity mean obeys

    m_t=-ik P_kappa<(kappa.u)w+u(kappa.w)>.

For `w_D=T_D+ik(q_D+z_D)+O_T(k^2)`, the zeroth translation term has
zero mean stress by periodic differentiation. The stationary corrector
therefore makes the second coefficient constant. Contracting it with
the unit transverse D gives

    R_D=<(kappa.u)^2>+<c,b_final>,
    c=(kappa.u)D+(u.D)kappa.

Here `kappa.D=0` is used explicitly. The physical observation is
`X_D=D+integral m_D`, as in the cited 0154/0180 definitions. It is not
the mean material displacement Y, whose integrated-current correction
remains part of the upstream Euler/Lin map.

I independently recomputed the exposed first-shell constants using
Cartesian cosine/sine vector coefficients, not 0188's Fourier/Haar
helper. The orthonormal-pair fourth tensor used is

    E[kappa_i kappa_j D_a D_b]
      =2 delta_ij delta_ab/15
                         -(delta_ia delta_jb+delta_ib delta_ja)/30.

The calculation gives

    h_shell=E(16t^2-47)/240,
    R_D=E(8t+13)/120,
    E||Pi_- c||^2=E/15.

Orthogonality of the remaining Fourier shells adds exactly R to the
first expression and nothing to the second. It follows independently
that

    h2-a=R+(120a-E)(120a+19E)/(960E), a=-R_D.

The explicit radical preparation in 0188 solves this equation, and
the strict rational R bound gives `E/240<a<E/120`. Thus actual
positive energy and actual restoring response agree; neither was
defined by the other. The review check also shows that omitting the
axial `-d alpha/2` return changes the energy constant, an exposing
mutation of this same calculation.

## 4. Actual initial phase and time scope

The finite-K D preparation has mean-zero physical velocity and

    pi_D=rho P_K[w_D+(Du)D]
        =ik rho P_K(q_D+z_D).

All microscopic Fourier coefficients are nonzero harmonics, and P_K
acts mode by mode without generating a zero harmonic. With
`eta_V(0)=0,w_V(0)=V`, one has `pi_V=rho V`. Pairing these actual
columns therefore gives `Omega_DD=Omega_VV=0,Omega_DV=rho I`.
This holds for the stated Bloch encoding, or its normalized real
conjugate pair with the usual fixed-cell nonresonant small-K choice.
There is no hidden cosine-average density factor. The preparation
changes initial circulation data, as explicitly allowed and stated;
it is not falsely identified with the Kelvin leaf of the constant D
generator. The initial physical kinetic mass of V is rho, not the
positive norm of a microscopic corrector.

The preceding actual mean equation now gives
`X_D(T)=D-a k^2 T^2 D/2+O_T(k^3)`. Smooth finite periodic initial data
and the full Euler/Lin finite-time estimates support this fixed-time
spatial expansion. It does not imply a uniform acoustic-time estimate.
Nor does the initial rho phase plus this D history fix the V history,
the physical Wronskian, or an autonomous observed Hamiltonian. The
proof retains the actual V response equation and its possible current
closure instead of making any such inference. This is the appropriate
positive attachment boundary, not a missing condition of its theorem.

## Evidence receipt and disposition

`verify_review.py` passed all 9 exact checks on its first execution;
`first.stdout` preserves that output and exit zero. Its first-shell
oracle uses an independent Cartesian phase representation and explicit
fourth moment, while reusing only the canonical CheckLedger. Targeted
Ruff and attempt-scoped diff checks pass. Existing source receipts are
reused; no full validation or numerical spectral design was performed.

Verified SHA256:

| Artifact | SHA256 |
| --- | --- |
| 0188/positive-correlated-preparation.md | ebb26daabdc86d8ae3a8167ea2030404c79690da5b751f8fb34d658105e4dee8 |
| 0188/initial-phase-energy.md | ad68f642da3a3ce85a9867f2de31fc90a20502102e773baf06bedcc9af30338b |
| 0188/verify_chebyshev.py | 284042c053e588d2893b6fd06122690906d5900f422bfdd69d05e944936dc81f |
| 0188/verify_stationary.py | e3c5c3c9314c44c1585d7b5d598a38232336c0da989d092252b42689bb2e9ea0 |
| 0188/chebyshev-first.stdout | d0f82c6a1fbcd3a426f7769e3a322fb75a507d8c0a40bf923746d10b4e94f7c7 |
| 0197/verify_review.py | a5baade595f82dcb32d21ea4ec26663a7bec903b180c1c668e2674aa8d6d2342 |
| 0197/first.stdout | 2606974cef34482d05500de9d35137257132c0e9bf1c94b7a8e5479042bf4e70 |

The single substantive review is complete with a positive verdict and
no correction request. The parent campaign remains active on its own
common-V, optical/current, same-field EPS and coupled continuum work;
none is claimed closed by this D-only result.
