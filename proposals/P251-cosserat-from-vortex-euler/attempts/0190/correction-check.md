# One correction check: 0191 resolves the transverse supplier finding

Reviewer: `/root/smooth_core_review`. This is the single correction check
following `0190/review.md`, restricted to 0191's new pairing,
normalization, tensor and directly affected local curvature representative.
0187/0181 remain unchanged; their prior proof and verification receipts
are reused. I originally derived the supplier-frame correction during
the independent review, as disclosed there; this check is not advertised
as another independent scientific pass on that same correction.

**Verdict: established as corrected.** The positive inherited-energy
column-mixture theorem survives, and the concrete geometric finding is
resolved. No load-bearing correction remains at this attachment boundary.
This verdict does not imply the parent same-field EPS/acoustic objective.

## Actual time reversal, phase, energy and current

I read the complete corrected proof, verifier and first-run output.
Writing `J=[[0,1],[-1,0]]`, the physical circular generator is
`A_s=s nu J` and the supplied form is `Omega_s=s M_raw nu J`.
Consequently `-Omega_s A_s=M_raw nu^2 I>0`: its sign is inherited from
the actual reversed mode, not fixed by taking an absolute value.
For `E_s=diag(A,s A/nu)`, direct multiplication gives

    E_s^T Omega_s E_s=M_raw A^2 J,
    E_s^T H_s E_s=M_raw A^2 diag(nu^2,1).

Both partners therefore have the same positive initial scalar phase and
physical energy. Equal half-weights neither cancel nor double these
quantities. The full actual tilt is
`Theta_s=n theta+s b theta_t/nu`; differentiating and using
`theta_tt=-nu^2 theta` gives its stated circular generator. The
conjugate component cancels between realizations, including parameter
derivatives. Since the actual 0181 rows are vector identities
`G_s=Delta_raw Theta_s,S_s=Delta_raw Theta_s,t`, the same cancellation
applies to both physical currents and retains the initial displacement
moment. This implements, rather than presumes, the required scalar
preparation on the transverse circular modes.

## Correct common-K normalization

The carrier shifts along `t`, whereas the observed scalar tilt axis is
`n`, with `n.t=0`. With common vector initial data projected on `n`,
the reconstructed observation is `3 E[n theta]`. The inherited energy
has no extra factor three. Thus its raw mass and literal current each
acquire exactly the Haar factor `1/3`:

    j=M_raw/3, Delta=Delta_raw/3.

Conditioning on `t` gives `E[nn^T|t]=(I-tt^T)/2`. Contracting with
the isotropic second/fourth moments independently yields

    3 E[nn^T(t.K)^2]=(2|K|^2 I-KK^T)/5.

The corrected factors `c_T=2/5,c_L=1/5` consequently enter the observed
histories, phase, physical energy and current together. The original
positive scalar variance/energy equations have the same common factor
in each channel, so their already derived two-slope solution requires
no change. Both curvatures remain positive. The statement continues to
describe positive probabilities with signed initial optical amplitudes,
not a signed probability law or a claim of universal intrinsic modulus.

## Directly affected curvature representative

The new local representative is also correct with the canonical
`micropolar.py` convention

    W=c_tr (tr G)^2+c_s |sym G|^2+c_a |skew G|^2.

For `d0=4C_T+3C_L`, its coefficients satisfy

    c_s+c_a=C_T, 2(c_s+c_tr)=C_L,
    c_s>0, c_a>0, 3c_tr+c_s=9C_L^2/(2d0)>0.

Thus all irreducible pointwise sectors are positive for every
`C_T,C_L>0`, including the corrected transverse/longitudinal ratio.
This removes reliance on the older optional equal-symmetric/skew
choice. For two representatives with the same bulk coefficients, let
`eta_N=c_s,old-c_s,new`. Their energy difference is exactly

    eta_N[(tr G)^2-tr(G^2)],

because their coefficient changes are `(eta_N,-eta_N,eta_N)` and
`tr(G^2)=|sym G|^2-|skew G|^2`. Its displayed divergence is therefore
the actual null-Lagrangian flux, not an energy term erased at a boundary.
The proof correctly scopes the `C_T,C_L` values from 0172/0182 to a
future join supplying those hypotheses; the optional canonical current
improvement does not assert equality of physical and canonical position.

## Receipt

The 12-check first-pass receipt was read and reused. It includes direct
pullbacks for both reversed partners, the actual generator, the new
Haar tensor, the canonical Fourier stiffness comparison, and the
coefficient/null-flux identities. The weaker trace-only part of the
one-third check is supplemented by the explicit action/current
contraction above. No full replay or new numerical design is needed
for this bounded exact correction.

Verified SHA256:

| Artifact | SHA256 |
| --- | --- |
| 0191/transverse-family-correction.md | 487bbef625e29180e602d04f87bd162dcef35b39aecf6183acabb96ebba7b165 |
| 0191/verify.py | edd8ac07bec3ff1432890dec792c27faf11cb98ea8883826786d289fb4964fff |
| 0191/first.stdout | 35a129b1522418a9dc6759e7ad02a873003718277786d591d2496cca0b0a41fe |

The 0190 review transaction is complete. The strongest supported result
is the positive actual smooth-column mixture, observed closure and
inherited physical energy/current through spatial order two on fixed
time windows, with the corrected transverse geometry. The previously
declared finite-packet, same-cell EPS and acoustic continuation work
remains with its respective active constructions.
