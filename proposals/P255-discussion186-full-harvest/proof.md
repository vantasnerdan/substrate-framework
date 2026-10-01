# P255 — independently checked results from discussion #186

This supplement to [P252](../P252-newton-sign-external-audit/audit.md) and [P254](../P254-m5-strand-scaling/proof.md) covers the **whole original** 2026-09-25 snapshot (74 top-level comments and eight replies). A [75th R26 comment](r26-addendum.md) arrived during PR delivery and was indexed and evaluated separately without broadening these six frozen theorem statements. The full Markdown, paginated GraphQL responses, and link inventory are saved locally; [comment-index.md](comment-index.md) and [script-references.md](script-references.md) expose public comment/link occurrences. This proof promotes only six narrow, independently falsifiable identities. Each action, field class, norm and source revision is declared; no M5 electron, physical charge, Newton force, continuum halo or stable full-field strand is promoted.

## C-M5H-001 — zero-cost spatial orbit, but *not* a global no-go

For the OpenWave R13-W *certified* action `L_cert=-4 I1-V4`, take `eta=diag(-1,1,1,1)`, `0<delta<1`, `D=diag(8,1,delta,0)`, and `M=R12(psi) D R12(psi)^T` with a **spatial** `(1,2)` rotation. The four trace powers of `M eta` equal `(-8)^p+1+delta^p`, so the vacuum trace potential vanishes. Every spatial jet is `A_i=psi_i M_psi`, hence every spatial commutator `F_ij=A_i eta A_j-A_j eta A_i` vanishes. For the specifically declared rigid `(2,3)` generator `G`, `A_0=[G,M]`; direct 4×4 multiplication gives the positive kinetic factor

```
4 tr(eta F_{0i} eta F_{0i}^T) / psi_i²
 = 8(1-delta)² [1-(1-delta²) cos² psi] > 0.
```

Choose `psi_w=Psi phi(z/w) chi(sqrt(x²+y²)/R)` for smooth compact bumps with nonzero longitudinal derivative, fixed `R,Psi`, `w>0`. Each is smooth and vacuum outside compact support. Its static spatial-curvature-plus-trace energy is exactly zero; its rigid-generator inertia `I_w` contains a strictly positive longitudinal integral `c/w` and hence diverges. For the **declared family and formal fixed-`J` Routhian** `E_J=E_stat+J²/(4I_w)` with `J!=0`, the infimum is zero and no finite-inertia member attains it. This family does **not** prove that the full eta-contracted static functional is globally nonnegative or that a global minimizer does not exist: boost directions can have negative signs. **Correction after review:** global constant spatial rotations preserve the unrestricted certified action, so `(2,3)` *is* a compact action symmetry. It rotates the nondegenerate exterior `D`, however; a time-dependent rigid rotation of this orbit is not admissible under the fixed-exterior boundary condition. The displayed formal `J` is not a proved Noether charge of a finite-energy fixed-boundary dynamical solution. Source [R13-W](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18256905), pinned OpenWave `0e487065` (`m5_32_r13w_w0.py`, `m5_32_r13w_common.py`). Independent 4×4 orbit calculation and positive-factor checks: `tests/test_m5_discussion186_harvest.py`.

## C-M5H-002 — R17 pseudoscalar-square basis identity

At OpenWave `1848124faf40861c4024521530721a82e889da30`, in the **R0 mixed eta/delta slot convention**, let four independent symmetric 4×4 jets `A_m` give `F_mn=A_m eta A_n-A_n eta A_m`. The source's contractions are

```
I1 = sum_{m<n} eta_m eta_n tr(eta F_mn eta F_mn^T),
I2 = sum_{m,n,a,b} F_mn[a,b] F_ab[m,n],
I3 = sum_{m,n,a,b} eta_m eta_b F_mn[a,b] F_ma[n,b],
X  = (1/2) sum_{m,n,a,b} epsilon_mnab eta_m eta_n F_mn[a,b].
```

An exact expansion over all **40 independent jet components** gives `X²=-2I1-I2+4I3`. The test recomputes both sides directly from the defining matrices and rejects coefficient `3` instead of `4`. Thus adding a *constant-coefficient* `X²` is a change of coefficients in the existing R0 quartic basis, not a new independent invariant. Field-dependent coefficients, a different index convention, and a positive/stable clock are not covered. The source `m5_32_r17_0_symbolic.py` [R17](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18358755) reports numerical least-squares span; this proof uses exact arithmetic instead. Its separately stated one-epsilon `X` divergence is not silently substituted for the square.

## C-M5H-003 — full spatial frame quartic and uniaxial limit

For the **positive spatial sector** of the OpenWave R22 eta action, set `M_sp=R diag(lambda_0,lambda_1,lambda_2)R^T`, where `R(x)∈SO(3)` is smooth and the three eigenvalues are **spatially constant**. Let `hat(W_i)=R^T ∂_iR`, `A_i=∂_i M_sp`, `F_ij=[A_i,A_j]` and `G_a=∏_{b≠a}(lambda_a-lambda_b)`. Then

```
4 sum_{i<j} tr(F_ij F_ij^T)
 = 8 sum_a G_a² sum_{i<j} [(W_i × W_j)^a]².
```

In the moving frame `A_i=[hat(W_i),diag(lambda)]`; expand the antisymmetric commutator and its Frobenius norm. The independent exact symbolic test keeps all three eigenvalues and six frame-rate components free and rejects the **opposite-gap** mutant. For `(1,delta,delta)`, only `G_0²=(1-delta)^4` remains: the uniaxial Faber identity is `u=8(1-delta)^4 sum_{i<j} [n·(∂_i n×∂_j n)]²` in this stack convention. The other author's sum-ordered convention sometimes reports `4 lambda^4|E_top|²`; compare definitions before coefficients. Local varying eigenvalues, time-space curvature, nonorientable frames across degeneracies, pair relaxation and a free-space `1/d` law are **not** proved. Pinned [R22](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18519432), OpenWave `f92ae6afa560e544eeb1fc5571c939c652acb5d4`, `m5_32_r22_0b_frame_identity.py`.

## C-M5H-004 — symmetric potential valley slope only at its stationary point

Let `V(lambda_1,...,lambda_4)` be `C²` and permutation invariant. At a fully symmetric stationary point `lambda_a=lambda_*`, write diagonal Hessian entry `A` and off-diagonal entry `B`. On `lambda=(nu+Delta,nu,nu,nu)`, the derivative of the `nu`-stationarity equation has `V_{nu,nu}=4(A+3B)` and `V_{nu,Delta}=A+3B`. If `A+3B!=0`, the implicit function theorem gives `nu'(0)=-1/4`. A nonseparable symmetric quartic test exercises this result; the same calculation for weights `(2,1,1,1)` gives `-2/5` and rejects dropping permutation symmetry. This does **not** imply slope `-1/4` at a nonstationary hierarchy point, for a nonsmooth ranked-eigenvalue maximum, or without the nondegenerate Hessian. Source [comment 18437206](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18437206); no linked executable for its named `claim_865_symmetric.py`.

## C-M5H-005 — one-director `E·B` obstruction, no all-model claim

For a smooth unit three-component director `n` on 3+1 spacetime, its canonical two-form `H_{mu nu}=n·(∂_mu n×∂_nu n)` is the pullback of the two-dimensional sphere's area form. Its tangent jets have only two independent components: `H_{mu nu}=u_mu v_nu-v_mu u_nu`. Thus `H∧H=0`, equivalently `H_{01}H_{23}-H_{02}H_{13}+H_{03}H_{12}=0`, or `E·B=0` in any consistent orientation convention. Exact symbolic four-jet test checks cancellation; **two independent** decomposable two-forms with nonoverlapping planes have a nonzero sum wedge, rejecting the extension to a sum of directors. Boost and eigenvalue-gradient terms, higher target dimension and full M5 curvature are outside the theorem. Source [v2.8/v2.9 exchange](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18557220), v2.9 [attached ZIP](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18557744) `three_forms.py`/`e_dot_b.py` (artifact, not a Git commit); both scripts exited zero on the downloaded archive and exhibited a two-form counterexample.

## C-M5H-006 — CP² projected Gram curvature, not field-only Maxwell

Take a unit complex `z∈C³`, rank-one projector `P=zz†`, horizontal derivatives `u_i=(I-P)∂_i z`, and `G_ij=u_i†u_j`. Put `F_ij=2 Im G_ij`. Direct block multiplication in the gauge `z=e_0`, `u_i⊥z`, gives

```
||[∂_iP,∂_jP]||_F²
 = 2 (Re G_ii Re G_jj - (Re G_ij)²) + (3/2) F_ij².
```

The first (real Gram determinant) term is necessary: a deterministic complex `C³` tangent counterexample yields a strictly nonzero residual if it is omitted. The test compares independent projector-matrix multiplication with this formula, not two copies of one expression. It is gauge invariant and geometric, not an SU(3) confinement or nuclear spectrum result. Source [v2.9 ZIP](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18557744) `cfn_identity.py`, executed with max reported identity residual `2.7e-14` on its own numerical sweep; the test's integer-complex example supplies independent exact evidence.

## Reproduction, explicitly excluding headline physics

`PYTHONPATH=src .venv/bin/python -m pytest tests/test_m5_discussion186_harvest.py -q` exercises the independent exact oracles and plausible wrong-coefficient/non-symmetric countercases. The source-linked scripts are provenance and independent execution where recorded in [harvest-decision.md](harvest-decision.md); a script's exit code or author checksum does not make its solver, asymptotic or physical inference an accepted result. [P254](../P254-m5-strand-scaling/proof.md) remains the sole accepted fixed-slot radial theorem for the **two separately normalized** Report-018/OpenWave R25 potentials.
