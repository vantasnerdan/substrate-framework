# 0163 — frozen response/action construction receipt

Owner `/root/construction_review`, 2026-09-05. The initial route and its
generalized-profile continuation were centrally registered and validated
by root before implementation. Base v0.176.0/dbf0c04; unchanged 0161
evidence was subsequently promoted by the parent, not by this attempt.

## Route verdicts and positive result

- Constant-curl axial-positive-acoustic candidate: **refuted**, with the
  mechanism derived from the full Euler current AND same material action.
  Its actual positive acceleration gives cosh on T/|k|; the exact phase
  energy has stiffness -2rho lambda(sqrt(lambda^2+k^2)-lambda) C_v.
  This is a verdict on that axial candidate, not on its optics, other
  directions, or the parent objective.
- Generalized-force-free ordered finite-C family: **established**. For
  each finite nonzero k, finite C chosen last controls all pressure/shear
  response and phase-action errors on the acoustic interval. No fixed-C
  assertion is inferred from that comparison alone.
- Stronger generalized-force-free FIXED-C acoustic candidate:
  **established**. For every sufficiently large fixed finite C, the
  actual smooth stationary same-array field has positive axial acoustic
  mean histories on t<=T/|k|, with O(|k|/lambda) normalized error, actual
  mean axial Doppler shift, and full material-action mass rho.

`evidence_scope: exact analytic full-pressure/current/corrector theorem,
with symbolic verification of the field, projection, boundary, action
and graph-smoothing identities`

The fixed-C coefficient is computed from the same field and its actual
stationary transport/Arnold responses:

    C_W=<a v tensor v>-<h F>-<W d(cT)>
        =C_v+O(delta^2),     C_v=3Psi^2 lambda^2 I/4.

No modulus is assigned from a desired oscillator. The nonzero conserved
translation rows, both triangular C3 cells and the hexagonal C6 cell,
full mean pressure, logarithmic transport domain, and graph-norm smooth
preparation are explicit. The leading physical mass and action come
from the same actual material phase, with all finite-k connections and
current corrections retained.

## Captured exact execution

All commands used `PYTHONPATH=src` and the existing interpreter
`/home/dan/substrate-framework/.venv/bin/python` from
`/tmp/pr199-completion`.

| Script | Captured successful output | Result |
| --- | --- | --- |
| verify.py | repaired-run.stdout | 15/15, exit 0 |
| generalized_verify.py | generalized-repaired.stdout | 16/16, exit 0 |
| fixed_profile_verify.py | fixed-profile-first.stdout | 17/17, exit 0 |

The exact polynomial-source extension is captured in
fixed-source-first.stdout and fixed-source-extension.stdout. The final
`fixed_profile_source.calculate(2)` still implements the first psi^2
calculation; its main entrypoint executes the subsequently added powers
1,3,4. Those values expose the nontrivial mean-frame calculation but
do not replace the general boundary proof.

The two native verifier failures are preserved with diagnoses: a
list/matrix subtraction in the first sign check, and expanded/factored
SymPy structural equality in the generalized action check. Their repairs
changed no equation, coefficient, profile, branch or tolerance. The new
fixed-profile verifier passed on its first execution. These are analytic
and symbolic checks, not sampled stability spectra or simulation evidence.

Ruff passes for all four new Python files. Direct inspection finds no
trailing whitespace in this attempt. Only this construction and its new
verifiers were run; no full repository validation or unrelated scientific
replay was performed. The generalized verifier imports the parent's new
euler_acoustic definitions rather than duplicating its planar phase API.

## Frozen load-bearing artifacts

| Artifact | SHA-256 |
| --- | --- |
| constant-curl-response.md | `c581163f3e88b1545aab5dfdee554882309578f4fb44457544c47e1c62940b7f` |
| ordered-force-free-transfer.md | `eebe156c62499a5637eb8a08c0c7b48ac623980cc1fb83042f5cd25f234b6820` |
| fixed-profile-acoustics.md | `1f77ed0f15a1b126952a052bdbea1008aee46843337c84330a230c4af7ff419b` |
| verify.py | `6fa5bf0d4396cc80ffececde2d880731a78d363025693aa563d8a84b28d35776` |
| generalized_verify.py | `90052904f0ad6a5c9fe84baca5c98ed0166d5be064f4a2c29240e03c20fd8f7d` |
| fixed_profile_verify.py | `9942c843ec7565656b96368e7eeadd45d824955adce581debf10ed4fff250ba9` |
| fixed_profile_source.py | `f204c50e4f90d7a1edd3497f24eefabf568778d743792b7be900e3a14446b889` |

## Remaining achievements

This attempt changes no canonical claim or API and makes no promotion or
parent-completion declaration. Individual review and reusable extraction
are the next claim-level steps. The shared generalized-field optical
construction, generic-direction/whole-field acoustic law, global closed
vortex-torus embedding and full coupled continuum remain active parent
obligations. The original compensated compact-core route 0156 remains
active rather than being replaced or declared exhausted. The constant-curl
generic-direction/O(3) candidate also remains explicitly registered.
