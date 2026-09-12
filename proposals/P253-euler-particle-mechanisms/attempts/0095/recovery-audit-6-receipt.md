# Recovery Audit 6 receipt

## Model provenance and quarantine boundary

The coordinator independently verified the preserved rollout record at
`/home/dan/.codex/sessions/2026/09/06/rollout-2026-09-06T12-49-28-01a07656-64cf-74e2-ac94-9ba86652daef.jsonl`,
line 16505, ordinal 16503, timestamp `2026-09-07T03:30:17.440Z`: model
`gpt-5.6-sol`, reasoning effort `high`, approval policy `never`, sandbox policy
`danger-full-access`, cwd `/tmp/substrate-particles`.  All intervening Luna-low
text and the old Recovery-3 `12/12` receipt were treated as scratch rather than
authority.  The old files remain append-only provenance.

The pre-recovery scientific byte boundary retained for comparison was:

```text
abed1b5443c750264da1050b4f225971552f97bea76edbceca86168b2d589167  derivation.md
853a7441b18ce4eb642c79c2b0e4173b1c0eb7b31a473050ad19ac2e237db4d2  source-audit.md
b9aff7e76c5334fc908c0fc9529e59906ed5a6efce452efaac39fbd8fa182c4b  p2_principal.py
3370705e3dca3bd72ff36af547269f1b8188928a08e8e96e5ba6c0610b76ecc3  verify_p2_principal.py
```

## Independent source and calculus rederivation

The cached Hattori--Fukumoto PDF was reread at SHA-256
`c6a35c44d55d26b8bb3e8fc9e55b29b638d0983b646b17eb522a4192b6444fb9`.
Equations (7)--(9) give `V=A'-A^2-BC`; equations (22)--(24) give the smooth-
center paired exponent `15/256`.  The post-equation-(20) warning was retained:
their approximate wavevector can fail to return and saturate the episode.
The replacement is an independently derived exact finite-Cao first integral
`phi=F(P)+ell theta`, so `k=dphi` obeys covector transport and returns exactly.
The `O(delta)` finite-direction change enters the diagonal detuning, which is
tuned by `F'`; the paired coefficient changes only at `O(delta^2)`.  Every
authoritative `255/1024` interpretation was removed; surviving occurrences are
explicitly labeled withdrawn/quarantined.

The six Recovery-6 bridges were rederived as follows.

1. **Packet scaling:** `xi_N=N^-1 d exp(iNphi)` has vorticity amplitude order
   one and velocity amplitude order `N^-1`; normalization is performed in the
   same global `H^s` graph norm at input and the local `H^s` observation at
   output.
2. **Coupled Maxwell:** Gauss fixes the longitudinal field at order `N^-1`;
   the constrained transverse block has a strict subluminal characteristic
   gap on this phase.  An explicit order-`N^-1` Sylvester transform cancels
   both slow-to-Maxwell and Maxwell-to-slow order-zero rows, with the remainder
   bounded from the weighted global graph core to the local observation.
3. **Same norm:** the exact nonlinear path is chosen two derivatives smoother
   after the finite packet frequency is fixed.  The difference quotient obeys
   the quasilinear remainder equation in the same local `H^s` target; no
   one-derivative-weaker solution-map claim is used.
4. **Exact Gauss path:** tag is transported by the same volume map and
   `E_tau=E_g-(g/epsilon_EM)grad(-Delta)^-1(chi_tau-chi_g)` enforces Gauss and
   zero charge exactly.
5. **Finite rows:** every linear frozen row is axisymmetric and vanishes by the
   packet's nonzero toroidal character.  There is no linear low-frequency
   restoration.  Smooth compact generators are graph-dense by the definition
   of the selected DA tangent; independence of the remaining continuous rows
   therefore makes their smooth-core map onto the finite row space.  Their IFT
   corrects only the exact path's `O(tau^2)` nonlinear defects.
6. **Weighted global to local:** `<x>^(2sigma)` is an `A_2` weight for
   `|sigma|<3/2`, giving weighted Leray/Riesz/Hodge boundedness.  The choice
   `sigma>1/2` separates the affine Coulomb monopole from the zero-charge
   dipole tangent.  Compact support and local restriction then give the
   fixed-time global-input/local-output estimate without propagating an
   instantaneous positive weight.

The charged conclusion is only: for each fixed sufficiently thin `delta`
there is `g_0(delta)>0`, and the returned gap persists for signed
`0<|g|<g_0(delta)` satisfying
`C_mon(delta)g^2<c_HF delta`.  No delta-uniform bound on `C_mon` was proved.

The failure-derived complementary route was also executed to its exact current
boundary: the charged straight-column Maxwell/radial equations and constrained
slow block are derived; the resonance angle, first paired coefficient
`b(tau)`, its near-zero nonvanishing consequence, and the ordered `d=2`
coefficient at isolated zeros are explicit.  A delta-uniform charged
column/free-boundary branch and first curvature graph jet remain the precise
construction dependency.

## Post-recovery core hashes

```text
1160396b65033db08b23fa7eb031c42e87007b9d05f8f547c0149de9c5287878  derivation.md
88dd15f73885dda6e6e1bc668ed2d29f9e165b2f674ba01994679d57262f7110  source-audit.md
98f0b73161825287abee443fc58970c6923dff16df2be88a249f74dcac190661  result.yaml
7f667706c4e22c9e36126466ef4f166ae0588e07920dac13b5aeddc10a49b870  validation.md
1d68d9634e5564e40c8748bfed7a0f879cb24af255a1cec45dd13352790cd502  p2_principal.py
ecac29802577e5f1d1e29eaaa874ab69f3909689df4136c0e0859e62a21bd66a  verify_p2_principal.py
c22eb7847075c4410f37c69cea4d336c342ed7610c83c0adcd3239ee624b97f9  src/substrate_framework/euler_p2_principal.py
c9e6fcce3498211efc709dc14f8890ce85cb446def1cb4da91e730e7b89ceea7  tests/test_euler_p2_principal.py
```

The verifier reported `14/14` exact checks, the focused API suite reported
`3 passed`, both exited `0`, and both used the repository interpreter.  Their
first successful outputs are preserved in the Recovery-6 receipts.  The two
preceding setup/structural-comparison failures are exposed separately in
`packaging-test-failures.md`.
