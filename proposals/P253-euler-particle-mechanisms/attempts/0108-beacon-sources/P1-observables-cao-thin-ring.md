# P1 exact observables on ONE carrier — fixed sufficiently thin Cao ring (beacon 0108)

Carrier: one fixed member of the Cao thin-ring family at fixed small core parameter `delta`,
uncharged base state `omega_g` with tag `chi_g` (0095/0104 scope, T9). Charged extension only
at signed `0 < |g| < g0(delta)` under `C_mon(delta)·g² < c_HF·delta` (0104 Unit H, T9).
No delta-uniform wedge; no complementary-regime claim.

## Exact observable inventory (all on this carrier + same state space)

| Observable | Exact statement on this carrier | Source (tool-cited) |
|---|---|---|
| Axial impulse | `I_z = pi·rho·kappa·R²`, physical (not velocity-integral) momentum row | T8 (`euler_cao_schur.py:53,87`) |
| Translation speed | `c = kappa·L/(4·pi·R)`, `L = log_inverse_core` | T8 (`euler_cao_schur.py:52,86`) |
| Chemical potential | `mu = 3/8·kappa·R·L/pi` | T8 (`euler_cao_schur.py:51,85`) |
| Parameter/moment Jacobians | `det d(mu,c)/d(k,r)`, `det d(kappa,I_z)/d(k,r)`, Schur `S = moment/parameter`, responses `dc/dr|…`, `dI/dr / dc/dr` | T8 (`euler_cao_schur.py:60–92`); module asserts only the source-derived leading jet, not the Fredholm theorem or charged branch |
| Far-field velocity | Dipole `[3r(r·I)−|r|²I]/(4·pi·|r|^5)` + absolute remainder `O(∫|y|²|ω|/|r|^4)` for `|r| > 2·core` | T8 (`euler_impulse.py:23–36`) |
| Pair cross energy | `rho·I_b·u_a(d)` at dipole order (leading kinetic cross term, not a force law) | T8 (`euler_impulse.py:39–51`; docstring: force needs a derived collective action) |
| Similarity weights | `u_AB(x,t) = A·u(B·x, A·B·t)`: `E ∼ a²/b³`, `I ∼ a/b³`, `Gamma ∼ a/b`, `H_helicity ∼ a²/b²`, topological charge invariant | T8 (`euler_scale_causality.py:46–62`) |
| Pressure far field | `p = rho·M_ij·∂_i∂_j(1/4·pi·|x|)`, `M_ij = ∫u_i·u_j`, leading term only | T8 (`euler_scale_causality.py:65–98`) |
| Linear observed growth | Fixed `J_loc`, `c_obs > 0` independent of circuit count: `||J_loc·S_g(jT_*)||_ess ≥ c_obs·lambda_+^j`, weighted global input → local same-`H^s` output; tube + `N0(j)` may depend on `j` | T9 (0104 Units E/F; fixed-time Egorov `O_T(N^-1)`, no growing-time uniform theorem) |
| First-curvature coefficient | `gamma_hat,eps/(delta·s) = 15·A^p/512 + O(s²) + o_eps(1)` on the compact regular tube; resonant exponent `gamma_hat = delta·|U_theta|·15/256 + …` | T9 (0104 Units A/C) |
| DA/tag/Gauss packet | `xi_N = N^-1·chi_tube·Re[d·e^{iNphi}]`, `q_N = curl(xi_N × bar_omega)`, tag zero-mean ⇒ electric starts at dipole; weight interval `1/2 < sigma < 3/2` | T9 (0104 Unit D; linear-input scope only) |

## Norms and constraints (declared, not invented)

- Input: weighted global complete-state graph domain `X_in^s` with `<x>^{2σ}` (`1/2<σ<3/2`), compact-curl fluid zero-Fourier row, monopole exclusion.
- Output: local same-`H^s` quotient under fixed `J_loc` (fluid + tag + field rows, fixed once, independent of `j`).
- Constraints enforced in the packet: linear material + Gauss law, Leray/Hodge continuation, Bogovskii toroidal correction.
- Total-momentum row `P_tot = rho_m·∫u + eps_EM·∫E×B` is NOT a continuous functional on the
  declared `σ < 3/2` topology for generic compact mean-zero vorticity perturbations
  (`uhat = O(1)`, `|x|^-3` tail, not `L¹`) — 0104 Unit G blocked (T9). This pass adds no repair;
  0107 Route A (`J_ren`) is the interrupted draft addressing it, not a result.

## Equal-coarse-state / memory check (P1 obligation)

- 0042/0045 (T4/T11): material observables use transported bounded region / finite-moment tag;
  forced response uses finite-rank Sobolev projection + complementary range; `P_obs, Q_obs`
  distinct from velocity projection; Mori–Zwanzig identity exact only where `Q_obs·L`
  generates the declared observable evolution. Missing internal variables / memory retained —
  no autonomous closure declared.
- Choi–Jeong exposing constraint (S9): the P1 observable set must detect filamentation-type
  geometric deformation coexisting with orbital stability. Current linear observed norm tracks
  the returned-phase hyperbolic channel only; filamentation detection is an open P1/P2 test,
  not a claimed capability.

## What is NOT claimed on this carrier

Fixed-row nonlinear closure (needs Unit G), delta-uniform charge wedge, complementary charge
regime, full-3D coercive Hessian, scale selection, charge/quantum/relativistic origins,
P2 persistence, P3 interaction law, any P4–P6 identification. All remain open by citation (T9).
