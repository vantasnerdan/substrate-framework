# P254: exact fixed-slot straight-strand tension, not a full M5 soliton

## Provenance and scope

This is a **conditional result for two distinct proposed actions**, not a new accepted equation of the framework's auxiliary-frame M5 action. Source snapshots: OpenWave `a73f6fcd71d6b45bf479e76a3999123b39e6d9d6`, `m5_32_r25_0_form.py` (discussion #186 R25); Mikulski `4823e5ec2368420f28229ac73067108ec01e6e23`, Report 018 `strand_bogomolny.py` and Report 017 `symbolic.py`. The latter uses a field-Euclidean positive matrix norm, whereas OpenWave's stationary, positive spatial-pair reduction comes from its own eta action. The expressions below are normalized **separately**. The two competing apparent claims, cubic versus quartic splitting power, cannot be compared as calculations of one action. Both author scripts were executed at these pinned commits; their numerical row generation is not independent proof.

Fix the exterior, all other slots and the pair mean `s0`. Set `S=s0 I+b(rho)[cos(m phi) sigma_z+sin(m phi) sigma_x]`, `m=2k` a nonzero integer, `b>=0`, `b(0)=0`, `b(infinity)=b0>0`. This includes a director half-winding `k=1/2` and full-winding `k=1`. The regular core, fixed slots, no gradients along the line, and positive **static transverse** curvature term are hypotheses, not consequences of either unrestricted 4D action. Energy is per unit length above the specified exterior vacuum. No mass calibration or pair interaction is inferred.

## Matrices, potentials, and normalization

Write `R=cos(m phi) sigma_z+sin(m phi) sigma_x`. Direct matrix multiplication gives `[S_rho,S_phi]=2m b b' J`, where `J=[[0,1],[-1,0]]`. Thus `F_xy=[S_rho,S_phi]/rho`, `Tr(F_xy F_xy^T)=8m^2(bb'/rho)^2`. Report 018's `2 Tr(F_xy F_xy^T)` has coefficient `A=16`; OpenWave's `4 Tr(F_xy F_xy^T)` on this static chart has `A=32`. The 4x4 eta embedding has the same spatial block when the other slots are fixed. These coefficients and contractions must not be interchanged.

For Report 018, the fixed exterior eigenvalues give `V_spec(b)=2(b0-b)^2`, `b0=beta/2`, and the reported 3x3 spectral model uses that potential as an input. For OpenWave let `q=b^2`, `q0=b0^2` and keep the other spectrum fixed. Direct differences of the mixed-matrix power traces `Delta_p=Tr(N^p)-c_p` for `p=1,2,3,4` are

```
Delta_1=0,
Delta_2=2(q-q0),
Delta_3=6 s0(q-q0),
Delta_4=(q-q0)[12 s0^2+2(q+q0)].
```

The R25 potential is `V_tr(q)=w(q0-q)^2 K(q)` with `w>0` and `K(q)=4+36 s0^2+[12 s0^2+2(q+q0)]^2`. In the R25 exterior `s0=b0=delta/2`. `K(q)>K_lead=4+36 s0^2+144 s0^4` for `q0>0`, `0<=q<=q0`: the earlier constant-`K_lead` expression is a **strict lower bound**, not the exact finite-splitting tension. The exact exterior stiffness is `K(q0)=4+36s0^2+(12s0^2+4q0)^2`, giving a corresponding upper comparison by constant stiffness.

## Square completion and sharp radial infimum

Put `s=rho^2/2` and `q=b^2`. Then `(bb'/rho)=q_s/2` and both **restricted** functionals are

```
T=2 pi integral_0^infinity [A m^2 q_s^2/4 + V(q)] ds.
```

Assume `A>0`, `q>=0`, `V(q)>=0`, `V(q0)=0` with `V(q)>0` below `q0`, absolutely continuous finite-energy `q`, and the pinned endpoints. The pointwise identity

```
A m^2 q_s^2/4 + V
  = (sqrt(A)|m| q_s/2 - sqrt(V))^2
    + sqrt(A)|m| q_s sqrt(V)
```

and a primitive of `sqrt(V)` give, *also for nonmonotone radial profiles*,

```
T >= 2 pi sqrt(A)|m| integral_0^q0 sqrt(V(q)) dq.             (1)
```

Equality requires `q_s=2 sqrt(V(q))/(sqrt(A)|m|)` almost everywhere. For both stated potentials this ODE starts with positive finite `q_s(0)` and approaches `q0` exponentially in `s` without reaching it at finite `s`; it furnishes a finite-energy equality profile in the continuous radial class. The `V_spec` square root has a nonanalytic core, so no claim of a smooth continuum matrix extremizer follows from this equality alone.

Substitution `q=b^2` into (1) yields

```
T_spec,inf = (sqrt(2) pi/3) |m| beta^3
           = (2 sqrt(2) pi/3) |k| beta^3                 (A=16),
T_tr,inf   = 2 pi sqrt(32w)|m| integral_0^q0
                 (q0-q) sqrt(K(q)) dq                    (A=32),
T_tr,inf   > pi sqrt(32 w K_lead)|m| b0^4               (b0,w>0).
```

At fixed `w>0`, `s0=b0=delta/2`, this exact trace integral obeys

```
T_tr,inf = [pi sqrt(2w)/2] |m| delta^4
           [1 + (9/8)delta^2 + O(delta^4)].
```

The cubic spectral power follows from an eigenvalue-distance potential linear in `b0-b` after the square root; the quartic trace power follows from a squared polynomial in `q0-q`. The absolute coefficients also differ because `A` differs. The R25 constant `K_lead` supplies the same leading power, but its Gaussian first-order profile does **not** attain the exact variable-`K` energy at finite `delta`.

## Origin regularity: infimum is not always an attained smooth minimum

In both models `V(0)>0`, so equality gives `q(s)=c s+o(s)`, `c=2 sqrt(V(0))/(sqrt(A)|m|)>0`, hence `b(rho)=sqrt(c/2) rho+o(rho)` for every winding. For `|m|>1`, `b(rho) cos(m phi)` and `b(rho) sin(m phi)` are not differentiable at the axis. In particular for `m=2`, the first diagonal displacement on `x>0` and `x<0` is proportional to `|x|`, whose one-sided derivatives disagree. A smooth matrix field with a nonzero angular `m` harmonic requires at least `b=O(rho^|m|)` at its core. Therefore **`m=2` has the same sharp infimum but does not attain it within smooth fields**; a finite lattice close to the bound cannot prove smooth attainment. The spectral `m=1` profile is C1 but generally not C2 at the origin because `sqrt(q)` creates a term of order `rho^2` in `b`; the trace `m=1` profile is smooth because its first-order equation is analytic in `q` near zero. Higher winding equality profiles fail already C1.

Sharpness for smooth winding `|m|>1` follows by an explicit recovery sequence. At `s_e=epsilon^2/2`, replace the equality profile for `s<=s_e` by `q_e(s)=q_BPS(s_e)(s/s_e)^|m|`. Then `b_e` has a polynomial angular core; `q_BPS(s_e)=c s_e+o(s_e)`, `q'_e=O(1)` and `V(q_e)=O(1)` there, so the energy discrepancy inside the disk is `O(s_e)=O(epsilon^2)`. Join derivatives across a thin annulus by a smooth cutoff with bounded derivative; it also costs `O(epsilon^2)`. For Report 018 `beta=.3,m=2`, an independent SciPy quadrature of the square-gap for disk radii `0.4,0.2,0.1,0.05` returned respectively `0.006609115673,0.001758505244,0.000454806699,0.0001157168372`; these finite values corroborate but are not needed for the exact recovery proof. The same local estimates apply to the trace potential.

## Conditional equal-tension partition corollary

For **isolated, straight, fixed-slot carriers in one declared action** with the same exterior and potential parameters, (1) gives `T_inf(k)=C|k|`: `C=2 sqrt(2) pi beta^3/3` for Report 018, while the R25 trace action has its own `C=4 pi sqrt(32w) integral_0^q0 (q0-q) sqrt(K(q)) dq`. If a unit charge's biaxial boundary forces total transverse index 2, the partitions `{2,2}`, `{2,1,1}` and `{1,1,1,1}` in half-index units each have `sum k=2`, hence the **sum of these separate strand infima** is `2C`. At `beta=0.3` the Report-018 sum is `0.1599437858` in its own units. The [post-harvest discussion](https://github.com/vantasnerdan/substrate-framework/discussions/186#discussioncomment-18604164) independently reports the `sqrt(2)` ratio of half-splitting radii; the pinned analytic `r_half` formula gives approximately `0.73912 sqrt(beta)` for `k=1/2` and `1.04528 sqrt(beta)` for `k=1`. The reported numerical coefficients `0.738` and `1.044` are close but not those analytic coefficients to the quoted precision.

This is **not** exact degeneracy of finite-separation three-dimensional particles: the `k=1` equality profile is not smooth at its axis, and the smooth sharp infimum requires the core recovery sequence above. Shared cores, inter-strand interactions, released slots, boundaries and potential choice are outside the additive calculation. A constrained `{2,2}` branch comparison and a separation/box ladder with charge, residual and mesh gates would discriminate core or interaction preference from boundary-driven splitting; neither the author's two finite-code split observations nor this corollary establishes a bound localized charge.

## Boundaries and next physical work

For the **R25 curvature-plus-trace-potential action itself**, the fixed-slot infimum is demonstrably **not** a lower bound on its larger static field class. Let the formerly fixed director eigenvalue `d=1` vary by a nonnegative smooth compactly supported bump in an annulus where `0<q<q0`. The new diagonal derivative is in the disjoint 33 block and commutes with every transverse-pair derivative, so `F_xy` and its curvature energy are unchanged, with **no director-gradient penalty in this action**. At `d=1`, direct differentiation of the four trace differences gives

```
dV_tr/dd = 2w(q-q0)[4+18s0+4(12s0^2+2(q+q0))] < 0
```

for `s0=b0>0`, `w>0`. Hence a small positive bump lowers the energy strictly. For `m=1` the trace-potential equality profile is smooth; for `m=2` a smooth core recovery sequence approaches its fixed-slot infimum while the same annular bump lowers energy by a definite amount. This counterexample rules out transferring (1) as a full-static-field bound in the specified action; adding a different kinetic term, constraint or action would require a new calculation. Allowing an unpinned axis `q(0)=q0` separately admits a zero-tension constant vacuum, so the regular-core bound then does not apply. No three-dimensional tension/energy comparison, loop radius, moving-frame hyperbolicity, source/response asymptotic, or two-body connected-minimizer result is supplied here. Those are actual next prerequisites before interpreting the radial result as an electron or confinement. #211 remains the broader R14–R25 referee frontier.

## Reproduction

At the target repository root, run `PYTHONPATH=src .venv/bin/python -m pytest tests/test_m5_strand.py -q`; `pair_curvature_xy` independently constructs the 2x2 commutator and the test constructs mixed 4x4 trace differences. The tests discriminate a halved/mismatched norm, changed winding, unpinned core, swapped potential and frozen-versus-variable stiffness. Run the author's two exact source scripts only at the pinned commits above; their passing status is corroboration rather than acceptance of the full physical conclusions.
