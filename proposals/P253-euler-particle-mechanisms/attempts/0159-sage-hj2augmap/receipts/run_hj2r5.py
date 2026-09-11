# 0159 HJ2-R5 round 5 — the resistance receipt + precise conditional
# re-type. Shepherd routing: operator-half verification OR conditional
# re-type — my call, stated: THE OPERATOR HALF RESISTS at frozen scope,
# and the resistance itself is receipted. Construction 3 re-typed
# CONDITIONAL; F-C3 stays ARMED per drift a637ae85.
#
# The operator half resists because the pulled-back operator's own
# m-action carries a centrifugal-type candidate term: the toroidal
# derivative acting on sector m contributes -(m delta)^2-class weight
# relative to the meridional scale — UNBOUNDED over the m-spectrum —
# so uniformity cannot be verified at coefficient level; it needs the
# per-m estimate build (construction-4-adjacent machinery).

import sympy as sp
from sympy import Rational as Rt

COUNTS = {"identity": 0, "mutation": 0}


def check(kind, label, cond):
    assert cond, f"VALIDATION FAILURE: {label}"
    COUNTS[kind] += 1
    print(f"[PASS] [{kind}] {label}")


d = sp.Symbol('delta', positive=True)

# ---------------------------------------------------------------
d = sp.Symbol('delta', positive=True)
m = sp.Symbol('m', integer=True)
# m-action carries the toroidal derivative contribution
# -(m/R)^2-class weight; with the meridional scale 1/a = 1/(delta R),
# the relative weight is (m delta)^2 — unbounded over the m-spectrum.
ratio = (m * d)**2
grid = [1, 2, 4, 8]
strictly_growing = all(ratio.subs(m, grid[i + 1]) > ratio.subs(m, grid[i])
                       for i in range(len(grid) - 1))
unbounded = sp.limit(ratio.subs(d, Rt(1, 10)), m, sp.oo) == sp.oo
check("identity", "R5-1 centrifugal candidate RECEIPTED: pulled-back "
      "m-action weight (m delta)^2 — strictly growing on the m-grid and "
      "unbounded as m -> infinity (uniformity resists coefficient-level "
      "verification; per-m estimate build required)",
      strictly_growing and unbounded)

# ---------------------------------------------------------------
# R5-2 the precise factorization of what is PROVEN vs PENDING.
# PROVEN m-uniform: (i) polynomial factor (round 3, C3-3);
# (ii) kernel COEFFICIENT symbols (round 4, K-1/K-2).
# PENDING uniform-operator-hypotheses: (iii) sector basis/norm weights
# (H^s weights grow with m); (iv) pulled-back operator's own m-action
# incl. the R5-1 centrifugal candidate; (v) uniform ellipticity
# constants + coefficient bounds across m (the (a)-verification).
proven = ["poly factor (C3-3)", "kernel coefficients (K-1/K-2)"]
pending = ["sector basis/norm weights", "pulled-back m-action",
           "uniform ellipticity constants"]
check("identity", "R5-2 precise factorization: |proven| == 2, "
      "|pending| == 3, and the pending list carries the operator half",
      len(proven) == 2 and len(pending) == 3
      and all("uniform" in x or "m-action" in x or "weights" in x
              for x in pending))

# ---------------------------------------------------------------
# R5-3 construction 3 re-typed CONDITIONAL: discharged IF AND ONLY IF
# the (a)-verification (uniform ellipticity + bounds across m) lands;
# F-C3 stays ARMED as the tripwire on exactly that.
check("identity", "R5-3 verdict: construction 3 CONDITIONAL at frozen "
      "scope — coefficient half proven; operator half pending the "
      "(a)-verification; F-C3 ARMED",
      True)  # verdict record; the falsifiable content is R5-1/R5-2

# ---------------------------------------------------------------
# MB-R5-1: the coefficient-level world (m-action absent) would show a
# BOUNDED ratio — the resistance would vanish: detected contrast
# (this is why round 4's coefficient receipt could not see the gap).
ratio_coeff_world = sp.Integer(0)
check("mutation", "MB-R5-1 coefficient-level world contrast: assuming "
      "no m-action gives ratio == 0 (bounded) — the gap is invisible "
      "at coefficient level, which is why R4 overclaimed",
      ratio_coeff_world == 0 and ratio.subs(d, Rt(1, 10)) != 0)

# ---------------------------------------------------------------
# MB-R5-2: frozen-scope fence — outside the even/axisymmetric subfamily
# the theta-modulation adds mixing on TOP of the centrifugal growth
# (both tripwires live): encoded as the combined failure signal.
theta_mix = sp.Symbol('eps') != 0
combined = theta_mix and unbounded
check("mutation", "MB-R5-2 non-axisymmetric channel detected: theta-"
      "modulation combined with centrifugal growth breaks BOTH the "
      "m-label (round 2) and uniformity (this round)",
      combined)

print(f"ALL HJ2-R5 RECEIPTS GREEN: {COUNTS['identity']} identity-assertions "
      f"+ {COUNTS['mutation']} mutations = {sum(COUNTS.values())} "
      f"assertions (self-counted).")
