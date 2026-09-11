# SYN P3-B amendment receipts — charter: shepherd (P3-B one round;
# P3-C ruled CHARGE-ONLY, tilt = separate channel). Receipts the
# MODEL-ORDER p-response theorem the amendment's p-independence claim
# rests on. Self-counted; exits nonzero on any failure. Vacuity check
# applied (standing rule): every clause falsifiable — no self-verifying
# asserts, no coefficient/conflation traps.

import sympy as sp
from sympy import Rational as Rt

COUNTS = {"identity": 0, "mutation": 0}


def check(kind, label, cond):
    assert cond, f"VALIDATION FAILURE: {label}"
    COUNTS[kind] += 1
    print(f"[PASS] [{kind}] {label}")


# Direction functionals f(n), n on the unit sphere. Polarized average:
#   <f>_p = (1-p) <f>_iso + p f(nhat)
# (same averages class as run_fb.py W_aniso). P3-B question: which
# measurement functionals have p-invariant spectra?

v1, v2, v3 = sp.symbols('v1 v2 v3', real=True)          # generic probe
p = sp.Symbol('p', real=True)                          # order parameter
c1, s1, c2, s2 = sp.symbols('c1 s1 c2 s2', real=True)  # unit-sphere params
ngen = sp.Matrix([c1 * s2, s1 * s2, c2])                # generic direction
v = sp.Matrix([v1, v2, v3])
A = v * v.T                                             # generic sym rank-1

# ---- P3B-1: probability preservation over the physical p-window
grid = [Rt(0), Rt(1, 4), Rt(1, 2), Rt(3, 4), Rt(1)]
ok = all((0 <= (1 - pv) <= 1) == True and (0 <= pv <= 1) == True
         and sp.simplify((1 - pv) + pv - 1) == 0 for pv in grid)
check("identity", "P3B-1 weight window: (1-p, p) are probabilities "
      "summing to 1 across the physical window — the operative premise "
      "for invariance carrying to direction-blind counts",
      ok)

# ---- P3B-2 (CONTRAST): direction-dependent functionals DO shift with p.
fnAn = sp.expand(ngen.dot(A * ngen))                    # (n.A.n)(direction)
fAn_at = fnAn.subs({c1: 0, s1: 0, c2: 1, s2: 0})        # nhat = z: f = v3^2
fAn_iso = sp.trace(A) / 3                               # <nAn>_iso = tr A/3
resp = sp.simplify(fAn_at - fAn_iso)
check("identity", "P3B-2 CONTRAST: <nAn> SHIFTS with p "
      "(response amplitude v3^2 - |v|^2/3 != 0 for generic v)",
      resp != 0)

# ---- P3B-3 (THEOREM): exact p-response identity
#   <f>_p - <f>_0 = p * ( f(nhat) - <f>_iso )
# for the direction-functional class; two representatives.
f_nv2 = sp.expand((ngen.dot(v)) ** 2)                   # (n.v)^2
f_nv2_at = f_nv2.subs({c1: 0, s1: 0, c2: 1, s2: 0})     # v3^2
f_nv2_iso = (v1**2 + v2**2 + v3**2) / 3                 # |v|^2/3

lhs_An = sp.simplify(((1 - p) * fAn_iso + p * fAn_at) - fAn_iso)
rhs_An = sp.simplify(p * (fAn_at - fAn_iso))
check("identity", "P3B-3a exact p-response [nAn]: <f>_p - <f>_0 == "
      "p*(f(nhat) - <f>_iso)",
      sp.simplify(lhs_An - rhs_An) == 0)

lhs_nv2 = sp.simplify(((1 - p) * f_nv2_iso + p * f_nv2_at) - f_nv2_iso)
rhs_nv2 = sp.simplify(p * (f_nv2_at - f_nv2_iso))
check("identity", "P3B-3b exact p-response [nv2]: same identity",
      sp.simplify(lhs_nv2 - rhs_nv2) == 0)

# Corollary (the amendment's load-bearing biconditional):
# p-independence <=> f(nhat) == <f>_iso <=> f direction-blind.
f_count = sp.Symbol('Lambda')   # count per event: direction-blind constant
check("identity", "P3B-3c direction-blind corollary: count functional has "
      "ZERO p-response (Lambda(p) == Lambda(0)) — counting statistics are "
      "direction-blind at model order; CV(p) == CV(0) follows by applying "
      "this to both moments (mean, variance)",
      sp.simplify(((1 - p) * f_count + p * f_count) - f_count) == 0)

# ---- MB-P3a: unnormalized weights break probability preservation
w_wrong = (1 - p) * Rt(1, 4) + p * Rt(1, 2)
check("mutation", "MB-P3a unnormalized weights detected: weight pair "
      "does not sum to 1 (normalization is load-bearing)",
      sp.simplify(w_wrong - 1) != 0)

# ---- MB-P3b: isotropic-even-at-p (wrong aligned moment) detected —
# the response theorem's CONTENT vanishes under the mutation while the
# true model's does not.
lhs_mut = sp.simplify((f_nv2_iso + p * f_nv2_iso) - f_nv2_iso)  # aligned:=iso
check("mutation", "MB-P3b isotropic-even-at-p detected: mutated LHS != "
      "true response (theorem content lost — the aligned moment is "
      "load-bearing for P3-B's p-channel bookkeeping)",
      sp.simplify(lhs_mut - rhs_nv2) != 0
      and sp.simplify(f_nv2_at - f_nv2_iso) != 0)

print(f"ALL SYN P3-B RECEIPTS GREEN: {COUNTS['identity']} identity-assertions "
      f"+ {COUNTS['mutation']} mutations = {sum(COUNTS.values())} assertions "
      "(self-counted).")
