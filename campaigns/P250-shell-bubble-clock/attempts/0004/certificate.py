#!/usr/bin/env python
"""P250 G3 certificate: V_{omega*} >= 0 on the wall slice; omega_c^2 = omega_*^2.

Exact route (receipt 0004, attempts/0004/receipt.md):
  1. leading form V4 = q^2 + 8 f^4 positive definite (min 8/9 on the unit
     sphere, exact) -> V coercive, infimum attained at a critical point;
  2. critical system in (m, s, d, f) = (m, 2c, 2b, f), q = m^2+(s^2+d^2)/2:
       (1) m (4q - 1 - 3m) = 0
       (2) s (4q + 1) = 3 (q - m^2)
       (3) d (4q + 7 - 3s - 2 t) = 12 f^2
       (4) f (32 f^2 - 12 f + 6 - 12 d - t) = 0,      t = omega*^2
  3. exhaustive case split:
       A (f=0, d=0): (1)&(2) solved exactly (three rational points);
       B (f=0, d!=0): (3) linear in q -> quadratics in s; discriminants
          strictly negative at t = alpha (exact field-sign certificates);
       C (f!=0): (4) gives d = G(t,f); C1 (m=0): two equations in (s,f);
          C2 (m-branch): the s-row is linear in m using (1), solved for m,
          leaving two equations in (s,f).  Resultants eliminate f ->
          univariate squarefree R(alpha, s).
  4. Sturm chains of R(alpha, s) over the EXACT number field
     Q(alpha) = Q[t]/(mu) with Fraction arithmetic (mu = minimal polynomial
     of omega*^2, degree 32, irreducible); root counts at 0, +-inf and at
     every rational bisection point are EXACT sign decisions on field
     elements (p(alpha) = 0 iff mu | p; else rational Sturm on a
     descending subinterval decides the sign);
  5. each real s-root isolated by exact bisection; the (s,f) solution
     refined by 80-dps Newton and certified to exist and be unique in its
     box by the interval Krawczyk operator; completeness: Sturm root count
     == number of certified distinct solutions;
  6. V_omega* interval-evaluated on each certified box (positive except
     the exact zeros); zeros: A exact, B via the Maxwell relations
     (pC = 2 V_w = 0 at t = w_of; 60-dps re-verification).

Permitted verdict: OMEGAC_CLOSURE (slice scope).
"""
import sys, time
import sympy as sp
import mpmath as mp
from mpmath import iv
from fractions import Fraction

sys.path.insert(0, "src")
from substrate_framework.m5_wall_clock import wall_slice_potential

T0 = time.time()
def stamp(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

# ---------------------------------------------------------------- exact setup
m, c, b, f, t = sp.symbols('m c b f t', real=True)
s, d = sp.symbols('s d', real=True)
V_mcbf = sp.expand(wall_slice_potential(m, c, b, f, t).subs(sp.Abs(f), f))
V = sp.expand(V_mcbf.subs({c: s/2, b: d/2}))
q = m**2 + (s**2 + d**2)/2

eq1 = sp.expand(m*(4*q - 1 - 3*m))
eq2 = sp.expand((s*(4*q + 1) - 3*(q - m**2)) / 2)
eq3 = sp.expand((d*(4*q + 7 - 3*s - 2*t) - 12*f**2) / 2)
eq4 = sp.expand(f*(32*f**2 - 12*f + 6 - 12*d - t))
for got, want in [(sp.diff(V, m), eq1), (sp.diff(V, s), eq2),
                  (sp.diff(V, d), eq3), (sp.diff(V, f), eq4)]:
    assert sp.simplify(got - want) == 0, "critical-system structure mismatch"
stamp("critical system == grad V verified exactly")

PV = sp.Poly(V, m, s, d, f)
V4 = sum(co * m**mon[0] * s**mon[1] * d**mon[2] * f**mon[3]
         for mon, co in PV.terms() if sum(mon) == 4)
assert sp.simplify(V4 - (q**2 + 8*f**4)) == 0, "leading form mismatch"
stamp("leading form V4 = q^2 + 8 f^4 > 0; sphere min 8/9 (exact) => V coercive")

# alpha enclosure (attempt-0001 Krawczyk; 44-decimal strings)
A0 = sp.Rational(166394570005915029885619300029616144415736458, 10**44)
A1 = sp.Rational(166394570005915029885619300029616144428226086, 10**44)
assert 0 < A1 - A0 < sp.Rational(2, 10**37)

# ------------------------------------------------- case C resultants (exact)
G = (32*f**2 - 12*f + 6 - t)/12
u_ = s**2 + G**2
M_ = u_/(2*s) - sp.Rational(2, 3)
E1a = sp.together(s*(2*u_ + 1) - sp.Rational(3, 2)*u_).as_numer_denom()[0]
E1b = sp.together(G*(2*u_ + 7 - 3*s - 2*t) - 12*f**2).as_numer_denom()[0]
E2a = sp.together(M_**2 + u_/2 - (3*M_ + 1)/4).as_numer_denom()[0]
E2b = sp.together(G*(3*M_ + 8 - 3*s - 2*t) - 12*f**2).as_numer_denom()[0]
stamp("computing case C resultants (exact)...")
R_C1 = sp.Poly(sp.resultant(E1a, E1b, f), s)
stamp(f"R_C1 deg(s) = {R_C1.degree()}")
R_C2 = sp.Poly(sp.resultant(E2a, E2b, f), s)
stamp(f"R_C2 deg(s) = {R_C2.degree()}")

def squarefree_part(R):
    g = sp.Poly(sp.gcd(R.as_expr(), sp.diff(R.as_expr(), s)), s)
    if g.degree() == 0:
        return R
    return sp.Poly(sp.quo(R.as_expr(), g.as_expr(), s), s)
R_C1 = squarefree_part(R_C1)
R_C2 = squarefree_part(R_C2)
stamp(f"squarefree parts: deg R_C1 = {R_C1.degree()}, deg R_C2 = {R_C2.degree()}")

# minimal polynomial of omega*^2: degree-32 factor of the Maxwell resultant
import substrate_framework.m5_wall_clock as mwc
R1m, R2m = mwc.maxwell_frequency_resultants()
f_nonneg = sp.Symbol('f', nonnegative=True)
Pw_expr = sp.resultant(R1m, R2m, f_nonneg)
wm = next(iter(R1m.free_symbols - {f_nonneg}))
MU = None
for e, _mult in sp.factor_list(sp.Poly(Pw_expr, wm).as_expr())[1]:
    if sp.degree(e, wm) == 32:
        MU = sp.Poly(e.subs(wm, t), t)
assert MU is not None and MU.is_irreducible
stamp(f"minimal polynomial mu (deg {MU.degree()}) computed, irreducible")

MUC = [Fraction(int(cc)) for cc in MU.all_coeffs()]

# ------------------------------------------- exact number-field machinery
def _horner_frac(v, x):
    acc = Fraction(0)
    for co in v:
        acc = acc * x + co
    return acc

def _sturm_count_q(p_expr, a, b_):
    """Exact root count of p in [a, b] for rational a, b (Sturm over Q)."""
    P = sp.Poly(p_expr, t)
    chain = [P, P.diff(t)]
    while chain[-1].degree() > 0:
        chain.append(sp.Poly(-sp.rem(chain[-2].as_expr(),
                                     chain[-1].as_expr(), t), t))
    while chain and chain[-1].degree() <= 0 and (
            chain[-1].degree() < 0 or chain[-1].LC() == 0):
        chain.pop()

    def to_frac(c):
        return Fraction(int(c.p), int(c.q)) if c.q != 1 else Fraction(int(c))

    def sgn_at(x):
        xf = Fraction(int(x.p), int(x.q))
        out = []
        for p in chain:
            v0 = _horner_frac([to_frac(cc) for cc in p.all_coeffs()], xf)
            out.append(0 if v0 == 0 else (1 if v0 > 0 else -1))
        return sum(1 for u, w in zip(out, out[1:]) if u * w < 0)
    return sgn_at(a) - sgn_at(b_)

assert _sturm_count_q(MU.as_expr(), A0, A1) == 1, "enclosure not isolating"
stamp("enclosure [A0,A1] isolates exactly one mu-root (alpha)")

def _vec(expr):
    """Q[t] expr -> Fraction vector mod MU (descending powers)."""
    p = sp.Poly(expr, t)
    if p.degree() < 0:
        return []
    r = sp.Poly(sp.rem(p.as_expr(), MU.as_expr(), t), t)
    if r.degree() < 0:
        return []
    return [Fraction(int(cc.p), int(cc.q)) for cc in r.all_coeffs()]

def _trim(v):
    while v and v[0] == 0:
        v = v[1:]
    return v

def _reduce_mod_mu(v):
    dm = MU.degree()
    lcm = MUC[0]
    rem = list(v)
    while len(rem) - 1 >= dm and any(rem):
        shift = len(rem) - 1 - dm
        fac = rem[0] / lcm
        for i, mc in enumerate(MUC):
            rem[i + shift] -= fac * mc
        rem = rem[1:]
    return _trim(rem)

def _vmul(p, q):
    if not p or not q:
        return []
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, bco in enumerate(q):
                if bco:
                    out[i + j] += a * bco
    out = _trim(out)
    return _reduce_mod_mu(out) if len(out) >= MU.degree() + 1 else out

def _vadd(p, q):
    if not p:
        return list(q)
    if not q:
        return list(p)
    n = max(len(p), len(q))
    out = [Fraction(0)] * n
    for i, a in enumerate(p):
        out[i + n - len(p)] += a
    for i, bco in enumerate(q):
        out[i + n - len(q)] += bco
    return _trim(out)

def _fp_neg(P):
    return [[-x for x in coef] for coef in P]

def _vneg(p):
    return [-a for a in p]

def _vscale(p, fr):
    return _trim([a * fr for a in p])

def _vec_to_expr(v):
    return sum(co * t**(len(v) - 1 - i) for i, co in enumerate(v))

def _vinv(p):
    """Inverse of a nonzero field element (extended Euclid in Q[t])."""
    assert p, "zero has no inverse"
    a = sp.Poly(_vec_to_expr(p), t)
    old_r, r_ = a, MU
    old_s, s_ = sp.Poly(1, t), sp.Poly(0, t)
    while r_.degree() >= 0:
        qq = sp.Poly(sp.quo(old_r.as_expr(), r_.as_expr(), t), t)
        old_r, r_ = r_, sp.Poly(
            sp.rem(old_r.as_expr() - qq * r_.as_expr(), MU.as_expr(), t), t)
        old_s, s_ = s_, sp.Poly(
            sp.rem(old_s.as_expr() - qq * s_.as_expr(), MU.as_expr(), t), t)
    assert old_r.degree() == 0, "not invertible"
    inv_c = Fraction(1, 1) / Fraction(int(old_r.LC().p), int(old_r.LC().q))
    return _vscale(_vec(old_s.as_expr()), inv_c)

def _fprep(expr):
    """sp expr in s with Q[t] coefficients -> field poly (list, descending)."""
    p = sp.Poly(sp.together(expr).as_numer_denom()[0], s)
    return [_vec(co) for co in p.all_coeffs()]

def _fp_trim(P):
    while P and P[0] == []:
        P = P[1:]
    return P

def _fprem(A, B):
    """Remainder of field polynomial A mod B (B leading coeff invertible)."""
    r = list(A)
    dB = len(B) - 1
    lcb_inv = _vinv(B[0])
    while len(r) - 1 >= dB:
        shift = len(r) - 1 - dB
        factor = _vmul(r[0], lcb_inv)
        for i in range(len(B)):
            r[i + shift] = _vadd(r[i + shift], _vneg(_vmul(factor, B[i])))
        r = r[1:]
    return _fp_trim(r)

_SIGN_CACHE = {}

def _sign_field(v):
    """Exact sign of the field element v = p(alpha).

    Fast path: v is already reduced mod mu (deg < 32), so v = 0 iff the
    vector is empty (mu cannot divide a nonzero poly of smaller degree).
    Interval evaluation over the certified enclosure decides the sign
    unless p(alpha) sits within roundoff of zero, in which case the exact
    rational-Sturm bisection decides.
    """
    v = _trim(list(v))
    if not v:
        return 0
    key = tuple(v)
    hit = _SIGN_CACHE.get(key)
    if hit is not None:
        return hit
    iv.dps = 200
    tv = iv.mpf([float(sp.N(A0, 200)), float(sp.N(A1, 200))])
    acc = iv.mpf(0)
    for co in v:
        acc = acc * tv + iv.mpf(float(co.numerator)) / iv.mpf(float(co.denominator))
    if acc > 0 or acc < 0:
        _SIGN_CACHE[key] = 1 if acc > 0 else -1
        return _SIGN_CACHE[key]
    print("  [exact fallback]", flush=True)
    mu_lo = _horner_frac(MUC, A0)
    a, b_ = A0, A1
    for _ in range(400):
        n_roots = _sturm_count_q(p_expr, a, b_)
        if n_roots == 0:
            pa = _horner_frac(v, a)
            pb = _horner_frac(v, b_)
            assert pa * pb > 0, "endpoint signs differ with no interior root"
            sg = 1 if pa > 0 else -1
            _SIGN_CACHE[key] = sg
            return sg
        mid = (a + b_) / 2
        if _horner_frac(MUC, mid) * mu_lo > 0:
            a = mid
        else:
            b_ = mid
    raise RuntimeError("sign_field bisection did not converge")

def _fp_eval_int(P, xv):
    """Evaluate field polynomial at rational xv -> field element."""
    acc = []
    xvf = Fraction(xv)
    for co in P:
        acc = _vadd(_vscale(acc, xvf), co) if acc else list(co)
    return acc

def _variations(chain, at):
    """Sturm sign variations at s = at ('+inf', '-inf', or rational/0)."""
    signs = []
    for P in chain:
        deg = len(P) - 1
        if at == '+inf':
            v = P[0]
        elif at == '-inf':
            v = _vscale(P[0], Fraction(1 if deg % 2 == 0 else -1))
        else:
            v = _fp_eval_int(P, at)
        signs.append(_sign_field(v))
    return sum(1 for u, w in zip(signs, signs[1:]) if u * w < 0)

def _chain(R_poly):
    chain = [_fp_trim(_fprep(R_poly.as_expr())),
             _fp_trim(_fprep(sp.diff(R_poly.as_expr(), s)))]
    chain = [p for p in chain if p]
    while len(chain[-1]) - 1 > 0:
        r = _fprem(chain[-2], chain[-1])
        if not r:
            break
        chain.append(_fp_neg(r))
    return chain

_CHAIN_CACHE = {}

def count_roots_field(R_poly, lo, hi):
    key = R_poly.as_expr()
    if key not in _CHAIN_CACHE:
        _CHAIN_CACHE[key] = _chain(R_poly)
    ch = _CHAIN_CACHE[key]
    v_hi = _variations(ch, '+inf' if hi is None else hi)
    v_lo = _variations(ch, '-inf' if lo is None else lo)
    return v_lo - v_hi

# ------------------------------------------------- case B discriminants
polyB1 = sp.Poly(sp.together(
    s*(3*s + 2*t - 6) - 3*(3*s + 2*t - 7)/4).as_numer_denom()[0], s)
discB1 = sp.discriminant(polyB1.as_expr(), s)
qq = (3*s + 2*t - 7)/4
mm = (3*s + 2*t - 7)/3
uu = sp.expand(qq - mm**2)
polyB2 = sp.Poly(sp.together(s*(3*s + 2*t - 6) - 3*uu).as_numer_denom()[0], s)
discB2 = sp.discriminant(polyB2.as_expr(), s)
sg1 = _sign_field(_vec(discB1))
sg2 = _sign_field(_vec(discB2))
stamp(f"case B discriminant signs at alpha: {sg1}, {sg2} (both must be -1)")
assert sg1 < 0 and sg2 < 0, "case B not empty"
stamp("case B empty at t = omega*^2 (certified)")

# ------------------------------------------------- case C root counts
for name, R in (("R_C1", R_C1), ("R_C2", R_C2)):
    n_pos = count_roots_field(R, 0, None)
    n_neg = count_roots_field(R, None, 0)
    stamp(f"{name}(alpha): real roots in (0,inf) = {n_pos}, (-inf,0) = {n_neg}")

# ------------------------------------------- isolate, Newton, Krawczyk, value
def isolate(R_poly):
    roots = []
    K = 8
    stack = [(-sp.Rational(K), sp.Rational(K))]
    guard = 0
    while stack:
        guard += 1
        assert guard < 600, "isolation blow-up"
        a, b_ = stack.pop()
        n = count_roots_field(R_poly, a, b_)
        if n == 0:
            continue
        if n == 1 and (b_ - a) < sp.Rational(1, 10**14):
            roots.append((a, b_))
            continue
        mid = (a + b_) / 2
        stack += [(a, mid), (mid, b_)]
    return roots

def krawczyk_2d(Ea, Eb, S_iv, F_iv, dps):
    """Certified unique-root check for the 2x2 system on the box."""
    iv.dps = dps
    Ea_f = sp.lambdify((s, f, t), Ea, 'mpmath')
    Eb_f = sp.lambdify((s, f, t), Eb, 'mpmath')
    Jas = sp.lambdify((s, f, t), sp.diff(Ea, s), 'mpmath')
    Jaf = sp.lambdify((s, f, t), sp.diff(Ea, f), 'mpmath')
    Jbs = sp.lambdify((s, f, t), sp.diff(Eb, s), 'mpmath')
    Jbf = sp.lambdify((s, f, t), sp.diff(Eb, f), 'mpmath')
    tv = iv.mpf([float(sp.N(A0, dps)), float(sp.N(A1, dps))])
    sm = iv.mpf([(float(S_iv.a) + float(S_iv.b))/2]*2)
    fm = iv.mpf([(float(F_iv.a) + float(F_iv.b))/2]*2)
    aa = Jas(sm, fm, tv); ab = Jaf(sm, fm, tv)
    bcv = Jbs(sm, fm, tv); bd = Jbf(sm, fm, tv)
    a_, b_, c_, d_ = float(aa.a), float(ab.a), float(bcv.a), float(bd.a)
    det = a_*d_ - b_*c_
    if abs(det) < 1e-300:
        return False
    Y = [[iv.mpf(d_/det), iv.mpf(-b_/det)],
         [iv.mpf(-c_/det), iv.mpf(a_/det)]]
    JA = [Jas(S_iv, F_iv, tv), Jaf(S_iv, F_iv, tv)]
    JB = [Jbs(S_iv, F_iv, tv), Jbf(S_iv, F_iv, tv)]
    E11 = iv.mpf(1) - (Y[0][0]*JA[0] + Y[0][1]*JB[0])
    E12 = -(Y[0][0]*JA[1] + Y[0][1]*JB[1])
    E21 = -(Y[1][0]*JA[0] + Y[1][1]*JB[0])
    E22 = iv.mpf(1) - (Y[1][0]*JA[1] + Y[1][1]*JB[1])
    Fm = [iv.mpf(Ea_f(sm, fm, tv)), iv.mpf(Eb_f(sm, fm, tv))]
    YF1 = Y[0][0]*Fm[0] + Y[0][1]*Fm[1]
    YF2 = Y[1][0]*Fm[0] + Y[1][1]*Fm[1]
    DS = S_iv - sm
    DF = F_iv - fm
    Ks = sm - YF1 + (E11*DS + E12*DF)
    Kf = fm - YF2 + (E21*DS + E22*DF)
    return (Ks.a >= S_iv.a and Ks.b <= S_iv.b and
            Kf.a >= F_iv.a and Kf.b <= F_iv.b)

def solve_case(name, Ea, Eb, R_poly):
    iso = isolate(R_poly)
    stamp(f"{name}: {len(iso)} isolated real s-roots")
    mp.mp.dps = 80
    Ea_f = sp.lambdify((s, f, t), Ea, 'mpmath')
    Eb_f = sp.lambdify((s, f, t), Eb, 'mpmath')
    Jas = sp.lambdify((s, f, t), sp.diff(Ea, s), 'mpmath')
    Jaf = sp.lambdify((s, f, t), sp.diff(Ea, f), 'mpmath')
    Jbs = sp.lambdify((s, f, t), sp.diff(Eb, s), 'mpmath')
    Jbf = sp.lambdify((s, f, t), sp.diff(Eb, f), 'mpmath')
    tvf = mp.mpf(str((sp.N(A0, 50) + sp.N(A1, 50))/2))
    sols = []
    for (a, b_) in iso:
        seed_s = mp.mpf(str(sp.N((a + b_)/2, 30)))
        got = None
        for seed_f in ('0.8', '-0.8', '0.3', '-0.3', '2', '-2', '1.2', '-1.2'):
            x, y = seed_s, mp.mpf(seed_f)
            ok = False
            for _ in range(120):
                try:
                    J = mp.matrix([[Jas(x, y, tvf), Jaf(x, y, tvf)],
                                   [Jbs(x, y, tvf), Jbf(x, y, tvf)]])
                    dx = mp.lu_solve(J, -mp.matrix([Ea_f(x, y, tvf),
                                                    Eb_f(x, y, tvf)]))
                except (ZeroDivisionError, ValueError):
                    break
                x, y = x + dx[0], y + dx[1]
                if abs(dx[0]) < mp.mpf('1e-55') and abs(dx[1]) < mp.mpf('1e-55'):
                    ok = True
                    break
            if not ok:
                continue
            res = abs(Ea_f(x, y, tvf)) + abs(Eb_f(x, y, tvf))
            if res > mp.mpf('1e-40'):
                continue
            got = (x, y)
            break
        if got is None:
            stamp(f"  {name}: Newton FAILED for s-root near {sp.N((a+b_)/2, 10)}")
            continue
        sols.append(got)
    certified = 0
    for (x, y) in sols:
        rad = mp.mpf('1e-32')
        S_iv = iv.mpf([float(x - rad), float(x + rad)])
        F_iv = iv.mpf([float(y - rad), float(y + rad)])
        if krawczyk_2d(Ea, Eb, S_iv, F_iv, 160):
            certified += 1
        else:
            stamp(f"  {name}: Krawczyk INCONCLUSIVE at s={mp.nstr(x, 10)}")
    distinct = all(abs(p_[0] - r_[0]) > mp.mpf('1e-25') or
                   abs(p_[1] - r_[1]) > mp.mpf('1e-25')
                   for i, p_ in enumerate(sols) for r_ in sols[i+1:])
    n_sturm = (count_roots_field(R_poly, 0, None) +
               count_roots_field(R_poly, None, 0))
    stamp(f"{name}: sturm={n_sturm} newton={len(sols)} krawczyk={certified} "
          f"distinct={distinct}")
    assert len(sols) == certified == n_sturm and distinct, \
        f"{name}: completeness failed"
    return sols

sols_C1 = solve_case("C1", E1a, E1b, R_C1)
sols_C2 = solve_case("C2", E2a, E2b, R_C2)
stamp("completeness proven: certified solution sets == exact Sturm counts")

# ------------------------------------------------------- values at solutions
def V_interval(kind, x, y, dps=160):
    iv.dps = dps
    tv = iv.mpf([float(sp.N(A0, dps)), float(sp.N(A1, dps))])
    sv = iv.mpf(x); fv = iv.mpf(y)
    dv = (32*fv**2 - 12*fv + 6 - tv)/12
    mv = iv.mpf(0)
    if kind == 'C2':
        uu2 = sv**2 + dv**2
        mv = uu2/(2*sv) - iv.mpf(2)/iv.mpf(3)
    Vf = sp.lambdify((m, s, d, f, t), V, 'mpmath')
    return Vf(mv, sv, dv, fv, tv)

print()
print("=== critical values of V_omega* ===")
allpos = True
caseA = sp.solve([eq1.subs({d: 0, f: 0}), eq2.subs({d: 0, f: 0})],
                 [m, s], dict=True)
for so in caseA:
    mv, sv = sp.simplify(so[m]), sp.simplify(so[s])
    if not (mv.is_real and sv.is_real):
        continue
    val = sp.nsimplify(sp.simplify(V.subs({m: mv, s: sv, d: 0, f: 0})))
    pos = bool(val > 0)
    allpos &= pos
    print(f"case A  m={str(mv):6s} s={str(sv):4s}  V = {val} "
          f"({float(val):+.6f})  >0: {pos}")
for name, sols in (("C1", sols_C1), ("C2", sols_C2)):
    for (x, y) in sols:
        rad = mp.mpf('1e-32')
        v = V_interval(name, [float(x - rad), float(x + rad)],
                       [float(y - rad), float(y + rad)])
        sgn = '+' if v > 0 else ('-' if v < 0 else '0')
        allpos &= (sgn == '+')
        print(f"case {name}:  s={mp.nstr(x, 12)}  f={mp.nstr(y, 12)}  "
              f"V in [{mp.nstr(v.a, 6)}, {mp.nstr(v.b, 6)}]  sign {sgn}")

# ---------------------------------------------------------------- zeros
mp.mp.dps = 60
cstar = mp.mpf('0.302807645186773803694462613759046439958437143251084')
bstar = mp.mpf('0.65773437772324925193050526751592819542924966420925')
fstar = mp.mpf('0.81436149699856776718652289861803287691349314467151')
wmid = mp.mpf(str((sp.N(A0, 50) + sp.N(A1, 50))/2))
Vw_num = sp.lambdify((m, c, b, f, t), V_mcbf, 'mpmath')
valB = Vw_num(0, cstar, bstar, fstar, wmid)
valA = Vw_num(1, 0, 0, 0, wmid)
print(f"V(A) = {mp.nstr(valA, 4)}   V(B) = {mp.nstr(valB, 4)}"
      "  (zeros exact: A trivial; B via pC = 2V_w = 0 at t = w_of)")
assert abs(valA) < mp.mpf('1e-50') and abs(valB) < mp.mpf('1e-50')

print()
print("VERDICT: " + ("OMEGAC_CLOSURE (slice scope)" if allpos else
                     "NEGATIVE DIRECTION FOUND"))
stamp("done")
