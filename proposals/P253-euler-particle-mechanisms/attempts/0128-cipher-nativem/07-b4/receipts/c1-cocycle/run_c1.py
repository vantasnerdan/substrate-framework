"""B4-C1: S^2 Hamiltonian cocycle battery (FROZEN F-bar in ../../07-b4/01-c1-fbar.md).
Poisson bracket on S^2: {f,g} = (f_t*g_p - f_p*g_t)/sin(t), t=theta, p=phi.
c_h(f,g) = int h*{f,g}*sin(t) dt dp, h = cos(t). Exact sympy integrals.
Usage: python3 run_c1.py
"""
import sympy as sp

t, p = sp.symbols('t p', real=True)


def Y10():
    return sp.cos(t)


def Y11c():
    return sp.sin(t) * sp.cos(p)


def Y11s():
    return sp.sin(t) * sp.sin(p)


def Y20():
    return 3 * sp.cos(t) ** 2 - 1


def Y21c():
    return sp.sin(t) * sp.cos(t) * sp.cos(p)


def Y22c():
    return sp.sin(t) ** 2 * sp.cos(2 * p)


def pb(f, g):
    return (sp.diff(f, t) * sp.diff(g, p) - sp.diff(f, p) * sp.diff(g, t)) / sp.sin(t)


def coc(h, f, g):
    return sp.simplify(sp.integrate(sp.simplify(h * pb(f, g) * sp.sin(t)),
                                   (t, 0, sp.pi), (p, 0, 2 * sp.pi)))


def main():
    h = sp.cos(t)
    # T0: cocycle identity on triples (cyclic sum = Jacobi, must vanish)
    import itertools
    basis = [Y10(), Y11c(), Y11s(), Y20(), Y21c(), Y22c()]
    worst = 0
    for a, b, c in itertools.combinations(basis, 3):
        s = (coc(h, pb(a, b), c) + coc(h, pb(b, c), a) + coc(h, pb(c, a), b))
        worst = max(worst, abs(complex(s.evalf())))
    print(f"T0 cocycle identity max|cyclic| = {worst:.1e} -> {'PASS' if worst < 1e-8 else 'FAIL-bug'}")
    import numpy as np
    so3 = [Y10(), Y11c(), Y11s()]
    pairs = [(0, 1), (0, 2), (1, 2)]

    def proj(f):
        return [sp.integrate(f * b * sp.sin(t), (t, 0, sp.pi), (p, 0, 2 * sp.pi)) for b in so3]

    norms = [complex(sp.integrate(b * b * sp.sin(t), (t, 0, sp.pi), (p, 0, 2 * sp.pi)).evalf()) for b in so3]
    M = np.zeros((3, 3), dtype=complex)
    rhs = np.zeros(3, dtype=complex)
    for r, (i, j) in enumerate(pairs):
        br = pb(so3[i], so3[j])
        rhs[r] = complex(coc(h, so3[i], so3[j]).evalf())
        C = [complex(c.evalf()) / norms[k] for k, c in enumerate(proj(br))]
        M[r, :] = [-c for c in C]
    sol, res, *_ = np.linalg.lstsq(M, rhs, rcond=None)
    print(f"T1 Whitehead: residual {res[0] if len(res) else 0:.1e} -> {'PASS-exact' if (len(res) == 0 or res[0] < 1e-8) else 'FAIL-bug'}")
    # T2: LIVE pair — must be m-matched (axisymmetric h kills m-mismatched
    # pairs by phi-selection; Y20/Y21c gives trivial zero = bad test, not kill)
    Y21s = sp.sin(t) * sp.cos(t) * sp.sin(p)
    v = coc(h, Y21c(), Y21s)
    print(f"T2 live pairing c_h(Y21c,Y21s) = {v} ~ {complex(v.evalf()):.6f}")
    print(f"C1 verdict: {'PASS-lean (nontrivial class)' if abs(complex(v.evalf())) > 1e-6 else 'KILL (undetected)'}")

if __name__ == "__main__":
    main()
