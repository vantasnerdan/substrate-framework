"""Finite exact algebra for the raw, non-Hilbert cubic time separation."""
import sympy as sp


def check():
    def matrix(name):
        return sp.Matrix(2, 2, sp.symbols(name + '0:4'))
    left, right, x, y, u, v = [matrix(n) for n in ('L', 'R', 'X', 'Y', 'U', 'V')]
    first, second = x * u.T, y * v.T
    direct = sum(left[a,b] * first[a,c] * second[b,d] * right[c,d]
                 for a in range(2) for b in range(2) for c in range(2) for d in range(2))
    separated = sum((x.T * left * y)[p,r] * (u.T * right * v)[p,r]
                    for p in range(2) for r in range(2))
    assert sp.expand(direct - separated) == 0
    # A W0 + W0 B^T is the simultaneous-time derivative at the origin.
    a, b, w = [matrix(n) for n in ('A', 'B', 'W')]
    t = sp.symbols('t')
    tangent = (sp.eye(2) + t*a) * w * (sp.eye(2) + t*b).T
    assert tangent.diff(t).subs(t,0) == a*w + w*b.T
    matrices = [left, right, x, y, u, v]
    substitutions = {}
    for j, m in enumerate(matrices):
        substitutions.update({symbol: sp.Integer(1+5*j+i) + sp.I*sp.Integer(2+i-j)
                              for i,symbol in enumerate(m)})
    exact = direct.subs(substitutions)
    mutations = {
        'right_Hilbert_instead_of_bilinear':
            sum((x.T*left*y)[p,r] * sp.conjugate((u.T*right*v)[p,r]) for p in range(2) for r in range(2)),
        'reverse_first_orientation': sum((x.T*left*y)[p,r]*(u*right*v)[p,r] for p in range(2) for r in range(2)),
        'missing_half': 2*separated,
        'missing_phase_mode': sum((x.T*left*y)[0,r]*(u.T*right*v)[0,r] for r in range(2)),
    }
    errors = {key: float(abs(complex(sp.expand(value.subs(substitutions) - exact))))
              for key,value in mutations.items()}
    assert min(errors.values()) > 0
    return {'status':'passed', 'exact_bilinear_identity':True,
            'stationarity_bridge':'(d_t+d_s) W(t,s) = exp(t A_q) [A_q W0_q + W0_q A_minus_q.T] exp(s A_minus_q.T); identically zero iff the raw defect is zero',
            'mutation_absolute_errors':errors,
            'normalization':'Pi=-i/2 sum_uv Luv(t) Ruv(s), canonical source normalization unchanged'}


if __name__ == '__main__':
    import json
    print(json.dumps(check(),sort_keys=True),flush=True)
