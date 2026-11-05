"""Independent entrywise SymPy antisymmetrization and all-operand differential."""
import itertools
import numpy as np
import sympy as sp
from second_born_functional import terms_from_vertices, evaluate, jvp, transpose


def run():
    rng = np.random.default_rng(4402)
    names = ('X','Y','U','V','greater1','greater2','lesser1','lesser2','reverse_greater','reverse_lesser')
    symbols = {name:sp.Matrix(2,2,lambda i,j:sp.Symbol(f'{name}_{i}{j}')) for name in names}
    flat = [v for matrix in symbols.values() for v in matrix]
    values = {s:sp.Rational(int(rng.integers(-5,6)),7)+sp.I*sp.Rational(int(rng.integers(-5,6)),11) for s in flat}
    directions = {s:sp.Rational(int(rng.integers(-4,5)),5)+sp.I*sp.Rational(int(rng.integers(-4,5)),13) for s in flat}
    x,y,u,v = (symbols[k] for k in names[:4])
    def component(f,h,r):
        out = sp.zeros(2)
        for i,j in itertools.product(range(2),repeat=2):
            out[i,j] = sp.expand(sum((x[i,b]*y[a,c]-x[i,c]*y[a,b])
                *(u[d,j]*v[e,z]-u[d,z]*v[e,j])*symbols[f][b,d]
                *symbols[h][c,e]*symbols[r][z,a]/2
                for a,b,c,d,e,z in itertools.product(range(2),repeat=6)))
        return out
    symbolic = dict(greater=component('greater1','greater2','reverse_lesser'),
                    lesser=component('lesser1','lesser2','reverse_greater'))
    symbolic['retarded'] = symbolic['greater']-symbolic['lesser']
    convert = lambda m, subs:np.array(m.subs(subs)).astype(complex)
    factors = [convert(symbols[k],values)[None] for k in names[:4]]
    df = [convert(symbols[k],directions)[None] for k in names[:4]]
    terms = terms_from_vertices(((factors[0],factors[1]),)*2,((factors[2],factors[3]),)*2)
    legs = {k:convert(symbols[k],values) for k in names[4:]}
    direction = {k:convert(symbols[k],directions) for k in names[4:]}
    actual = evaluate(terms,legs)
    tangent = jvp(terms,legs,direction,factor_direction=[df]*4)
    dual = {o:{'total':rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))} for o in symbolic}
    pullback = transpose(terms,legs,dual)
    lhs = 0j
    for order, matrix in symbolic.items():
        expected = convert(matrix,values)
        dmatrix = matrix.applyfunc(lambda e:sum(sp.diff(e,s)*directions[s] for s in flat))
        expected_jvp = convert(dmatrix,values)
        assert np.allclose(actual[order]['total'],expected,rtol=2e-11,atol=2e-11)
        assert np.allclose(tangent[order]['total'],expected_jvp,rtol=2e-11,atol=2e-11)
        lhs += np.sum(dual[order]['total']*expected_jvp)
    rhs = sum(np.sum(pullback['legs'][k]*direction[k]) for k in legs)
    rhs += sum(np.sum(g*d) for row in pullback['factors'] for g,d in zip(row,df,strict=True))
    conjugating = sum(np.vdot(pullback['legs'][k],direction[k]) for k in legs)
    assert abs(lhs-rhs)<2e-10*max(1,abs(lhs))
    assert abs(lhs-conjugating)>1e-7
    objective = sum(sp.sympify(complex(dual[order]['total'][i,j]))*matrix[i,j]
                    for order,matrix in symbolic.items() for i,j in itertools.product(range(2),repeat=2))
    for index,name in enumerate(names):
        expected_dual = np.array([[complex(sp.diff(objective,symbols[name][i,j]).subs(values))
                                  for j in range(2)] for i in range(2)])
        actual_dual = (sum(row[index][0] for row in pullback['factors'])
                       if index < 4 else pullback['legs'][name])
        assert np.allclose(actual_dual,expected_dual,rtol=2e-10,atol=2e-10)
    direct_only = actual['greater']['direct']
    assert np.linalg.norm(direct_only-actual['greater']['total'])>1e-7
    mutated = dict(legs,reverse_lesser=legs['reverse_lesser'].conj().T)
    assert np.linalg.norm(evaluate(terms,mutated)['greater']['total']-actual['greater']['total'])>1e-7
    return dict(schema='p254-a44-symbolic-oracle-v1',bilinear_dual_defect=float(abs(lhs-rhs)),
                mutations_exposed=['direct-only','conjugating-transpose','inferred-reversed-leg'],
                independent_layer='entrywise antisymmetric dense vertices, exact complex rational polynomial derivatives')
