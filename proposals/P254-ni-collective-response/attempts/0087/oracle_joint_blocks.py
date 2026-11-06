"""Narrow new IF/FI/FF boundary oracle. Explicit execution only; no old suite.

Run: python oracle_joint_blocks.py --output /fresh/path/oracle.json
Requires SymPy; exercises independent non-Hermitian vertex and contour legs.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
import sympy as sp
from joint_blocks import sigma_correction, sigma_rectangular
from electron_ion_functional import EdgePair, evaluate_components
from history_banks import binding, save


def exact(left, electron, right, boson, rows, cols):
    m = left.shape[1]
    return sp.Matrix(m, m, lambda a, d: sp.I*sum(
        sp.sympify(left[i,a,b])*sp.sympify(electron[b,c])*
        sp.sympify(right[j,c,d])*sp.sympify(boson[i,j])
        for i in rows for j in cols for b in range(m) for c in range(m)))


def run(output):
    rng = np.random.default_rng(87)
    def raw(shape):
        return rng.integers(-4, 5, size=shape)+1j*rng.integers(-4, 5, size=shape)
    left, right = raw((4,2,2)), raw((4,2,2))
    electron, other = raw((2,2)), raw((2,2))
    dg, dl = raw((4,4)), raw((4,4))
    d = 2
    def symraw(x):
        return np.array([sp.Integer(int(z.real))+sp.I*sp.Integer(int(z.imag))
                         for z in x.flat], dtype=object).reshape(x.shape)
    sl, se, sr, sd = map(symraw, (left, electron, right, dg))
    full_exact = exact(sl, se, sr, sd, range(4), range(4))
    ionic_exact = exact(sl, se, sr, sd, range(d), range(d))
    blocks_exact = (exact(sl, se, sr, sd, range(d), range(d,4))
                    + exact(sl, se, sr, sd, range(d,4), range(d))
                    + exact(sl, se, sr, sd, range(d,4), range(d,4)))
    if any(sp.expand(x) != 0 for x in full_exact-ionic_exact-blocks_exact):
        raise AssertionError('exact block identity failed')
    pair = EdgePair(left, right, left, right, electron, other, other,
                    electron, raw((2,2)), raw((2,2)), dg, dl)
    ion = replace(pair, left_t=left[:d], right_t=right[:d],
                  left_s=left[:d], right_s=right[:d],
                  d_greater=dg[:d,:d], d_lesser=dl[:d,:d])
    full = evaluate_components(pair, n=1, t=1, s=0)
    ionic = evaluate_components(ion, n=1, t=1, s=0)
    errors, mutations = {}, {}
    for order, g, boson in (('greater', electron, dg), ('lesser', other, dl)):
        actual = sigma_correction(left, g, right, boson, d)
        errors[order] = float(np.max(np.abs(actual-(full['sigma_'+order]-ionic['sigma_'+order]))))
        omitted = sigma_rectangular(left[d:], g, right[d:], boson[d:,d:])
        wrong_order = sigma_correction(left, g, right.swapaxes(-1,-2), boson, d)
        wrong_cross = sigma_correction(left, g, right, boson.T, d)
        mutations[order] = dict(omitted_cross=float(np.linalg.norm(actual-omitted)),
            wrong_reverse_band_order=float(np.linalg.norm(actual-wrong_order)),
            transposed_boson=float(np.linalg.norm(actual-wrong_cross)))
    # An exactly representable tiny correction beside a huge II term: the
    # direct blocks must retain what subtraction rounds out in binary64.
    l = np.ones((2,1,1), np.complex128)
    r = np.ones_like(l)
    g = np.ones((1,1), np.complex128)
    boson = np.array([[2.**60,1.],[2.,4.]], np.complex128)
    stable = sigma_correction(l,g,r,boson,1)[0,0]
    naive = (sigma_rectangular(l,g,r,boson)
             - sigma_rectangular(l[:1],g,r[:1],boson[:1,:1]))[0,0]
    if (max(errors.values()) > 1e-10 or
            min(v for group in mutations.values() for v in group.values()) <= 1e-9 or
            stable != 7j or naive == stable):
        raise AssertionError(dict(errors=errors, mutations=mutations,
                                  stable=str(stable), naive=str(naive)))
    result = dict(schema='p254-a87-joint-block-oracle-v1', status='passed',
        exact_sympy_identity=True, owning_A41_full_minus_ionic_errors=errors,
        rejected_mutations=mutations,
        finite_sensitive=dict(exact=[0,7], direct=[stable.real,stable.imag],
                              subtractive=[naive.real,naive.imag]),
        scope='Only new rectangular IF/FI/FF identity; no old physical suite',
        implementation=[binding(Path(__file__)), binding(Path(__file__).with_name('joint_blocks.py'))])
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    save(output, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), indent=2))
