"""Paid-table number-channel contraction; not an Ndot/error enclosure.

Compare Sigma[G_previous] with its own input and the newly solved output.
The finite exact checks establish the real-part cancellation/sign convention
for Hermitian scalar electron-boson and two-spin local Born examples only.
They do not re-prove the admitted full material functional or the time solver.
"""
from pathlib import Path
import argparse
import json
import sys
import numpy as np
from scipy.linalg.blas import zdotu
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0056'))
from dressed_reference import greater, authenticate, binding, save


def finite_convention():
    ar, ai, br, bi, cr, ci, gt, gs = sp.symbols('ar ai br bi cr ci gt gs', real=True)
    a, b, d = ar+sp.I*ai, br+sp.I*bi, cr+sp.I*ci
    sg = sp.I*gt*a*gs*d
    sl = sp.I*gt*b*gs*(-sp.conjugate(d))
    value = sg*(-sp.conjugate(b))-sl*(-sp.conjugate(a))
    exact = sp.simplify(sp.re(value))
    assert exact == 0
    wrong_reverse = sp.simplify(sp.re(sg*sp.conjugate(b)-sl*(-sp.conjugate(a))))
    witness = wrong_reverse.subs(dict(zip((ar, ai, br, bi, cr, ci, gt, gs), (1, 2, 3, 5, 7, 11, 2, 3))))
    assert witness != 0
    a1, b1 = cr+sp.I*ci, gt+sp.I*gs
    ee_g = -a*a1*sp.conjugate(b1)
    ee_l = -b*b1*sp.conjugate(a1)
    ee = ee_g*(-sp.conjugate(b))-ee_l*(-sp.conjugate(a))
    ee_exact = sp.simplify(sp.re(ee))
    assert ee_exact == 0
    return dict(electron_boson_real_contraction=str(exact),
        two_spin_local_born_real_contraction=str(ee_exact),
        reversed_lesser_sign_mutation_witness=str(witness),
        scope='Finite exact sign examples, not a full material Ward or time-accuracy proof.')


def contraction(sigma_g, sigma_l, dataset):
    with np.load(dataset, allow_pickle=False) as data:
        ret, less, times = data['electron_ret'], data['electron_less'], data['times']
    nt = len(times)
    result = np.empty((nt, nt), complex)
    crosscheck = 0.
    for a in range(nt):
        for b in range(nt):
            reverse_l = less[b, a]
            reverse_g = greater(ret, less, b, a)
            value = zdotu(sigma_g[a, b].ravel(), reverse_l.swapaxes(-1, -2).ravel())
            value -= zdotu(sigma_l[a, b].ravel(), reverse_g.swapaxes(-1, -2).ravel())
            direct = (np.einsum('kij,kji->', sigma_g[a, b], reverse_l)
                      -np.einsum('kij,kji->', sigma_l[a, b], reverse_g))
            crosscheck = max(crosscheck, float(abs(value-direct)))
            result[a, b] = value
    return dict(real=result.real.tolist(), imag=result.imag.tolist(),
        max_real_absolute=float(np.max(np.abs(result.real))),
        max_imag_absolute=float(np.max(np.abs(result.imag))),
        independent_bilinear_trace_crosscheck_maxabs=crosscheck)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--components', type=Path, required=True)
    p.add_argument('--output-state', type=Path, required=True)
    p.add_argument('--result', type=Path, required=True)
    args = p.parse_args()
    oracle = finite_convention()
    catalog = json.loads(args.components.read_bytes())
    old = authenticate(catalog['input_dataset'])
    with np.load(authenticate(catalog['archive']), allow_pickle=False) as data:
        sg, sl = data['sigma_greater'], data['sigma_lesser']
    result = dict(schema='p254-paid-number-channel-contraction-v1',
        definition='sum_slot Tr[Sigma_greater(t,s) G_lesser(s,t) - Sigma_lesser(t,s) G_greater(s,t)]; real/imag raw, nonconjugating trace.',
        components=binding(args.components), input_state=binding(old),
        output_state=binding(args.output_state), implementation=binding(__file__),
        finite_exact_convention=oracle,
        matched_functional_input=contraction(sg, sl, old),
        newly_solved_output=contraction(sg, sl, args.output_state),
        interpretation='Distinguishes matched-functional algebra from an unclosed state update; not a discrete continuity equation, converged fixed point, temporal accuracy or physical loss.')
    save(args.result, result)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
