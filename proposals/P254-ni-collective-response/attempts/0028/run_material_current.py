#!/usr/bin/env python3
"""Bounded full-row optical current/contact consumer of accepted A25 data.

CLI reads the source-registration bindings, restores accepted background
arrays verbatim, loads the original spin-1 Gamma basis and acts on at most
eight explicitly selected saved columns. No Native backend or DFT is loaded.
Root may instead import consume(hamiltonian, vectors, ...) to reuse its already
restored accepted background and saved/reconstructed full-PW factors directly.
The output is unweighted action data, NOT a populated material current/noise.
"""
from pathlib import Path
import argparse
import io
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
for directory in (ROOT / 'src', HERE.parent / '0025'):
    sys.path.insert(0, str(directory))

import numpy as np
from substrate_framework.qe_material_current import MaterialOpticalCurrent
from material_run_contract import runtime_context, implementation


def consume(hamiltonian, vectors, *, direction, step_bohr_inv):
    """Return raw arrays + measured refinement, bounded to eight full-row columns.

    All three levels are saved, not just a preferred value. Deltas are reported
    without assigning an unjustified tolerance or extrapolating across cubic
    knots. Both derivative/contact must settle before physical interpretation.
    """
    vectors = np.asarray(vectors)
    if vectors.ndim != 2 or not 1 <= vectors.shape[1] <= 8:
        raise ValueError('one through eight explicit full-row columns required')
    operator = MaterialOpticalCurrent(hamiltonian)
    arrays = dict(input_columns=vectors, direction=np.asarray(direction, dtype=float))
    levels = []
    previous = None
    for index, (step, linear, contact) in enumerate(operator.refinement_actions(
            vectors, direction, step_bohr_inv=step_bohr_inv)):
        arrays[f'linear_ry_bohr_level{index}'] = linear
        arrays[f'contact_ry_bohr2_level{index}'] = contact
        arrays[f'overlap_linear_level{index}'] = vectors.conj().T @ linear
        arrays[f'overlap_contact_level{index}'] = vectors.conj().T @ contact
        # Full-row ordered product, not a product of band-projected matrices.
        arrays[f'ordered_linear_product_level{index}'] = linear.conj().T @ linear
        entry = dict(step_bohr_inv=step, linear_norm=float(np.linalg.norm(linear)),
                     contact_norm=float(np.linalg.norm(contact)))
        if previous is not None:
            entry.update(linear_delta=float(np.linalg.norm(linear-previous[0])),
                         contact_delta=float(np.linalg.norm(contact-previous[1])))
        levels.append(entry)
        previous = linear, contact
    shift = step*np.asarray(direction)
    zero = operator.field_action(vectors, np.zeros(3))
    plus = operator.field_action(vectors, shift)
    minus = operator.field_action(vectors, -shift)
    direct_linear = (plus-minus)/(2*step)
    direct_contact = ((plus-zero)+(minus-zero))/step**2
    arrays.update(direct_h_linear_ry_bohr=direct_linear, direct_h_contact_ry_bohr2=direct_contact)
    summary = dict(schema='P254-material-optical-current-v1',
        scope='selected unweighted columns, all physical PW rows; optical q=0 only',
        ispin=hamiltonian.ispin, npw=vectors.shape[0], ncolumn=vectors.shape[1],
        direction=np.asarray(direction).tolist(), refinement=levels,
        finest_direct_H_linear_delta=float(np.linalg.norm(linear-direct_linear)),
        finest_direct_H_contact_delta=float(np.linalg.norm(contact-direct_contact)),
        convergence='not automatically certified; inspect both action refinement series and interpolation knots',
        field_family='H_Ry(a)=(p+a)^2+V_saved+sum B(p+a) D B(p+a)^dagger; atomic placement phase fixed',
        units=dict(linear='Ry bohr', contact='Ry bohr^2', shift='bohr^-1'),
        charge_current_atomic_units='-linear/2; no k weight, state contraction or volume included',
        sources=list(hamiltonian.factors.sources),
        limitations=['not finite-q transverse coupling', 'not an electromagnetic mode or loss model',
                    'fixed saved local potential and finite PW rows',
                    'no physical temperature, effective ionic charge, linewidth or cavity supplied',
                    'piecewise cubic interpolation may not possess a second derivative at knots'])
    return arrays, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, default=HERE/'optical-current-sources.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--first-column', type=int, default=0, help='zero-based saved column')
    parser.add_argument('--columns', type=int, default=2)
    parser.add_argument('--direction', type=float, nargs=3, default=(1., 0., 0.))
    parser.add_argument('--step-bohr-inv', type=float, required=True)
    args = parser.parse_args()
    if args.first_column < 0 or not 1 <= args.columns <= 8:
        parser.error('nonnegative first column and 1..8 columns required')
    args.output.mkdir(parents=True, exist_ok=False)
    context = runtime_context()
    with (args.output/'runtime-context.json').open('x') as stream:
        json.dump(context, stream, indent=2)
    from substrate_framework.qe_saved_native import read_bound_bytes, load_wavefunctions
    from substrate_framework.qe_material_reference import load_material_inputs
    from substrate_framework.qe_material_hamiltonian import build_pw_factors, MaterialHamiltonian
    from material_reference_output import restore_background
    registration = json.loads(args.registration.read_text())
    binding = registration['accepted_material']
    metadata_bytes, _ = read_bound_bytes(binding['background_metadata'])
    metadata = json.loads(metadata_bytes)
    archive_bytes, _ = read_bound_bytes(metadata['archive'])
    with np.load(io.BytesIO(archive_bytes), allow_pickle=False) as archive:
        background_arrays = {key: archive[key] for key in archive.files}
    del archive_bytes
    sources = binding['inputs']
    inputs = load_material_inputs(sources['xml'], sources['density'], sources['upfs'])
    background = restore_background(inputs, background_arrays)
    del background_arrays
    wave = load_wavefunctions(binding['wavefunction'])
    stop = args.first_column+args.columns
    if stop > wave.nbnd:
        parser.error('selected columns exceed original saved band count')
    factors = build_pw_factors(inputs, wave)
    hamiltonian = MaterialHamiltonian(background, factors, wave.ispin, block_size=args.columns)
    arrays, summary = consume(hamiltonian, wave.coefficients[:, args.first_column:stop],
                              direction=args.direction, step_bohr_inv=args.step_bohr_inv)
    summary.update(original_wavefunction=dict(wave.source),
                   selected_columns_zero_based=list(range(args.first_column, stop)),
                   accepted_background_metadata=binding['background_metadata'])
    summary['runtime_context'] = context
    summary['implementation_sources'] = implementation((Path(__file__),))
    # Exclusive publication avoids silently overwriting a previous numerical run.
    with (args.output/'optical-current.npz').open('xb') as stream:
        np.savez(stream, **arrays)
    with (args.output/'optical-current.json').open('x') as stream:
        json.dump(summary, stream, indent=2)
        stream.write('\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
