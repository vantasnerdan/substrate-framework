#!/usr/bin/env python3
"""Expose column-order, complex-D and whole-output memory errors in paneling."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
import numpy as np
from scipy.linalg import expm
import sympy as sp
from joint_column_mean import ColumnLimits
from moving_exponential_panels import panel_plan, bounded_panel_columns
from substrate_framework.column_frame_krylov import memory_bound_bytes
from substrate_framework.streamed_material_state import file_binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    epsilon, a, b, c = sp.symbols('epsilon a b c', positive=True)
    allocation = sp.simplify(sum((epsilon*sp.sqrt(w/(a+b+c)))**2 for w in (a,b,c))-epsilon**2)
    assert allocation == 0
    rng = np.random.default_rng(225392)
    h = rng.normal(size=(7,7))+1j*rng.normal(size=(7,7))
    frame = rng.normal(size=(7,11))+1j*rng.normal(size=(7,11))
    frame[:,3] = (.2+.3j)*frame[:,0]+.7*frame[:,1]
    frame[:,7] = 0
    raw_d = rng.normal(size=(11,11))+1j*rng.normal(size=(11,11))
    workspace = frame.nbytes+memory_bound_bytes(7,3,7)
    limits = ColumnLimits(7,workspace,1e-9)
    actual, receipt = bounded_panel_columns(lambda x:h@x,frame,.017,limits=limits,
        norm_bound=float(np.linalg.norm(h,2)),
        growth_bound=float(np.linalg.norm((h-h.conj().T)/(2j),2)),identity='independent-dense-nonHermitian')
    expected = expm(-1j*.017*h)@frame
    error = float(np.linalg.norm(actual-expected))
    density_error = float(np.max(np.abs(actual@raw_d@actual.conj().T-expected@raw_d@expected.conj().T)))
    wrong = actual[:,::-1]@raw_d@actual[:,::-1].conj().T
    mutation = float(np.max(np.abs(wrong-expected@raw_d@expected.conj().T)))
    assert error < limits.absolute_error_target and density_error < 1e-7 and mutation > 1
    assert np.array_equal(actual[:,7],np.zeros(7))
    assert receipt['panel_plan']['panel_count'] > 1
    native_limits = ColumnLimits(6,1073741824,1e-5)
    native = panel_plan((18613,346),native_limits)
    assert native['maximum_numerical_bytes'] <= native_limits.workspace_bytes
    output_bytes = 16*18613*346
    wrong_width = max(w for w in range(1,347) if memory_bound_bytes(18613,w,6)<=native_limits.workspace_bytes)
    wrong_memory = output_bytes+memory_bound_bytes(18613,wrong_width,6)
    assert wrong_memory > native_limits.workspace_bytes
    result = dict(schema='P254-independent-column-panel-oracle-v1',status='passed',
        exact_error_allocation_residual=str(allocation),frame_error=error,
        complex_raw_D_density_error=density_error,reordered_column_density_separation=mutation,
        finite_panel_plan=receipt['panel_plan'],native_panel_plan=native,
        omitted_full_output_memory_mutation_bytes=wrong_memory,
        source_bindings=[file_binding(path) for path in (Path(__file__),HERE/'moving_exponential_panels.py',
            HERE/'gramian_column_mean.py',HERE/'gramian_krylov_bounds.py',
            HERE/'joint_column_mean.py',ROOT/'src/substrate_framework/column_frame_krylov.py')],
        scope='Finite independent matrix-function/complex-D and actual-shape workspace evidence; not Native propagation, RSS or interval certification')
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(result))


if __name__=='__main__':
    main()
