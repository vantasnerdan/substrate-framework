"""Genuine public saved-output numeric consumer for coupled-driver handoff."""
import argparse
import json
from pathlib import Path
import numpy as np
from supplied_vie2 import load_solution
from command_receipts import binding


def consume(path):
    result = load_solution(path)
    factor = -1j if int(result['statistics']) == -1 else 1j
    prepared = factor * result['less'][0,0]
    final = factor * result['less'][-1,-1]
    return {'input':binding(path),'statistics':int(result['statistics']),
            'times':result['times'].tolist(),'initial_trace':float(np.trace(prepared).real),
            'final_trace':float(np.trace(final).real),
            'final_matrix_real':final.real.tolist(),'final_matrix_imag':final.imag.tolist(),
            'causal_upper_retarded_max':float(np.max(np.abs(result['ret'][np.triu_indices(len(result['times']),1)]))),
            'scope':'finite supplied-kernel numerical output, not NiH/material calibration'}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('solution')
    parser.add_argument('output')
    args = parser.parse_args()
    report = consume(args.solution)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
