"""Root-only execution entrypoint: physical factors, CAR, SymPy, material bridge."""
import argparse
import json
from pathlib import Path
from time import perf_counter
# Reuse the existing A29/A41 import roots, without loading material archives.
from coulomb_material_inputs import HERE, file_binding
from oracle_coulomb_factors import run as factors
from oracle_second_born_fock import run as fock
from oracle_second_born_sympy import run as symbolic
from oracle_material_contraction import run as material


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    started = perf_counter()
    results = [check() for check in (factors,fock,symbolic,material)]
    receipt = dict(schema='p254-a44-independent-oracles-v1',results=results,
        seconds=perf_counter()-started,
        implementation=[file_binding(p) for p in sorted(HERE.glob('*.py'))],
        scope='executed encoded finite-model identities, not actual material rate/dynamics or source-independent physical adequacy')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(receipt,stream,indent=2,allow_nan=False); stream.write('\n')
    print(json.dumps(dict(receipt=file_binding(args.output),seconds=receipt['seconds'])),flush=True)


if __name__ == '__main__':
    main()
