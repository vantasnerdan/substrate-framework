"""Raw-byte provenance for the new matched clamped consumers and finite oracle."""
from pathlib import Path

from material_run_contract import binding, bound_json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REGISTRATION = HERE/'matched-clamped-registration.json'


def source_bindings():
    registration = bound_json(binding(REGISTRATION))
    for record in registration['settled_sources']:
        if binding(Path(record['path'])) != record:
            raise ValueError('settled matched-clamped source changed: '+record['path'])
    paths = [REGISTRATION, *[HERE/name for name in registration['new_files']]]
    paths += [HERE.parent/'0025'/name for name in
              ('material_run_contract.py', 'material_reference_output.py')]
    paths += [ROOT/'src/substrate_framework'/name for name in
              ('material_modes.py', 'material_control_reference.py',
               'material_force_lattice.py', 'physical_bilinear_preparation.py',
               'inertial_fermion_dynamics.py')]
    return registration['settled_sources']+[binding(path) for path in paths]


def require_oracle(path):
    sources = source_bindings()
    record = binding(path)
    result = bound_json(record)
    if (result.get('schema') != 'P254-matched-clamped-oracle-v1'
            or result.get('status') != 'passed' or result.get('source_bindings') != sources):
        raise ValueError('passed new-join oracle bound to these exact source bytes required')
    return record, sources
