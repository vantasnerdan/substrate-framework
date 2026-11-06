"""Compound A25/A27 authority and completed-boundary recovery for Bloch pilot."""
from pathlib import Path
from material_run_contract import bound_json, recover_block
from material_reference_output import load_arrays
from .material_dyad_output import load_accepted_background, REFERENCE_SHA256
from .qe_bloch_grid import UNITARY_POLICY

SCHEMA = 'material-bloch-coverage-v1'
CONTINUATION = 'cbb9f3698d8b432fb3825c640487095f9abef15004d36aeb1c311a38174cb92e'
IMPLEMENTATION = 'b59c25bd13d015ae17b2d81c8b98567f1358bacdc783b5c235f2f186aef268bc'
ORACLE = 'b66f8b8c7f0921340ad2a1b6e9b34d2a9f96b7d4b4a43e8c601f040bf6ddebb9'
PROBE = 'first-folded-spin1-bands1-2-H32x+iNi1y-number-cosb1+isinb2-dyad1-2-v1'
RECEIPT_FIELDS = ('recipe_sha256','archive','arrays','metadata','readback')


def wave_bindings(inputs):
    return {(int(spin),int(ik)): value for spin, table in inputs['wavefunctions'].items() for ik,value in table.items()}


def validate_recipe(recipe):
    inputs, settings = recipe['inputs'], recipe['settings']
    if set(inputs) != {'xml','density','upfs','wavefunctions','reference_terminal'}:
        raise ValueError('exact bound material and accepted background inputs required')
    if set(inputs['upfs']) != {'Ni','H'} or set(inputs['wavefunctions']) != {'1','2'}:
        raise ValueError('both original species and explicit spins required')
    if any(set(table) != {str(i) for i in range(1,11)} for table in inputs['wavefunctions'].values()):
        raise ValueError('all20 original WFC bindings required')
    reference_binding = inputs['reference_terminal']
    if reference_binding['sha256'] != REFERENCE_SHA256 or str(reference_binding['bytes']) != '1000830':
        raise ValueError('accepted original A25 background terminal required')
    reference_recipe = bound_json(reference_binding)['recipe']
    reduced = background_inputs(inputs)
    if any(reduced[key] != reference_recipe['inputs'][key] for key in ('xml','density','upfs','wavefunctions')):
        raise ValueError('source admission inputs differ from accepted original A25 material')
    if set(settings) != {'continuation_registration','implementation_registration','oracle_receipt','symmetry_policy','probe_recipe'}:
        raise ValueError('exact continuation/implementation/oracle/policy/probe settings required')
    for key,digest in (('continuation_registration',CONTINUATION),('implementation_registration',IMPLEMENTATION),('oracle_receipt',ORACLE)):
        if settings[key]['sha256'] != digest:
            raise ValueError('wrong settled A27 authority: '+key)
        bound_json(settings[key])
    implementation = bound_json(settings['implementation_registration'])
    spatial = next(item for item in implementation['settled_inputs'] if Path(item['path']).name == 'bloch-spatial-grounding.json')
    source_report = bound_json(spatial)
    table = source_report['source_access']['existing_byte_closure']['wfc_identity_table']
    for item in table:
        for spin, key in ((1,'up_sha256'),(2,'down_sha256')):
            original = inputs['wavefunctions'][str(spin)][str(item['ik'])]
            if original['sha256'] != item[key] or str(original['bytes']) != item['bytes_each']:
                raise ValueError('WFC binding differs from settled original source closure')
    if (settings['symmetry_policy'] != UNITARY_POLICY or settings['probe_recipe'] != PROBE
            or recipe.get('library_binding') is not None):
        raise ValueError('registered unitary/probe recipe, no Native backend, required')


def validate_grant(grant):
    validate_recipe(grant)
    allowed = {f'source-{spin}-{ik}' for spin in (1,2) for ik in range(1,11)} | {'coverage-probe'}
    receipts = grant['resume_receipts']
    if not isinstance(receipts,dict) or not set(receipts) <= allowed or (not grant['resume'] and receipts):
        raise ValueError('explicit completed-boundary receipt bindings required; fresh uses {}')


def background_inputs(inputs):
    return {**inputs,'wavefunctions':{spin:table['1'] for spin,table in inputs['wavefunctions'].items()}}


def restore_accepted(inputs):
    return load_accepted_background(background_inputs(inputs))


def source_metadata(recipe, spin, ik):
    return dict(schema=SCHEMA, boundary='complete-original-source', ispin=spin, source_ik=ik,
                inputs=recipe['inputs'], settings=recipe['settings'],
                original_binding=recipe['inputs']['wavefunctions'][str(spin)][str(ik)],
                coefficient_storage='original bound WFC only; all columns identity retained',
                cutoff_policy='strict-interior complete; original native boundary retained; 1e-8 Ry')


def final_metadata(recipe, background_receipt):
    return dict(schema=SCHEMA,boundary='complete-coverage-and-selected-probe', inputs=recipe['inputs'],
                settings=recipe['settings'],source_count=20,slot_count=128,
                background_receipt=background_receipt,
                scope='complete lazy spin/k representation; selected folded nonGamma q=0 raw probe, not preparation',
                units=dict(h0='Ry',ionic='Ry/bohr',number='dimensionless',density='bohr^-3',density_gradient='bohr^-4'),
                state=None,physical_temperature=None,finite_q=None)


def recover_authorized(run, key, metadata):
    """Only an exact caller-bound completed receipt licenses cache consumption."""
    target = run.target
    exists = (target/(key+'.json')).exists() or (target/(key+'.npz')).exists()
    admitted = run.grant['resume_receipts']
    if not exists:
        if key in admitted:
            raise ValueError('admitted completed checkpoint is absent: '+key)
        return None
    if not run.grant['resume'] or key not in admitted:
        raise ValueError('existing checkpoint lacks explicit resume authority: '+key)
    if Path(admitted[key]['path']).resolve() != target/(key+'.json'):
        raise ValueError('resume receipt path differs')
    expected = bound_json(admitted[key])
    if expected['metadata'] != metadata:
        raise ValueError('resume metadata differs from current immutable recipe')
    actual = recover_block(target,key,run.identity,metadata)
    if actual != expected:
        raise ValueError('resume archive receipt differs from explicit authority')
    run.record(key,actual,'completed-cache-reuse')
    return load_arrays(actual), actual
