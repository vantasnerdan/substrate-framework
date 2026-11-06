"""Authenticate coherent readout provenance without rewriting producer receipts."""
from hashlib import sha256
import json
from pathlib import Path

from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import file_binding

HERE = Path(__file__).resolve().parent
COMMITTED_SCHEMA = 'P254-committed-joint-q-readout-v1'


def binding_key(binding):
    return str(Path(binding['path']).resolve()), int(binding['bytes']), binding['sha256']


def same_binding(a, b):
    return binding_key(a) == binding_key(b)


def authenticate(binding, snapshots=()):
    """Accept original bytes, or an explicitly bound byte-identical source snapshot."""
    path = Path(binding['path'])
    if path.is_file() and same_binding(file_binding(path), binding):
        return
    matches = [item['snapshot'] for item in snapshots
               if same_binding(item['original'], binding)]
    if len(matches) != 1:
        raise ValueError('bound source bytes differ without one exact snapshot: '+str(path))
    snapshot = matches[0]
    if (int(snapshot['bytes']) != int(binding['bytes'])
            or snapshot['sha256'] != binding['sha256']
            or not same_binding(file_binding(Path(snapshot['path'])), snapshot)):
        raise ValueError('historical source snapshot bytes differ: '+str(path))


def canonical(value):
    return json.dumps(value, sort_keys=True, allow_nan=False)


def source_key(bindings):
    keys = [binding_key(b) for b in bindings]
    if not keys or len({k[0] for k in keys}) != len(keys):
        raise ValueError('unique explicit producer source bindings required')
    return sorted(keys)


def committed_seed(record, request):
    """Validate the native consumer and its historical forward producer separately.

    Return only the coherent kick/value, not a synthetic material-v2 receipt.
    Checkpoint headers and proof receipts are read; large propagated arrays are
    not consumed again by this scalar-only spatial reader.
    """
    if (record.get('schema') != COMMITTED_SCHEMA or record.get('status') != 'complete'
            or record.get('complete') is not True or 'error' in record
            or 'quantum' in record or 'coherent' in record
            or request.get('mode') != 'forward'
            or record.get('observable_selection') != 'total-force'
            or record.get('units') != 'Ry/bohr'
            or record.get('observable_metadata') is not None
            or record.get('fixed_original_D') is not True):
        raise ValueError('complete native committed total-force readout required')
    historical = record['historical_producer_request']
    if (historical['metadata'] != request or bound_json(historical['binding']) != request
            or request['identity'] != record['historical_producer_identity']):
        raise ValueError('historical forward request or producer identity differs')
    consumer_request = record['request']
    for key, expected in (('checkpoint', record['committed_checkpoint']),
                          ('historical_producer_request', historical['binding']),
                          ('physical_operator_oracle', record['physical_operator_oracle']),
                          ('readout_oracle', record['readout_oracle'])):
        if not same_binding(consumer_request[key], expected):
            raise ValueError('native consumer input binding differs: '+key)
    if (consumer_request['schema'] != COMMITTED_SCHEMA
            or consumer_request['historical_producer_identity'] != request['identity']
            or consumer_request['consumer_identity'] != record['consumer_identity']):
        raise ValueError('native consumer and historical producer identities differ')
    checkpoint = bound_json(record['committed_checkpoint'])
    if (checkpoint.get('schema') != 'streamed-joint-response-v1'
            or checkpoint.get('complete') is not True or checkpoint.get('kind') != 'forward'
            or checkpoint['identity'] != request['identity'] or checkpoint['q'] != record['q']
            or checkpoint['time'] != record['time_ry_inverse']
            or checkpoint['time'] != request['times'][-1]
            or checkpoint['provenance'] != record['checkpoint_provenance']
            or not same_binding(checkpoint['provenance']['baseline'], record['baseline'])):
        raise ValueError('native scalar is not bound to this complete forward checkpoint')
    operator = consumer_request['operator']
    for key in ('case', 'initialization', 'model_identity', 'q', 'atom', 'axis', 'dynamics'):
        if operator[key] != record[key]:
            raise ValueError('native operator selection differs: '+key)
    if (operator['observable'] != 'total-force' or operator['field_axis'] is not None
            or operator['observable_metadata'] is not None
            or operator['probe_identity'] != record['observable']
            or operator['time'] != record['time_ry_inverse']
            or operator['block_size'] != request['block_size']
            or not same_binding(operator['history'], record['baseline'])
            or operator['specification']['scalar'] != record['scalar']
            or not same_binding(operator['specification']['aggregate'], record['matched_aggregate'])
            or sha256(canonical(operator['specification']).encode()).hexdigest() != record['model_identity']
            or source_key(operator['sources']) != source_key(record['source_bindings'])):
        raise ValueError('native scalar/operator/model/source contract differs')
    identity = sha256(canonical(dict(model=record['model_identity'], q=record['q'],
        dynamics=record['dynamics'], sources=record['source_bindings'])).encode()).hexdigest()
    if identity != record['consumer_identity']:
        raise ValueError('native consumer source identity differs')
    for key, schema in (('physical_operator_oracle', 'P254-joint-partitioned-q-oracle-v1'),
                        ('readout_oracle', 'P254-joint-q-readout-oracle-v1')):
        proof = bound_json(record[key])
        if proof.get('schema') != schema or proof.get('status') != 'passed':
            raise ValueError('native bound operator/readout proof is not passed')
    return dict(readout=record['readout'], amplitude=request['amplitude'], cell=request['cell'])


def validate_committed_split(record, request, pair, pair_request):
    seed = committed_seed(record, request)
    if (pair_request.get('mode') != 'pair' or pair.get('status') != 'complete'
            or 'error' in pair or 'coherent' in pair):
        raise ValueError('native coherent readout requires a genuine successful pair receipt')
    for key in ('case', 'model_identity', 'scalar', 'q', 'dynamics', 'observable',
                'all_explicit_edges', 'fixed_original_D'):
        if record[key] != pair[key]:
            raise ValueError('native forward/pair scientific selection differs: '+key)
    for key in ('initialization', 'baseline', 'matched_aggregate'):
        if not same_binding(record[key], pair[key]):
            raise ValueError('native forward/pair original input differs: '+key)
    if request['times'] != pair['time_grid'] or request['second_times'] != pair['second_time_grid']:
        raise ValueError('native forward/pair integration grids differ')
    ignored = {'identity', 'mode', 'left_checkpoint', 'left_checkpoint_sha256',
               'right_checkpoint', 'right_checkpoint_sha256'}
    if ({key: value for key, value in request.items() if key not in ignored}
            != {key: value for key, value in pair_request.items() if key not in ignored}):
        raise ValueError('native forward/pair scientific request or original kick differs')
    # These three roles only read a completed forward state. The material-v2
    # pair body never calls them. All propagation, covariance, preparation and
    # producer-wrapper sources must still match exactly; no general source
    # equivalence or old/current producer-identity substitution is admitted.
    readout_only = {str(HERE/name) for name in (
        'partitioned_material_readout.py', 'run_joint_q_readout.py',
        'mixed_material_driver_support.py')}
    physical = lambda bindings: source_key([b for b in bindings
                                            if str(Path(b['path']).resolve()) not in readout_only])
    if physical(record['source_bindings']) != physical(pair['source_bindings']):
        raise ValueError('native forward/pair physical producer sources differ')
    return seed
