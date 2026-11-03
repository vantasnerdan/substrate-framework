"""Authenticate genuine original and cumulative-continuation pair envelopes.

Both delegate the same physical pair calculation. A continuation retains its
own schema and grant provenance; it is never relabeled an original result.
"""
from substrate_framework.material_control_reference import bound_json
from substrate_framework.partitioned_joint_slots import digest
from joint_spatial_provenance import authenticate, same_binding, source_key

MATERIAL_SCHEMA = 'P254-source-local-q-material-v1-result'
CONTINUATION_SCHEMA = 'P254-source-local-material-cumulative-continuation-v1-result'
PAIR_SCHEMAS = frozenset((MATERIAL_SCHEMA, CONTINUATION_SCHEMA))
ALLOCATION_SCHEMA = 'P254-source-local-material-cumulative-continuation-v1-allocation'


def validate_pair_identity(pair, owner, snapshots):
    if pair['schema'] == MATERIAL_SCHEMA:
        expected = digest(dict(schema='P254-source-local-q-material-v1', owner=pair['owner']))
        if pair['material_identity'] != expected:
            raise ValueError('source-local material owner identity differs')
        return
    if pair['schema'] != CONTINUATION_SCHEMA:
        raise ValueError('unknown source-local material pair producer')
    grant = bound_json(pair['allocation'])
    if (grant.get('schema') != ALLOCATION_SCHEMA
            or bound_json(pair['owned_allocation_copy']) != grant
            or not same_binding(grant['original_owner'], pair['owner'])
            or not same_binding(grant['previous_terminal'], pair['previous_terminal'])
            or grant['original_started_unix'] != owner['started_unix']
            or pair['original_started_unix'] != owner['started_unix']
            or pair['cumulative_caps'] != grant['caps']
            or pair.get('original_owner_rewritten') is not False
            or pair.get('resume') is not True):
        raise ValueError('continuation grant/original owner/clock differs')
    bound_json(pair['previous_terminal'])
    sources = pair['continuation_source_bindings']
    if source_key(sources) != source_key(grant['source_bindings']):
        raise ValueError('continuation producer sources differ from its bound grant')
    for source in sources:
        authenticate(source, snapshots)
