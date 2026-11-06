"""Bind already-paid original bare factors; never produce U, means or controls."""
from pathlib import Path
import argparse
import json
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/'0063')]
from history_storage import authenticate, binding
from saved_mean import SavedMean, BARE_SCHEMA, REFERENCE


def produce(output):
    output = Path(output).resolve()
    prereg = json.loads((HERE/'preregistration.json').read_bytes())
    expected = {item['path']: item for item in prereg['verified_raw_input_bindings']}
    bundle = expected[prereg['genuine_reference']]
    archive = expected[prereg['bare_source_factors']]
    for item in (bundle, archive):
        authenticate(item)
    result = Path(archive['path']).parent/'result.json'
    paid = json.loads(result.read_bytes())
    if paid['bare_source_factors'] != archive:
        raise ValueError('paid original factors differ from frozen source')
    spec = dict(schema=BARE_SCHEMA, status='complete', reference_source_kind=REFERENCE,
        reference_bundle=bundle, factor_result=binding(result), archive=archive,
        times_inverse_ry=paid['times_inverse_ry'],
        source_factor_semantics=paid['source_factor_semantics'],
        represented_control=None,
        interpretation='Literal prescribed bare zero-iterate; no common mean feedback or Native history solve.')
    with output.open('x') as stream:
        json.dump(spec, stream, indent=2)
        stream.write('\n')
    with SavedMean(output) as source:
        if not source.is_bare or source.reference_identity != bundle:
            raise ValueError('saved source identity differs')
    return binding(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    print(json.dumps(produce(parser.parse_args().output)), flush=True)
