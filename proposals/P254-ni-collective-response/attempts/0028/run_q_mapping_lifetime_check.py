#!/usr/bin/env python3
"""Verify q mapping lifetime on the complete real seed without replaying physics."""
import argparse
import ast
import json
import os
from pathlib import Path
import re
import resource
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from run_joint_q_material import bindings, producer_identity
from substrate_framework.material_control_reference import bound_json
from substrate_framework.partitioned_joint_slots import array_identity
from substrate_framework.partitioned_q_storage import QBudget, bounded_array_identity, release_mapped_pages
from substrate_framework.streamed_material_state import authenticate, file_binding


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('retained-oracle', 'producer-result', 'checkpoint', 'propagation-request'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    return parser.parse_args()


def checked(path, digest):
    binding = file_binding(path)
    if binding['sha256'] != digest:
        raise ValueError('explicit input digest differs: '+str(path))
    return binding, bound_json(binding)


def unchanged_numerics(before, after, module):
    old, new = ast.parse(before), ast.parse(after)
    if module == 'partitioned_q_storage.py':
        helpers = [n for n in new.body if isinstance(n, ast.FunctionDef) and n.name in ('release_mapped_pages', 'bounded_array_identity')]
        if len(helpers) != 2:
            return False
        for node in helpers:
            new.body.remove(node)
        for node in list(new.body):
            if isinstance(node, ast.Import) and len(node.names) == 1 and node.names[0].name in ('mmap', 'numpy'):
                new.body.remove(node)
            if isinstance(node, ast.ImportFrom) and node.module == 'partitioned_joint_slots':
                node.names = [alias for alias in node.names if alias.name != 'array_identity']
    else:
        for node in new.body:
            if isinstance(node, ast.ImportFrom) and node.module == 'partitioned_joint_slots':
                node.names.insert(0, ast.alias(name='array_identity'))
            if isinstance(node, ast.ImportFrom) and node.module == 'partitioned_q_storage':
                node.names = [alias for alias in node.names if alias.name not in ('release_mapped_pages', 'bounded_array_identity')]
        release_loop = ast.dump(ast.parse('for value in x.values():\n    release_mapped_pages(value)').body[0], include_attributes=False)
        class Normalize(ast.NodeTransformer):
            loops = 0
            names = 0
            def visit_For(self, node):
                if ast.dump(node, include_attributes=False) == release_loop:
                    self.loops += 1
                    return None
                return self.generic_visit(node)
            def visit_Name(self, node):
                if node.id == 'bounded_array_identity':
                    node.id = 'array_identity'
                    self.names += 1
                return node
        normalization = Normalize()
        new = normalization.visit(new)
        if normalization.loops != 1 or normalization.names != 8:
            return False
    return ast.dump(old, include_attributes=False) == ast.dump(new, include_attributes=False)


def mapped_rss(value):
    address = value.ctypes.data
    selected = False
    for line in Path('/proc/self/smaps').read_text().splitlines():
        match = re.match(r'^([0-9a-f]+)-([0-9a-f]+) ', line)
        if match:
            selected = int(match[1], 16) <= address < int(match[2], 16)
        elif selected and line.startswith('Rss:'):
            return 1024*int(line.split()[1])
    raise ValueError('actual array mapping not present in own smaps')


def execute(args, budget):
    retained_binding, retained = checked(args.retained_oracle, args.retained_oracle_sha256)
    producer_binding, producer = checked(args.producer_result, args.producer_result_sha256)
    checkpoint_binding, checkpoint = checked(args.checkpoint, args.checkpoint_sha256)
    request_binding, request = checked(args.propagation_request, args.propagation_request_sha256)
    current = bindings()
    original_identity = producer_identity(producer['model_identity'], producer['q'],
                                         producer['dynamics'], retained['source_bindings'])
    if original_identity != checkpoint['identity']:
        raise ValueError('seed is not the original total-force q producer')
    old = {item['path']: item for item in retained['source_bindings']}
    changed = {item['path'] for item in current if old.get(item['path'], {}).get('sha256') != item['sha256']}
    modules = ('partitioned_q_storage.py', 'partitioned_q_response.py')
    expected = {str(ROOT/'src/substrate_framework'/name) for name in modules}
    if (retained['status'] != 'passed' or changed != expected
            or producer['source_bindings'] != retained['source_bindings']
            or producer['budget']['hamiltonian_actions'] != 0
            or checkpoint['kind'] != 'coherent-initial' or checkpoint['time'] != 0.
            or request['identity'] != checkpoint['identity'] or request['dual'] is not False):
        raise ValueError('exact prior q sources and genuine unpropagated coherent seed required')
    checks, prior_sources = [], []
    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        if not passed:
            raise AssertionError(name)
    for module in modules:
        source = ROOT/'src/substrate_framework'/module
        prior = HERE/'source-history'/(module.removesuffix('.py')+'.pre-mapped-release-001.py')
        binding = file_binding(prior)
        if binding['sha256'] != old[str(source)]['sha256']:
            raise ValueError('retained original numerical source differs')
        prior_sources.append(binding)
        check(module+'_numerics_unchanged_outside_page_lifetime', unchanged_numerics(prior.read_bytes(), source.read_bytes(), module))
    maps = []
    for binding in checkpoint['legs']:
        budget.check()
        authenticate(binding)
        maps.append(np.load(binding['path'], mmap_mode='r', allow_pickle=False))
    # Checkpoint files interleave ket/bra; values() and request identities do not.
    maps = maps[::2]+maps[1::2]
    original_identity = array_identity(maps[0])
    before_rss = mapped_rss(maps[0])
    released_identity = bounded_array_identity(maps[0])
    after_rss = mapped_rss(maps[0])
    check('actual_identity_unchanged_and_no_release_mutation_detected',
          original_identity == released_identity and before_rss > 0 and after_rss == 0)
    array_identity(maps[0])
    release_mapped_pages(maps[0])
    check('retained_original_rk4_input_releases_resident_pages', mapped_rss(maps[0]) == 0)
    identities = []
    for index, value in enumerate(maps):
        budget.check()
        identities.append(bounded_array_identity(value))
        budget.check_processes(os.getpid(), [])
        if (index+1) % 32 == 0:
            print('Authenticated and hashed actual mapped legs: '+str(index+1)+'/'+str(len(maps)), flush=True)
    for name in ('u', 'p', 'coordinate', 'momentum'):
        identities.append(bounded_array_identity(np.asarray([complex(*z) for z in checkpoint['small'][name]])))
    check('every_actual_seed_array_identity_preserved_under_same_memory_cap', identities == request['initial'])
    with TemporaryDirectory(dir=args.output) as temporary:
        path = Path(temporary)/'writable.npy'
        a = np.arange(4096, dtype=float).reshape(64, 64)
        expected_values = a+1j*(a+7)
        writable = np.lib.format.open_memmap(path, mode='w+', dtype='complex128', shape=a.shape)
        writable[:] = expected_values
        check('dirty_shared_mapping_identity_and_persisted_values_preserved',
              bounded_array_identity(writable) == array_identity(expected_values)
              and mapped_rss(writable) == 0 and np.array_equal(np.load(path), expected_values))
        private = np.load(path, mmap_mode='c')
        private[0, 0] += 11-7j
        expected_private = expected_values.copy(); expected_private[0, 0] += 11-7j
        check('private_copy_on_write_values_not_discarded_or_written_back',
              bounded_array_identity(private) == array_identity(expected_private)
              and np.array_equal(private, expected_private) and np.array_equal(np.load(path), expected_values))
        writable._mmap.close(); private._mmap.close()
    new_identity = producer_identity(producer['model_identity'], producer['q'], producer['dynamics'], current)
    result = dict(schema='P254-joint-partitioned-q-oracle-v1', status='passed',
        extension=retained['extension'], source_bindings=current, source_compatibility=retained['source_compatibility'],
        retained_proof=retained_binding, retained_numerical_checks_not_rerun=retained['retained_numerical_checks_not_rerun'],
        retained_history_ownership_checks_not_rerun=retained['new_check_count'],
        checks=checks, new_check_count=len(checks), proof_sources=[file_binding(Path(__file__).resolve()), *prior_sources],
        original_producer=producer_binding, original_seed=checkpoint_binding, original_propagation_request=request_binding,
        mapped_array_count=len(maps), mapped_payload_bytes=sum(value.nbytes for value in maps),
        first_mapping_rss_before_bytes=before_rss, first_mapping_rss_after_bytes=after_rss,
        original_producer_identity=checkpoint['identity'], current_producer_identity=new_identity,
        model_identity=producer['model_identity'], q=producer['q'], dynamics=producer['dynamics'],
        material_result=False, numerical_accuracy_certified=False,
        scope='All actual seed file bindings and array identities consumed. Only resident page lifetime changes; private-copy mappings remain untouched. Original numerical, field and history proofs retained, not replayed. This is software/runtime recovery evidence, not a propagated q result.')
    for value in maps:
        value._mmap.close()
    return result


def main():
    args = parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=1, max_seconds=args.max_seconds, max_bytes=1048576,
        max_memory_bytes=args.max_memory_bytes, worker_memory_bytes=args.max_memory_bytes)
    budget.reserve(524288)
    started = perf_counter()
    result = dict(schema='P254-q-mapping-lifetime-check-v1', status='failed')
    try:
        result = execute(args, budget)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result.update(seconds=perf_counter()-started, budget=budget.snapshot(),
        process_peak_rss_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        cpu_affinity=sorted(os.sched_getaffinity(0)))
    result['budget']['time_cap'] = 'Cooperative checks and external systemd hard wall cap; standalone reader has no worker watchdog.'
    output = args.output/'result.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output, flush=True)
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()
