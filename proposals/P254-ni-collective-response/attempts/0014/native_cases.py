"""Prepare isolated same-functional SCF cases using the existing QE stage API."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import threading
import time
import xml.etree.ElementTree as ET
from decimal import Decimal

import numpy as np
from ase.io import read, write
from ase.io.formats import get_ioformat
from ase.units import create_units
from catdft import qe_stage
from substrate_framework.periodic_curvature import conjugate_forces

HERE = Path(__file__).resolve().parent
REG = json.loads((HERE / 'registration.json').read_text())
ROOT = Path('/home/administrator/projects/catalyst-nc-operator')
WORK = ROOT / 'runs/2026-10-19_collective-response/coherence-action-0014'
PYTHON = '/home/administrator/projects/catalyst-dft/.venv/bin/python'
BOHR = create_units('2006')['Bohr']


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def setup():
    # Preload the actual selected lazy output parser before source inventory.
    get_ioformat('espresso-out')
    template = json.loads(Path(REG['inputs']['template_spec']).read_text())
    primitive = read(REG['inputs']['primitive_pw_in'], format='espresso-in')
    primitive.calc = None
    paths = set(template['input_files_sha256']) | {str(p) for p in qe_stage._source_paths()}
    bindings = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in sorted(paths)}
    save(WORK / 'current-source-bindings.json', {
        'bindings': bindings,
        'changed_from_historical_template': [p for p, h in template['input_files_sha256'].items()
                                             if bindings[p] != h],
        'scope': 'Current source epoch; no retrospective expansion of G10 admission.'})
    return template, primitive, bindings


def prepare(template, primitive, bindings, q, name, amplitudes, seed=None):
    copies = 1 if q == 'gamma' else 2
    atoms = primitive.repeat((copies, 1, 1))
    signs = [1] if copies == 1 else [1, -1]
    for r, sign in enumerate(signs):
        atoms.positions[2*r, 0] += sign * amplitudes[0] * BOHR
        atoms.positions[2*r+1, 1] += sign * amplitudes[1] * BOHR
    directory = WORK / q / name
    geometry = directory / 'structure.extxyz'
    directory.mkdir(parents=True, exist_ok=True)
    write(geometry, atoms, format='extxyz')
    spec = {k: template[k] for k in ('schema_version', 'execution')}
    spec.update(id='p254-action14-' + q + '-' + name, directory=str(directory),
                structure={'path': str(geometry)}, settings=dict(template['settings']))
    if copies == 2:
        spec['settings'].update(kpts=[1, 2, 2], nbnd=40)
    if seed:
        spec['initial_state'] = (dict(seed) if isinstance(seed, dict) else
                                {'kind': 'qe75-collinear-dat-save-v1', 'directory': str(seed)})
    inputs = dict(bindings)
    inputs[str(geometry)] = hashlib.sha256(geometry.read_bytes()).hexdigest()
    spec['input_files_sha256'] = inputs
    if seed:
        from catdft import qe_scf_seed
        for path in qe_scf_seed.input_files(spec, atoms):
            inputs[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
    path = directory / 'case.json'
    save(path, spec)
    return path, signs


def _sample_native(pid, stop, records):
    while not stop.wait(.02):
        queue = [pid]
        while queue:
            parent = queue.pop()
            try:
                children = Path(f'/proc/{parent}/task/{parent}/children').read_text().split()
                queue.extend(int(p) for p in children)
                if Path(f'/proc/{parent}/comm').read_text().strip() != 'pw.x':
                    continue
                tasks = {p.name: sorted(os.sched_getaffinity(int(p.name)))
                         for p in Path(f'/proc/{parent}/task').iterdir()}
                records[str(parent)] = tasks
            except (FileNotFoundError, ProcessLookupError):
                continue


def run_case(path, signs):
    argv = [PYTHON, '-u', '-m', 'catdft.qe_stage', '--spec', str(path)]
    directory = path.parent
    invocation = directory / 'invocations' / f'{time.time_ns()}-{os.getpid()}'
    invocation.mkdir(parents=True, exist_ok=False)
    save(invocation / 'invoke.json', {'argv': argv, 'started_ns': str(time.time_ns()),
                                     'timeout_seconds': 90})
    started = time.monotonic()
    process = subprocess.Popen(argv, env=os.environ.copy(), stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stop, affinity = threading.Event(), {}
    sampler = threading.Thread(target=_sample_native, args=(process.pid, stop, affinity))
    sampler.start()
    try:
        stdout, stderr = process.communicate(timeout=90)
    except subprocess.TimeoutExpired as exc:
        (invocation / 'wrapper.stdout').write_bytes(exc.stdout or b'')
        (invocation / 'wrapper.stderr').write_bytes(exc.stderr or b'')
        save(invocation / 'execution.json', {'status': 'case_timeout', 'native_affinity': affinity})
        raise SystemExit('Case budget reached; owner exit lets systemd stop its whole cgroup.')
    finally:
        stop.set()
        sampler.join()
    (invocation / 'wrapper.stdout').write_bytes(stdout)
    (invocation / 'wrapper.stderr').write_bytes(stderr)
    save(invocation / 'execution.json', {'returncode': process.returncode,
        'wall_seconds': time.monotonic()-started, 'native_affinity': affinity})
    print(stdout.decode(), end='', flush=True)
    if process.returncode:
        raise RuntimeError(stderr.decode() + f' SCF failed: {path}')
    if affinity and any(mask != [6] for tasks in affinity.values() for mask in tasks.values()):
        raise RuntimeError('Observed native task outside the declared singleton CPU6')
    result = json.loads((directory / 'result.json').read_text())
    xml_path = Path(result['data']['artifacts']['data-file-schema.xml']['path'])
    xml = ET.fromstring(xml_path.read_bytes())
    energy = Decimal(xml.findtext('output/total_energy/etot').strip()) * 2 / len(signs)
    forces = np.array([float(v) for v in xml.findtext('output/forces').split()]).reshape(-1, 3) * 2
    record = {'id': result['job_id'], 'energy_ry_per_primitive': str(energy),
        'conjugate_force_ry_per_bohr': conjugate_forces(forces, signs).tolist(),
        'magnetization_per_primitive': float(xml.findtext('output/magnetization/total'))/len(signs),
        'xml': str(xml_path), 'xml_sha256': hashlib.sha256(xml_path.read_bytes()).hexdigest(),
        'result': str(directory / 'result.json'), 'native_affinity_observed': bool(affinity),
        'native_timing_seconds': result['data']['native_xml']['native_timing_seconds']}
    save(directory / 'action-observables.json', record)
    print('ACTION14_CASE_COMPLETE ' + result['job_id'], flush=True)
    return record
