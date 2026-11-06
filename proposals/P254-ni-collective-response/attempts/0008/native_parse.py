"""Independent text extraction; never calls the production ASE/QE parser."""
import re
import shlex
import numpy as np

NUMBER = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'


def one(pattern, text):
    matches = re.findall(pattern, text, re.M)
    if len(matches) != 1:
        raise ValueError(f'Expected singleton {pattern!r}, found {len(matches)}')
    return matches[0]


def numeric(token):
    return float(token.replace('D', 'E').replace('d', 'e'))


def card(lines, header, count):
    start = lines.index(header) + 1
    return [shlex.split(line) for line in lines[start:start+count]]


def parse_input(text):
    lines = text.splitlines()
    positions = card(lines, 'ATOMIC_POSITIONS angstrom', 37)
    species = card(lines, 'ATOMIC_SPECIES', 2)
    values = dict(re.findall(r'^\s*([\w_]+(?:\(\d+\))?)\s*=\s*([^,\n]+)', text, re.M))
    return {'values': {k: v.strip().strip("'") for k,v in values.items()},
            'symbols': [r[0] for r in positions],
            'positions': np.array([[numeric(x) for x in r[1:4]] for r in positions]),
            'position_widths': [len(r) for r in positions],
            'cell': np.array(card(lines, 'CELL_PARAMETERS angstrom', 3), float),
            'species': {r[0]: {'mass': numeric(r[1]), 'pseudo': r[2]} for r in species},
            'kpoints': list(map(int, card(lines, 'K_POINTS automatic', 1)[0]))}


def parse_extxyz(text):
    lines = text.splitlines()
    n = int(lines[0])
    header = dict(re.findall(r'(\w+)=("[^"]*"|\S+)', lines[1]))
    fields = header['Properties'].split(':')
    if fields != ['species','S','1','pos','R','3','masses','R','1']:
        raise ValueError('Unexpected physical geometry properties or constraint column')
    rows = [line.split() for line in lines[2:]]
    if len(rows) != n or any(len(r) != 5 for r in rows):
        raise ValueError('Frame/atom count mismatch')
    return {'symbols': [r[0] for r in rows],
            'positions_angstrom': [[numeric(x) for x in r[1:4]] for r in rows],
            'masses_amu': [numeric(r[4]) for r in rows],
            'cell_angstrom': np.array(header['Lattice'].strip('"').split(),float).reshape(3,3).tolist(),
            'pbc': [x == 'T' for x in header['pbc'].strip('"').split()], 'constraint_count': 0}


def parse_stdout(text):
    # Stop BEFORE verbose force decomposition: the last atom-force block is NOT total force.
    section = text.split('Forces acting on atoms (cartesian axes, Ry/au):')
    if len(section) != 2:
        raise ValueError('Expected exactly one total-force section')
    force_rows = re.findall(r'atom\s+(\d+) type\s+(\d+)\s+force =\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')', section[1].split('The non-local contrib.')[0])
    geom = text.split('positions (alat units)')[1].split('Crystallographic axes')[0]
    atoms = re.findall(r'(\d+)\s+(Ni|H)\s+tau\(\s*(\d+)\) = \(\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')', geom)
    axes = re.findall(r'a\(\d\) = \(\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')', text)
    stress_part = text.split('total   stress  (Ry/bohr**3)')[1]
    stress_rows = [line.split() for line in stress_part.splitlines()[1:4]]
    return {'forces': np.array([[numeric(x) for x in r[2:]] for r in force_rows]),
            'force_ids': [int(r[0]) for r in force_rows], 'types': [int(r[1]) for r in force_rows],
            'symbols': [r[1] for r in atoms], 'atom_ids': [(int(r[0]),int(r[2])) for r in atoms],
            'fractional': np.array([[numeric(x) for x in r[3:]] for r in atoms]),
            'axes': np.array(axes,float), 'alat_bohr': numeric(one(r'celldm\(1\)=\s*('+NUMBER+r')',text)),
            'energy_ry': numeric(one(r'^!\s+total energy\s*=\s*('+NUMBER+r') Ry',text)),
            'iterations': int(one(r'convergence has been achieved in\s+(\d+) iterations',text)),
            'version': one(r'Program PWSCF v\.([^\s]+)',text),
            'done_count': len(re.findall(r'JOB DONE\.',text)),
            'pressure_kbar': numeric(one(r'P=\s*('+NUMBER+r')',stress_part)),
            'stress_ry': np.array([r[:3] for r in stress_rows],float),
            'stress_kbar': np.array([r[3:] for r in stress_rows],float),
            'moment': numeric(re.findall(r'total magnetization\s*=\s*('+NUMBER+r')',text)[-1]),
            'absolute_moment': numeric(re.findall(r'absolute magnetization\s*=\s*('+NUMBER+r')',text)[-1])}
