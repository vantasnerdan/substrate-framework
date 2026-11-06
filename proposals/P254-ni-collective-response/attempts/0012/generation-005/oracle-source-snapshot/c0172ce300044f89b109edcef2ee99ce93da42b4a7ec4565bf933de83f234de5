"""Read native input and total-force decimal tokens independently of ASE."""
import re
import shlex
import numpy as np

NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?'
FORCE=re.compile(r'^\s*atom\s+(\d+)\s+type\s+(\d+)\s+force\s*=\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s*$',re.M)


def singleton(pattern,text):
    matches=re.findall(pattern,text,re.M)
    if len(matches)!=1:
        raise ValueError(f'Expected singleton {pattern}, found {len(matches)}')
    return matches[0]


def number(token):
    return float(token.replace('D','E').replace('d','e'))


def input_text(text):
    lines=text.splitlines()
    pairs=re.findall(r'^\s*([\w_]+(?:\(\d+\))?)\s*=\s*([^,\n]+)',text,re.M)
    if len(dict(pairs))!=len(pairs):
        raise ValueError('Duplicate native setting')
    values={k:v.strip().strip("'") for k,v in pairs}
    nat,ntyp=int(values['nat']),int(values['ntyp'])
    def card(header,count):
        if lines.count(header)!=1:
            raise ValueError('Missing/duplicate native card '+header)
        start=lines.index(header)+1
        return [shlex.split(line) for line in lines[start:start+count]]
    positions=card('ATOMIC_POSITIONS angstrom',nat)
    species=card('ATOMIC_SPECIES',ntyp)
    if any(len(row)!=4 for row in positions) or any(len(row)!=3 for row in species):
        raise ValueError('Incomplete or constrained native geometry')
    return {'values':values,'symbols':[r[0] for r in positions],
            'position_tokens':[r[1:] for r in positions],
            'positions':np.array([[number(x) for x in r[1:]] for r in positions]),
            'cell':np.array(card('CELL_PARAMETERS angstrom',3),float),
            'species':species,'kpoints':list(map(int,card('K_POINTS automatic',1)[0]))}


def stdout_text(text,nat,symbols):
    sections=text.split('Forces acting on atoms (cartesian axes, Ry/au):')
    if len(sections)!=2 or 'The non-local contrib.' not in sections[1]:
        raise ValueError('Missing/duplicate total-force boundary')
    total=sections[1].split('The non-local contrib.',1)[0]
    rows=FORCE.findall(total)
    if len(rows)!=nat or len(re.findall(r'\batom\b',total))!=nat:
        raise ValueError('Incomplete/malformed total-force table')
    ids=[int(r[0]) for r in rows]; types=[int(r[1]) for r in rows]
    expected=[1 if symbol=='Ni' else 2 for symbol in symbols]
    if ids!=list(range(1,nat+1)) or types!=expected:
        raise ValueError('Force atom order/type mismatch')
    geom=text.split('positions (alat units)')
    if len(geom)!=2:
        raise ValueError('Ambiguous stdout geometry')
    atoms=re.findall(r'(\d+)\s+(Ni|H)\s+tau\(\s*(\d+)\)\s*=\s*\(\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')',geom[1].split('Crystallographic axes')[0])
    if [(int(r[0]),r[1],int(r[2])) for r in atoms]!=[(i,s,i) for i,s in enumerate(symbols,1)]:
        raise ValueError('Stdout atom identity/order mismatch')
    axes=re.findall(r'a\(\d\) = \(\s*('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')',text)
    if len(axes)!=3 or text.count('JOB DONE.')!=1 or 'convergence NOT achieved' in text:
        raise ValueError('Missing complete converged native stdout')
    return {'force_tokens':[list(r[2:]) for r in rows], 'force_ids':ids,'force_types':types,
            'energy_ry':singleton(r'^!\s+total energy\s*=\s*('+NUMBER+r')\s+Ry',text),
            'total_moment':re.findall(r'total magnetization\s*=\s*('+NUMBER+r')',text)[-1],
            'absolute_moment':re.findall(r'absolute magnetization\s*=\s*('+NUMBER+r')',text)[-1],
            'iterations':int(singleton(r'convergence has been achieved in\s+(\d+) iterations',text)),
            'scf_accuracy_ry':re.findall(r'estimated scf accuracy\s*<\s*('+NUMBER+r')\s+Ry',text)[-1],
            'alat_bohr':singleton(r'celldm\(1\)=\s*('+NUMBER+r')',text),
            'axes':np.array(axes,float),
            'fractional':np.array([[number(x) for x in r[3:]] for r in atoms]),
            'version':singleton(r'Program PWSCF v\.([^\s]+)',text),'job_done':True}
