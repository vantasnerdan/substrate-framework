"""Final periodic distances, native moments and pressure without seed-site labels."""
from itertools import combinations
from pathlib import Path
import re
import numpy as np
import sympy as sp
from substrate_framework.material_primary_boundary import cubic_pair_distance, image_enumeration_distance


def distance_block(ledger,c):
    positions=np.array(c['final']['positions_angstrom']); cell=np.array(c['final']['cell_angstrom']); length=cell[0,0]
    pairs=list(combinations(range(31,37),2))+[(i,j) for i in range(31,37) for j in range(31)]
    result=[]
    for i,j in pairs:
        exact=cubic_pair_distance(positions[i],positions[j],length)
        independent=image_enumeration_distance(positions[i],positions[j],cell)
        value=float(exact)
        ledger.check(f'periodic pair {i+1}/{j+1}: exact wrap vs 27-image SciPy norm',abs(value-independent)<2e-14)
        result.append({'i_1based':i+1,'j_1based':j+1,'pair':'HH' if j>=31 else 'HNi','distance_angstrom':value,'exact_decimal_distance':str(exact)})
    # Boundary-crossing witness exposes a missing minimum-image map, not just symmetry.
    witness=next(r for r in result if abs(np.linalg.norm(positions[r['i_1based']-1]-positions[r['j_1based']-1])-r['distance_angstrom'])>1.)
    i,j=witness['i_1based']-1,witness['j_1based']-1
    ledger.mutation_sensitive('periodic distance is not unwrapped or doubled',lambda d:abs(d-image_enumeration_distance(positions[i],positions[j],cell))<2e-14,witness['distance_angstrom'],[np.linalg.norm(positions[i]-positions[j]),2*witness['distance_angstrom']])
    summaries=[]
    for i in range(32,38):
        neighbours=sorted((r for r in result if r['i_1based']==i and r['pair']=='HNi'),key=lambda r:r['distance_angstrom'])
        summaries.append({'H_index_1based':i,'ordered_Ni_neighbours':neighbours,
            'nearest_shell_count_within_1e-6_angstrom_of_min':sum(r['distance_angstrom']<=neighbours[0]['distance_angstrom']+1e-6 for r in neighbours),
            'shell_tolerance_angstrom':1e-6})
    return {'all_HH_pairs':[r for r in result if r['pair']=='HH'],'H_Ni_by_H':summaries,
            'scope':'Minimum periodic pair distances from exact saved final coordinates; shell counts use explicitly stated descriptive 1e-6 A equality tolerance, not an assumed bond criterion.'}


def native_block(ledger,c):
    out=c['native']; u=c['units']; text=Path(c['scf']['data']['artifacts']['pw.stdout']['path']).read_text()
    lines=text.splitlines()
    independent={}
    for label in ('total magnetization','absolute magnetization'):
        independent[label]=float([line for line in lines if line.strip().startswith(label)][-1].split('=')[1].split()[0])
    ledger.mutation_sensitive('native final total magnetic moment extraction',lambda m:m==independent['total magnetization'],out['moment'],[0.,out['moment']+1.])
    ledger.mutation_sensitive('native final absolute magnetic moment extraction',lambda m:m==independent['absolute magnetization'],out['absolute_moment'],[0.,out['moment']])
    local=np.array([float(x) for x in re.findall(r'atom\s+\d+ \(R=[^\n]+magn=\s*([-\d.]+)',text)])
    ledger.check('all site sphere-integrated moments copied at source precision',local.shape==(37,) and np.max(abs(local-c['snapshot']['magmoms']))<=1e-12)
    trace=float(np.trace(out['stress_kbar'])/3)
    ledger.mutation_sensitive('native pressure equals printed stress trace not zero equilibrium',lambda p:abs(p-trace)<=0.01,out['pressure_kbar'],[0.,-out['pressure_kbar']])
    # SI pressure route checks both stdout unit columns; rounding in Ry is propagated.
    factor=u['rydberg_ev']*u['elementary_charge_si']/(u['bohr_angstrom']*1e-10)**3/1e8
    stress_allowance=0.5e-8*factor+0.005
    ledger.check('native Ry stress independently converts to kbar columns',np.max(abs(out['stress_ry']*factor-out['stress_kbar']))<=stress_allowance)
    # ASE stress is negative QE stress and expressed in eV/A^3.
    stored=-out['stress_ry']*u['rydberg_ev']/u['bohr_angstrom']**3
    ledger.mutation_sensitive('SCF stress snapshot convention',lambda s:np.max(abs(s-stored))<1e-13,c['snapshot']['stress'],[-c['snapshot']['stress'],np.zeros((3,3))])
    return {'total_moment_bohr_magneton_per_cell':out['moment'],'absolute_moment_bohr_magneton_per_cell':out['absolute_moment'],
            'moment_print_half_quantum':0.005,'site_sphere_moments_bohr_magneton':local.tolist(),
            'site_moment_print_half_quantum':0.00005,'pressure_kbar':out['pressure_kbar'],
            'pressure_print_half_quantum_kbar':0.005,'stress_ry_per_bohr3':out['stress_ry'].tolist(),
            'stress_kbar':out['stress_kbar'].tolist(),'stress_column_rounding_allowance_kbar':stress_allowance,
            'scope':'Native moment on one converged spin initialization and clamped-cell pressure, not unique spin ground state or pressure equilibrium.'}


def run(ledger,c):
    return {'periodic_distances':distance_block(ledger,c),'native_magnetism_pressure':native_block(ledger,c)}
