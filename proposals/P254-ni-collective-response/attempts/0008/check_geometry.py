"""Bind ordered physical geometry, mobile host, isotope convention and model."""
import copy
import re
from pathlib import Path
import numpy as np
from native_parse import parse_input, parse_stdout, parse_extxyz, one, NUMBER


def run(ledger, c):
    p=c['published']['data']; s=c['scf']['data']; a=Path(s['attempt_directory'])
    inp=parse_input((a/'pw.in').read_text()); out=parse_stdout((a/'pw.stdout').read_text())
    c['input'],c['native']=inp,out
    final=parse_extxyz(Path(p['artifacts']['final.extxyz']['path']).read_text())
    initial=parse_extxyz(Path(c['relax']['structure']['path']).read_text())
    ledger.mutation_sensitive('physical geometry full publication and constraint identity',lambda x:x==final,p['geometry'],
        [{**p['geometry'],'constraint_count':1},{**p['geometry'],'symbols':['H']+p['geometry']['symbols'][1:]},
         {**p['geometry'],'masses_amu':[58.6934]*37},
         {**p['geometry'],'positions_angstrom':[[0.,0.,0.]]+p['geometry']['positions_angstrom'][1:]}])
    seed=parse_extxyz(Path(c['spec']['structure']['path']).read_text())
    ledger.check('force seed exactly physical final file',seed==final)
    ledger.check('37 ordered atoms 31 Ni then six H in input stdout and publication',
        all(v==['Ni']*31+['H']*6 for v in (final['symbols'],initial['symbols'],inp['symbols'],out['symbols'],s['symbols'])))
    ledger.check('stdout atom and force identifiers retain ordering',out['force_ids']==list(range(1,38))
        and out['atom_ids']==[(i,i) for i in range(1,38)] and out['types']==[1]*31+[2]*6)
    ledger.check('all 111 host/H coordinate components mobile in SCF input',inp['position_widths']==[4]*37 and final['constraint_count']==0 and final['pbc']==[True]*3)
    ledger.check('input cell and masses retain initial clamped model',final['cell_angstrom']==initial['cell_angstrom']==inp['cell'].tolist()
        and final['masses_amu']==initial['masses_amu'] and final['cell_angstrom']==(np.eye(3)*7.048).tolist())
    ledger.check('input coordinates exactly physical publication',np.array_equal(inp['positions'],final['positions_angstrom']))
    expected_mass=[58.6934]*31+[1.00782503223]*6
    ledger.mutation_sensitive('natural Ni mean and neutral atomic H1 masses',lambda x:x==expected_mass,final['masses_amu'],
        [[58.6934]*31+[2.01410177812]*6,[58.6934]*31+[1.007276466621]*6,[58.0]*31+[1.00782503223]*6])
    isotope=Path(c['supplement']['isotope_source']).read_text()
    ledger.check('H1 and D atomic masses sourced, not nuclear masses',all(str(x) in isotope for x in (1.00782503223,2.01410177812)))
    ledger.check('native species card masses and pseudos agree',all(inp['species'][sym]['mass']==expected_mass[i] and inp['species'][sym]['pseudo']==c['spec']['execution']['pseudopotentials'][sym] for sym,i in [('Ni',0),('H',31)]))
    values=inp['values']; settings=c['spec']['settings']
    expected={'ecutwfc':'40.0','ecutrho':'480.0','degauss':'0.02','conv_thr':'1e-10','electron_maxstep':'200',
              'nspin':'2','starting_magnetization(1)':'0.06','starting_magnetization(2)':'0.0',
              'calculation':'scf','ibrav':'0','nat':'37','ntyp':'2','input_dft':'PBE','smearing':'mv','occupations':'smearing'}
    ledger.mutation_sensitive('native fixed-cell primary model settings',lambda x: all(x.get(k)==v for k,v in expected.items()),values,[{**values,'ecutwfc':'50.0'},{**values,'nspin':'1'},{**values,'calculation':'vc-relax'}])
    ledger.check('native automatic Gamma-centered k444 grid',inp['kpoints']==[4,4,4,0,0,0] and settings['kpts']==[4]*3)
    text=(a/'pw.stdout').read_text()
    for sym,i in [('Ni',0),('H',31)]:
        mass=float(one(r'^\s*'+sym+r'\s+(?:18\.00|1\.00)\s+('+NUMBER+r')\s+',text))
        ledger.check(sym+' stdout rounded mass agrees with exact card',abs(mass-expected_mass[i])<=0.5e-5)
    bohr=c['units']['bohr_angstrom']; cell=out['axes']*out['alat_bohr']*bohr
    positions=out['fractional']*out['alat_bohr']*bohr
    # celldm(1) 6 decimals; tau 7 decimals; axes 6 decimals. Cubic axes are exact by the bound input.
    quantum_cell=0.5e-6*bohr
    quantum_pos=0.5e-7*out['alat_bohr']*bohr+np.max(abs(out['fractional']))*quantum_cell+1e-14
    ledger.check('stdout geometry publication reproduced at floating precision',np.max(abs(cell-np.array(s['cell_angstrom'])))<1e-13 and np.max(abs(positions-np.array(s['positions_angstrom'])))<1e-13)
    ledger.mutation_sensitive('native rounded geometry matches exact force geometry',lambda x:np.max(abs(x-inp['positions']))<=quantum_pos,positions,[positions+np.eye(37,3)*1e-3])
    ledger.check('stdout fixed cell agrees within source print precision',np.max(abs(cell-inp['cell']))<=quantum_cell)
    from scf_artifact import parse
    snapshot=parse(Path(s['artifacts']['final.extxyz']['path']).read_text())
    c['snapshot']=snapshot
    ledger.check('published SCF snapshot geometry matches stdout at eight decimals',snapshot['symbols']==out['symbols'] and np.max(abs(snapshot['positions']-positions))<=0.5e-8+1e-14 and np.max(abs(snapshot['cell']-cell))<1e-13 and snapshot['pbc']=='"T T T"')
    c['final']=final
    return {'atom_count':len(final['symbols']),'Ni_count':final['symbols'].count('Ni'),'H_count':final['symbols'].count('H'),
            'physical_constraint_count':final['constraint_count'],'mobile_coordinates':int(inp['positions'].size),
            'masses_amu':final['masses_amu'],'cell_angstrom':final['cell_angstrom'],
            'native_stdout_cell_angstrom':cell.tolist(),'native_position_max_discrepancy_angstrom':float(np.max(abs(positions-inp['positions']))),
            'source_print_position_bound_angstrom':float(quantum_pos),'source_print_cell_bound_angstrom':quantum_cell,
            'mass_model':'uniform natural-Ni mean, neutral atomic H-1; no isotope disorder, nuclear mass or exchange model',
            'native_input_values':values,'kpoints':inp['kpoints']}
