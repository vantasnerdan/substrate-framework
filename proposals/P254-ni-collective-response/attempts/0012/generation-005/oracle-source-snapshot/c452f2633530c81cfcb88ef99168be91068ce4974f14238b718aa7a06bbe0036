"""Validate native model identity separately from scientific force decisions."""
import re
from pathlib import Path
import numpy as np
from native_text import number
from source_graph import digest,fingerprint


def frame(text):
    lines=text.splitlines(); n=int(lines[0])
    header=dict(re.findall(r'(\w+)=("[^"]*"|\S+)',lines[1]))
    fields=header['Properties'].split(':'); offsets={}; width=0
    for i in range(0,len(fields),3):
        offsets[fields[i]]=(width,int(fields[i+2])); width+=int(fields[i+2])
    rows=[line.split() for line in lines[2:]]
    if len(rows)!=n or any(len(row)!=width for row in rows):
        raise ValueError('Malformed single extxyz frame')
    def column(name):
        start,count=offsets[name]
        return [row[start:start+count] for row in rows]
    return {'symbols':[x[0] for x in column('species')],
            'positions':np.array(column('pos'),float),
            'cell':np.array(header['Lattice'].strip('"').split(),float).reshape(3,3),
            'masses':np.array(column('masses'),float).ravel() if 'masses' in offsets else None,
            'forces':np.array(column('forces'),float) if 'forces' in offsets else None,
            'pbc':header['pbc'].strip('"').split(),'properties':set(offsets)}


def envelope_valid(result,spec):
    declared=spec['input_files_sha256']; published=result['data']['input_files_sha256']
    # Frozen qe_io.verify_input_files publishes resolved filesystem identities.
    source_links=(set(published)=={str(Path(p).resolve()) for p in declared} and
                  all(published.get(str(Path(p).resolve()))==sha for p,sha in declared.items()))
    return (result['schema_version']==1 and result['status']=='complete' and
            result['job_id']==spec['id'] and result['input_sha256']==fingerprint(spec) and
            result['data']['returncode']==0 and result['data']['job_done'] is True and
            result['data']['scf_converged'] is True and
            source_links)


def artifacts_valid(artifacts,spec):
    expected={'command.json','final.extxyz','pw.in','pw.stdout','pw.stderr','spec.json'}
    return (set(artifacts)==expected and all(Path(row['path']).is_relative_to(spec['directory'])
            and Path(row['path']).name==name and digest(row['path'])==row['sha256']
            for name,row in artifacts.items()))


def input_matches(native,geometry,settings,execution):
    values=native['values']
    numeric={'ecutwfc':'ecutwfc_ry','ecutrho':'ecutrho_ry','degauss':'degauss_ry',
             'conv_thr':'conv_thr_ry','electron_maxstep':'electron_maxstep','nspin':'nspin',
             'starting_magnetization(1)':'starting_magnetization_ni'}
    fixed={'calculation':'scf','restart_mode':'from_scratch','ibrav':'0','nat':'37','ntyp':'2',
           'input_dft':'PBE','occupations':'smearing','smearing':'mv','starting_magnetization(2)':'0.0',
           'tprnfor':'.true.','tstress':'.true.','verbosity':'high'}
    return (all(values.get(k)==v for k,v in fixed.items()) and
            set(values)==set(fixed)|set(numeric)|{'prefix','pseudo_dir','outdir'} and
            all(number(values[k])==settings[v] for k,v in numeric.items()) and
            values['pseudo_dir']==execution['pseudo_dir'] and
            native['kpoints']==settings['kpts']+[0,0,0] and
            native['symbols']==geometry['symbols'] and
            np.array_equal(native['positions'],geometry['positions']) and
            np.array_equal(native['cell'],geometry['cell']) and
            [row[0] for row in native['species']]==['Ni','H'] and
            all(row[2]==execution['pseudopotentials'][row[0]] for row in native['species']) and
            np.array_equal(np.array([number(native['species'][0 if s=='Ni' else 1][1])
                            for s in native['symbols']]),geometry['masses']))


def command_valid(command,result,spec):
    execution=spec['execution']; argv=command['argv']
    expected=['/usr/bin/mpirun','-np',str(execution['mpi_ranks']),'--bind-to','core',
              execution['pw_executable'],'-nk',str(execution['kpoint_pools']),'-in','pw.in']
    return (argv==expected and argv==result['data']['command'] and
            command['cwd']==result['data']['attempt_directory'] and
            command['environment']==result['data']['environment'] and
            all(command['environment'][k]=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')) and
            command['mpi_launcher_sha256']==digest(Path(argv[0]).resolve()))
