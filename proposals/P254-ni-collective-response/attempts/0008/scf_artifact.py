"""Independent published SCF extxyz snapshot parser, including stored precision."""
import re
import numpy as np


def parse(text):
    lines=text.splitlines()
    header=dict(re.findall(r'(\w+)=("[^"]*"|\S+)',lines[1]))
    if header['Properties']!='species:S:1:pos:R:3:forces:R:3:magmoms:R:1':
        raise ValueError('Unexpected SCF artifact property layout')
    rows=[line.split() for line in lines[2:]]
    if len(rows)!=int(lines[0]) or any(len(r)!=8 for r in rows):
        raise ValueError('SCF snapshot count/width mismatch')
    return {'symbols':[r[0] for r in rows], 'positions':np.array([r[1:4] for r in rows],float),
            'forces':np.array([r[4:7] for r in rows],float),'magmoms':np.array([r[7] for r in rows],float),
            'cell':np.array(header['Lattice'].strip('"').split(),float).reshape(3,3),
            'stress':np.array(header['stress'].strip('"').split(),float).reshape(3,3),
            'energy':float(header['energy']),'free_energy':float(header['free_energy']), 'pbc':header['pbc']}
