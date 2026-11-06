"""Finite analytic NC model and literal QE records for NEW A27 oracles.

This is not a material surrogate result: execution tests implementation
identities only. No Ni/H file is opened or output imitated.
"""
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path
from types import SimpleNamespace as NS
import struct
import xml.etree.ElementTree as ET
import numpy as np
from substrate_framework.qe_saved_reference import SavedKPoint, SavedReference, SavedSpecies
from substrate_framework.qe_saved_symmetry import parse_saved_symmetries


class AnalyticNC:
    """Defined rank-one s-channel Gaussian separable NC test operator."""
    def __init__(self):
        self.upf = NS(nh=1, angular_l=np.array([0]), channel_radial=np.array([0]),
                      channel_lm=np.array([0]), channel_l=np.array([0]), expanded_d=np.array([[1.3]]))

    def beta(self, radius):
        return np.exp(-.2*np.asarray(radius)**2)[:, None]

    def vloc(self, radius):
        return -.7*np.exp(-.17*np.asarray(radius)**2)


def source_rows(k, cutoff):
    """Literal candidate enumeration, independent of production bounding box."""
    return np.array([m for m in product(range(-3, 4), repeat=3)
                     if sum((float(m[i])+k[i])**2 for i in range(3)) < cutoff], dtype=np.int32).T


def bind(path, data):
    Path(path).write_bytes(data)
    return dict(path=str(path), bytes=str(len(data)), sha256=sha256(data).hexdigest())


def record(data):
    return struct.pack('<i', len(data)) + data + struct.pack('<i', len(data))


def wfc_bytes(ik, spin, k, m, c):
    return (record(struct.pack('<i3diid', ik, *k, spin, 0, 1.))
            + record(struct.pack('<4i', m.shape[1], m.shape[1], 1, c.shape[1]))
            + record(np.eye(3, dtype='<f8').tobytes(order='F'))
            + record(m.astype('<i4').tobytes(order='F'))
            + b''.join(record(c[:, b].astype('<c16').tobytes()) for b in range(c.shape[1])))


def fixture(directory):
    # Identity first, then independently enumerated cubic group. The 3-cycle
    # is deliberately not symmetric, exposing XML C-order interpretation.
    rotations = [np.eye(3, dtype=int)]
    for perm in permutations(range(3)):
        for signs in product((-1, 1), repeat=3):
            s = np.eye(3, dtype=int)[list(perm)]*np.array(signs)[:, None]
            if not any(np.array_equal(s, old) for old in rotations):
                rotations.append(s)
    root = ET.Element('fixture')
    flags = ET.SubElement(ET.SubElement(root, 'input'), 'symmetry_flags')
    ET.SubElement(flags, 'no_t_rev').text = 'false'
    sym = ET.SubElement(ET.SubElement(root, 'output'), 'symmetries')
    for key, value in dict(nsym=48, nrot=48, colin_mag=-1, space_group=0).items():
        ET.SubElement(sym, key).text = str(value)
    for index, s in enumerate(rotations):
        op = ET.SubElement(sym, 'symmetry')
        ET.SubElement(op, 'info', name=str(index)).text = 'crystal_symmetry'
        ET.SubElement(op, 'rotation', dims='3 3').text = ' '.join(map(str, s.ravel(order='F')))
        ET.SubElement(op, 'fractional_translation').text = '0 0 0'
        ET.SubElement(op, 'equivalent_atoms', nat='2').text = '1 2'
    xml = ET.tostring(root)
    xml_binding = bind(Path(directory)/'symmetry.xml', xml)
    saved = parse_saved_symmetries(root, xml_binding, 2, noncolin=False)
    representatives = [(0,0,0), (0,0,.25), (0,0,-.5), (0,.25,.25), (0,.25,-.5),
                       (0,-.5,-.5), (.25,.25,.25), (.25,.25,-.5), (.25,-.5,-.5), (-.5,-.5,-.5)]
    counts = [1,6,3,12,12,3,8,12,6,1]
    points, bindings = [], {}
    for ik, (k, count) in enumerate(zip(representatives, counts), 1):
        k = np.array(k)
        m = source_rows(k, 2.3)
        energies = (np.array([-.4, .8]), np.array([-.1, 1.2]))
        raw = (np.array([1.02, -.02]), np.array([.9, .1]))
        points.append(SavedKPoint(ik, m.shape[1], count/64, k, k, energies, energies, raw))
        for spin in (1, 2):
            rows = m if spin == 1 else m[:, ::-1]
            label = np.arange(rows.shape[1])[:, None]
            c = (np.cos(.31*label+np.array([.3, 1.1])) + 1j*np.sin(.71*label+spin))/np.sqrt(rows.shape[1])
            bindings[spin, ik] = bind(Path(directory)/f'wfc-{spin}-{ik}.dat', wfc_bytes(ik, spin, k, rows, c))
    a = 2*np.pi
    ref = SavedReference(source=xml_binding, creator_version='7.5', source_energy_units='Ha', nat=2,
        species=(SavedSpecies('X', 1., 'analytic'),), atom_names=('X','X'), atom_indices=(1,2),
        positions_bohr=np.array([[0,a/2]]*3), cell_bohr=np.eye(3)*a, alat_bohr=a,
        reciprocal_2pi_alat=np.eye(3), reciprocal_bohr_inv=np.eye(3), fft_grid=(9,9,9),
        fft_smooth=(9,9,9), fft_box=(9,9,9), ngm=125, ngms=125, npwx=100,
        gamma_only=False, lsda=True, noncolin=False, spinorbit=False, nbnd_by_spin=(2,2),
        nks=10, kpoints=tuple(points), nelec=2., occupations_kind='smearing', smearing='mv',
        degauss_ha=.01, degauss_ry=.02, fermi_energies_ha=np.array([0.]),
        fermi_energies_ry=np.array([0.]), total_energies_ha={}, total_energies_ry={},
        ecutwfc_ha=1.15, ecutrho_ha=8., monkhorst_pack=dict(nk1=4,nk2=4,nk3=4,k1=0,k2=0,k3=0),
        saved_symmetries=saved)
    density = NS(miller=np.array(list(product(range(-2,3), repeat=3))).T)
    inputs = NS(reference=ref, density=density, tables={'X': AnalyticNC()}, sources=(dict(xml_binding),))
    return inputs, bindings, root, rotations
