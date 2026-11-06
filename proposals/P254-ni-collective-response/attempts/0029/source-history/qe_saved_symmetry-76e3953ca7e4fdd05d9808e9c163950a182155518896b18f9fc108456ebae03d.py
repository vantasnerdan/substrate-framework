"""Lossless QE XML symmetry metadata; effective scalar policy is separate.

qexsd_copy:277 reshapes the serialized rotation in Fortran order. Omitted
info/time_reversal is not an XML false: qexsd_init:291 and symm_base's zero
default justify the separate scalar interpretation only for colin_mag=-1.
"""
from dataclasses import dataclass
from types import MappingProxyType
import numpy as np


@dataclass(frozen=True)
class SavedSymmetry:
    index: int
    rotation: np.ndarray
    fractional_translation: np.ndarray | None
    equivalent_atoms: np.ndarray | None
    info: str
    raw_info_attributes: object
    raw_rotation_attributes: object
    raw_time_reversal: bool | None
    effective_time_reversal: bool | None
    effective_policy: str


@dataclass(frozen=True)
class SavedSymmetries:
    nsym: int
    nrot: int
    colin_mag: int | None
    space_group: int
    operations: tuple[SavedSymmetry, ...]
    raw_flags: object
    source: object


def _numbers(node, size, integer=False):
    values = np.array((node.text or '').split(), dtype=float)
    if values.size != size or not np.isfinite(values).all():
        raise ValueError('invalid saved symmetry array')
    if integer:
        if not np.equal(values, np.rint(values)).all():
            raise ValueError('noninteger saved rotation/atom permutation')
        values = values.astype(np.int64)
    values.setflags(write=False)
    return values


def parse_saved_symmetries(root, source, nat, *, noncolin):
    """Consume the existing XML tree, preserving all nrot operations in order."""
    node = root.find('output/symmetries')
    if node is None:
        return None
    nsym, nrot = (int(node.findtext(key)) for key in ('nsym', 'nrot'))
    raw_colin = node.findtext('colin_mag')
    colin = int(raw_colin) if raw_colin is not None else None
    operations = []
    for index, operation in enumerate(node.findall('symmetry'), 1):
        info, rotation = operation.find('info'), operation.find('rotation')
        if info is None or rotation is None:
            raise ValueError('missing saved symmetry info/rotation')
        if rotation.attrib.get('order', 'F') != 'F':
            raise ValueError('QE symmetry requires Fortran matrix order')
        if tuple(map(int, rotation.attrib.get('dims', '3 3').split())) != (3, 3):
            raise ValueError('saved rotation must be 3 by 3')
        s = _numbers(rotation, 9, True).reshape((3, 3), order='F')
        ft, atoms = operation.find('fractional_translation'), operation.find('equivalent_atoms')
        f = _numbers(ft, 3) if ft is not None else None
        irt = _numbers(atoms, nat, True) if atoms is not None else None
        if irt is not None and (int(atoms.attrib.get('nat', nat)) != nat or
                               sorted(irt.tolist()) != list(range(1, nat + 1))):
            raise ValueError('saved equivalent_atoms is not a full atom permutation')
        raw = info.attrib.get('time_reversal')
        if raw not in (None, 'true', 'false'):
            raise ValueError('invalid saved symmetry time_reversal')
        raw_bool = None if raw is None else raw == 'true'
        default = not noncolin and colin in (None, -1)
        effective = raw_bool if raw is not None else (False if default else None)
        policy = 'explicit XML attribute' if raw is not None else (
            'QE7.5 scalar t_rev=0; qexsd_init:291/symm_base:63' if default else 'unresolved')
        operations.append(SavedSymmetry(index, s, f, irt, (info.text or '').strip(),
                          MappingProxyType(dict(info.attrib)), MappingProxyType(dict(rotation.attrib)),
                          raw_bool, effective, policy))
    if not 1 <= nsym <= nrot or len(operations) != nrot:
        raise ValueError('saved symmetry counts disagree')
    flags = root.find('input/symmetry_flags')
    raw_flags = {} if flags is None else {child.tag: (child.text or '').strip() for child in flags}
    return SavedSymmetries(nsym, nrot, colin, int(node.findtext('space_group')),
                           tuple(operations), MappingProxyType(raw_flags), source)
