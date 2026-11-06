"""Consume actual A22/Saved QE inputs for Gamma, finite q and raw remainder.

No Native/cloud work. No WFCs. No screened modes or repaired electronic blocks.
The optional complete grid calls the existing force-lattice Fourier API.
Construction worker has not executed this consumer.
"""
import argparse
from dataclasses import fields
from pathlib import Path
import hashlib
import json

import numpy as np

from substrate_framework.material_frozen_lattice import load_material_frozen_lattice
from substrate_framework.material_force_lattice import (
    q_grid_to_force_lattice, force_lattice_to_q_grid, lattice_residuals)
from substrate_framework.qe_saved_native import read_bound_bytes
from material_run_contract import runtime_context, implementation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sources', type=Path, default=Path(__file__).with_name('ionic-q-sources.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--complete-grid', action='store_true')
    args = parser.parse_args()
    # New output only; this is not an overwrite/resume or allocation system.
    args.output.mkdir(parents=True, exist_ok=False)
    context = runtime_context()
    with (args.output/'runtime-context.json').open('x') as stream:
        json.dump(context, stream, indent=2)
    manifest_raw=args.sources.read_bytes()
    manifest=json.loads(manifest_raw)
    for source in manifest['producer_sources'].values():
        read_bound_bytes(source)
    data=load_material_frozen_lattice(args.sources)
    frozen, ref = data.lattice, data.reference
    qfrac=np.array([0.,0.,.25])
    q=ref.reciprocal_bohr_inv @ qfrac
    gamma_terms=frozen.ionic.terms(np.zeros(3))
    qterms=frozen.ionic.terms(q)
    minus=frozen.ionic.matrix(-q)
    ionic_q=qterms.total.reshape(3*ref.nat,3*ref.nat)
    frozen_q=frozen.saved_gamma+(ionic_q-frozen.ionic_gamma)
    remainder=frozen.remainder_report()
    arrays=dict(saved_gamma=frozen.saved_gamma, ionic_gamma=frozen.ionic_gamma,
                ionic_q=ionic_q, ionic_minus_q=minus, frozen_q=frozen_q,
                electronic_remainder=frozen.electronic_remainder,
                electronic_onsite_blocks=remainder.pop('onsite_blocks_ry_per_bohr2'),
                atom_masses_amu=data.atom_masses_amu, q_fractional=qfrac,
                q_cart_bohr_inv=q, valence_charges=frozen.ionic.charges)
    for prefix, terms in (('gamma',gamma_terms),('q',qterms)):
        for field in fields(terms):
            arrays[prefix+'_'+field.name]=getattr(terms,field.name)
    # These are inspection outputs, not a demand that frozen electronic curvature
    # obey an ASR on its own or that the imported decimal residual vanish.
    report=dict(scope='actual material saved-reference frozen ionic-q extension, not screened response',
        sources=manifest, source_manifest=dict(path=str(args.sources.resolve()),bytes=str(len(manifest_raw)),
            sha256=hashlib.sha256(manifest_raw).hexdigest()),
        atom_names=list(ref.atom_names), atom_masses_amu=data.atom_masses_amu.tolist(),
        ecutrho_ry=2*ref.ecutrho_ha, signed_G_count=ref.ngm,
        producer_alpha_printed=data.producer_alpha, alpha_binary_reconstructed=data.computed_alpha,
        real_cutoff_bohr=frozen.ionic.real_cutoff,
        q_fractional=qfrac.tolist(), q_cart_bohr_inv=q.tolist(),
        remainder=remainder,
        ionic_gamma_hermiticity=float(np.linalg.norm(frozen.ionic_gamma-frozen.ionic_gamma.conj().T)),
        ionic_q_hermiticity=float(np.linalg.norm(ionic_q-ionic_q.conj().T)),
        ionic_q_opposite_conjugacy=float(np.linalg.norm(minus-ionic_q.conj())),
        ionic_gamma_translation=float(np.linalg.norm(gamma_terms.total.sum(axis=2))),
        gamma_saved_recovery_max=float(np.max(np.abs(frozen.matrix(np.zeros(3))-frozen.saved_gamma))),
        ionic_q_minus_gamma_norm=float(np.linalg.norm(ionic_q-frozen.ionic_gamma)),
        term_norms={prefix:{field.name:float(np.linalg.norm(getattr(terms,field.name))) for field in fields(terms)}
                    for prefix,terms in (('gamma',gamma_terms),('q',qterms))},
        limitations=manifest['limitations'])
    report['runtime_context'] = context
    report['implementation_sources'] = implementation((Path(__file__),))
    if args.complete_grid:
        mesh=ref.monkhorst_pack
        if any(mesh['k'+str(i)] != 0 for i in (1,2,3)):
            raise ValueError('unshifted admitted BvK grid required')
        shape=tuple(mesh['nk'+str(i)] for i in (1,2,3))
        grid=frozen.complete_q_grid(shape)
        lattice=q_grid_to_force_lattice(grid)
        arrays.update(frozen_q_grid=grid, frozen_force_lattice=lattice)
        opposite=grid
        for axis,n in enumerate(shape):
            opposite=np.take(opposite,(-np.arange(n))%n,axis=axis)
        report['complete_grid']=dict(shape=shape, slots=int(np.prod(shape)),
            force_lattice=lattice_residuals(lattice),
            folded_opposite_conjugacy=float(np.linalg.norm(opposite-grid.conj())),
            fourier_roundtrip_max=float(np.max(np.abs(force_lattice_to_q_grid(lattice)-grid))),
            scope='finite periodic frozen force lattice; finite G boundary residuals retained')
    archive=args.output/'ionic-q-material.npz'
    with archive.open('xb') as stream:
        np.savez_compressed(stream,**arrays)
    raw=archive.read_bytes()
    report['archive']=dict(path=str(archive.resolve()),bytes=str(len(raw)),sha256=hashlib.sha256(raw).hexdigest())
    with (args.output/'result.json').open('x',encoding='utf-8') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.output),remainder=remainder,
        gamma_saved_recovery_max=report['gamma_saved_recovery_max'],
        ionic_q_minus_gamma_norm=report['ionic_q_minus_gamma_norm']),indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
