"""Small-coordinate Gamma response algebra, in original atom/xyz order.

All stages are independent complex analytic coefficients, ordered
mean, plus, minus, mixed. Matter coordinates are already physical per-cell;
Maxwell mode coordinates are global. No Fourier rescaling is applied here.
Kinetic output is the mixed RESPONSE of the minimally coupled kinetic
component, not an allocation of total energy, heat, covariance or finite dose.
"""
import numpy as np

from mixed_excitation_energy import _bilinear

STAGES = ('mean', 'plus', 'minus', 'mixed')
DIRECTIONS = STAGES[1:]


def vector(value, size, label):
    """Require one finite complex-preserving vector without changing ordering."""
    value = np.asarray(value, dtype=np.complex128)
    if value.shape != (size,) or not np.isfinite(value).all():
        raise ValueError(label+' requires a finite vector of length '+str(size))
    return value


def canonical_stages(mean, small, *, natoms, nfield):
    """Read saved canonical records; never infer minus from plus or zero mixed."""
    result = {'mean': dict(u=vector(mean['u_bohr'], 3*natoms, 'mean.u'),
                          p=vector(mean['p_canonical'], 3*natoms, 'mean.p'),
                          coordinate=vector(mean['field_coordinate'], nfield, 'mean.coordinate'),
                          momentum=vector(mean['field_momentum'], nfield, 'mean.momentum'))}
    if set(small) != set(DIRECTIONS):
        raise ValueError('explicit independent plus/minus/mixed small stages required')
    for direction in DIRECTIONS:
        stage = small[direction]
        result[direction] = {name: vector(getattr(stage, name), size, direction+'.'+name)
                             for name, size in (('u', 3*natoms), ('p', 3*natoms),
                                                ('coordinate', nfield), ('momentum', nfield))}
    return result


def mechanical_stages(stages, *, masses, charges, polarizations):
    """Return k=p-Za and v=k/M; E[mode,xyz] maps GLOBAL modes to shift."""
    if tuple(stages) != STAGES:
        raise ValueError('canonical stage order must be mean,plus,minus,mixed')
    charges = np.asarray(charges, dtype=float)
    masses = np.asarray(masses, dtype=float)
    polarizations = np.asarray(polarizations, dtype=float)
    if (charges.ndim != 1 or masses.shape != (3*len(charges),)
            or polarizations.ndim != 2 or polarizations.shape[1] != 3
            or not np.isfinite(charges).all() or not np.isfinite(masses).all()
            or not np.isfinite(polarizations).all() or np.any(masses <= 0)):
        raise ValueError('original finite Cartesian masses, charges and mode/xyz polarizations required')
    result = {}
    for name in STAGES:
        s = stages[name]
        u = vector(s['u'], len(masses), name+'.u').reshape(-1, 3)
        p = vector(s['p'], len(masses), name+'.p').reshape(-1, 3)
        coordinate = vector(s['coordinate'], len(polarizations), name+'.coordinate')
        momentum = vector(s['momentum'], len(polarizations), name+'.momentum')
        shift = coordinate @ polarizations
        kinetic_momentum = p-charges[:, None]*shift
        result[name] = dict(u=u, p=p, shift=shift, coordinate=coordinate,
                            momentum=momentum, kinetic_momentum=kinetic_momentum,
                            velocity=kinetic_momentum/masses.reshape(-1, 3))
    return result


def kinetic_atomwise(stages, *, masses, ncells):
    """N/2 bilinear(k,k,sum_xyz(a*b/M)) per original atom, without dagger.

    Includes BOTH mean-times-mixed and BOTH independent first cross terms.
    Summing atoms is the existing ionic_kinetic_minimal scalar response.
    """
    if type(ncells) is not int or ncells < 1:
        raise ValueError('positive original cell count required')
    k = tuple(stages[name]['kinetic_momentum'] for name in STAGES)
    m = np.asarray(masses, dtype=float).reshape(-1, 3)
    return ncells*_bilinear(k, k, lambda a, b: np.sum(a*b/m, axis=1))/2


def displacement_weights(u):
    """Diagnostic |u_direction|^2 by atom; NOT a bilinear response or covariance.

    The magnitude diagnostic is explicitly real/nonanalytic, while original
    complex u,p,v and the analytic kinetic response remain unrepaired.
    Undefined fractions for an identically zero direction are reported as None.
    """
    weight = np.sum(np.abs(u)**2, axis=1)
    total = float(np.sum(weight))
    return dict(atom_bohr2=weight, total_bohr2=total,
                atom_fraction=None if total == 0 else weight/total)


def atom_table(stages, kinetic, *, atom_ids, species, positions):
    """37-atom (or finite oracle) rows and additive species/all-atom tables.

    positions has source layout (xyz,atom); only its transpose maps to rows.
    Species u,p,v are vector SUMS, not averages or collective normalized modes.
    """
    atom_ids, species = np.asarray(atom_ids), np.asarray(species)
    positions = np.asarray(positions)
    n = len(atom_ids)
    if positions.shape != (3, n) or species.shape != (n,) or np.asarray(kinetic).shape != (n,):
        raise ValueError('original IDs/species, Cartesian (3,natoms) positions and atom kinetic required')
    weights = {c: displacement_weights(stages[c]['u']) for c in DIRECTIONS}
    rows = []
    for i, atom in enumerate(atom_ids):
        rows.append(dict(atom_original_one_based=int(atom), species=str(species[i]),
            position_bohr_xyz=positions[:, i], kinetic_mixed_global_ry=kinetic[i],
            response={c: dict(u_bohr_xyz=stages[c]['u'][i], p_bohr_inverse_xyz=stages[c]['p'][i],
                velocity_bohr_ry_xyz=stages[c]['velocity'][i],
                displacement_weight_bohr2=weights[c]['atom_bohr2'][i],
                displacement_fraction=None if weights[c]['atom_fraction'] is None
                else weights[c]['atom_fraction'][i]) for c in DIRECTIONS}))
    totals = {}
    for label, mask in [('all', np.ones(n, bool))]+[(s, species == s) for s in dict.fromkeys(species.tolist())]:
        totals[str(label)] = dict(atom_count=int(np.count_nonzero(mask)),
            kinetic_mixed_global_ry=np.sum(kinetic[mask]),
            response={c: dict(u_sum_bohr_xyz=np.sum(stages[c]['u'][mask], axis=0),
                p_sum_bohr_inverse_xyz=np.sum(stages[c]['p'][mask], axis=0),
                velocity_sum_bohr_ry_xyz=np.sum(stages[c]['velocity'][mask], axis=0),
                displacement_weight_bohr2=float(np.sum(weights[c]['atom_bohr2'][mask])),
                displacement_fraction=None if weights[c]['atom_fraction'] is None
                else float(np.sum(weights[c]['atom_fraction'][mask]))) for c in DIRECTIONS})
    return dict(atoms=rows, totals=totals,
                global_field={c: {key: stages[c][key] for key in ('coordinate', 'momentum', 'shift')}
                              for c in DIRECTIONS})
