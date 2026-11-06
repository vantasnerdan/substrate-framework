"""NEW finite oracle source; run only with Root's source-bound execution grant.

PYTHONPATH=src python proposals/P254-ni-collective-response/attempts/0027/run_zero_temperature_oracle.py

Fixture numbers are dimensionless, not Ni/H data. Comparator enumerates Fock
configurations and grand Gibbs weights with its own common-mu solve, then
compares the zero-T selector. Selected CAR moments test the NEW preparation
seam only, inheriting rather than replaying the accepted factor-algebra proof.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np

from substrate_framework.material_preparation import (
    FullGridSector, MaterialPopulationSpecification, MaterialPreparationPolicy,
    SourceBinding, instantiate_material_sector,
)
from substrate_framework.physical_factor_preparation import (
    ActionFactors, FactorFunctional, FactorPreparation,
    coefficients_from_orthonormal_populations,
)
from substrate_framework.zero_temperature_preparation import (
    FillingPolicy, ModeLevel, fill_zero_temperature,
)

TOL = 3e-11  # small matrices/finite enumeration; not a material tolerance
POLICY = FillingPolicy('analytic fixture energies', 'two particles / two cells',
                       'explicit independent replicas / two cells',
                       'exact supplied energy equality only', 'fixture energy', 'particles/cell')


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def close(actual, expected, label):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape, label + ': incompatible shapes')
    require(np.isfinite(actual).all() and np.isfinite(expected).all(), label + ': nonfinite operands')
    error = float(np.max(np.abs(actual-expected), initial=0))
    require(error < TOL, label + ': ' + str(error))
    return error


def occupations(n):
    return np.array([[(state >> j) & 1 for j in range(n)] for state in range(1 << n)])


def gibbs_fock(energies, beta, target):
    """Direct many-body diagonal Gibbs ensemble; no selector/FD formula reused."""
    configurations = occupations(len(energies))
    particle_count = configurations.sum(axis=1)
    energy = configurations @ energies

    def ensemble(mu):
        log_weight = -beta * (energy - mu * particle_count)
        weight = np.exp(log_weight - log_weight.max())
        return weight / weight.sum()

    left, right = -8.0, 8.0
    for _ in range(160):
        middle = (left + right) / 2
        if ensemble(middle) @ particle_count < target:
            left = middle
        else:
            right = middle
    mu = (left + right) / 2
    return configurations, ensemble(mu), mu


def bilinear(matrix):
    """Literal signed creation/annihilation on bit-labelled Fock configurations."""
    n = len(matrix)
    result = np.zeros((1 << n, 1 << n), complex)
    for state in range(1 << n):
        for j in range(n):
            if not (state >> j) & 1:
                continue
            removed = state ^ (1 << j)
            sign_j = (-1) ** (state & ((1 << j) - 1)).bit_count()
            for i in range(n):
                if (removed >> i) & 1:
                    continue
                sign_i = (-1) ** (removed & ((1 << i) - 1)).bit_count()
                result[removed | (1 << i), state] += sign_i * sign_j * matrix[i, j]
    return result


def functional(psi, matrices, identity):
    probes = {name: ActionFactors(matrix @ psi, matrix.conj().T @ psi,
              name, name + '^dagger', 'fixture units', 'fixture full PW', 'analytic Fock oracle')
              for name, matrix in matrices.items()}
    return FactorFunctional.from_actions(psi, probes, factor_identity=identity,
                                         pw_basis_identity='fixture full PW')


def filling_oracle():
    multiplicities = np.array([1, 2, 1, 2])
    energies = np.array([-2.0, 0.0, 0.0, 3.0])
    spins = (1, 1, 2, 2)
    levels = tuple(ModeLevel(str(i), spins[i], energies[i], multiplicities[i] / 2)
                   for i in range(4))
    filling = fill_zero_temperature(levels, 1.0, policy=POLICY)
    expected = np.array([1.0, 1 / 3, 1 / 3, 0.0])
    residuals = {'populations': close(filling.populations, expected, 'fractional scalar block'),
                 'neutral': close(filling.weighted_number, 1, 'weighted neutral number'),
                 'unequal_spin': close([n for _, n in filling.spin_numbers], [5 / 6, 1 / 6], 'free M')}
    expanded = np.repeat(energies, multiplicities)
    errors = []
    for beta in (8, 16, 32):
        configurations, weights, mu = gibbs_fock(expanded, beta, 2)
        marginal = weights @ configurations
        errors.append(float(np.max(np.abs(marginal - np.repeat(expected, multiplicities)))))
    require(errors[-1] < TOL and errors[0] > errors[-1], 'independent Gibbs limit convergence')
    residuals['gibbs_errors_beta_8_16_32'] = errors
    residuals['scaled_multiplier'] = close(32 * mu, np.log(0.5), 'common boundary multiplier slope')
    gapped = fill_zero_temperature(levels, 0.5, policy=POLICY)
    close(gapped.populations, [1, 0, 0, 0], 'gapped filling')
    require(gapped.chemical_potential_interval == (-2.0, 0.0), 'gapped multiplier interval')
    gap_configs, gap_weights, _ = gibbs_fock(expanded, 32, 1)
    close(gap_weights @ gap_configs, [1, 0, 0, 0, 0, 0], 'gapped Gibbs comparator')
    near = (ModeLevel('a', 1, 0., 1.), ModeLevel('b', 2, 1e-10, 1.))
    close(fill_zero_temperature(near, .4, policy=POLICY).populations, [.4, 0], 'no near clustering')
    for target, expected_edge in ((0., np.zeros(4)), (3., np.ones(4))):
        close(fill_zero_temperature(levels, target, policy=POLICY).populations, expected_edge, 'capacity edge')
    count = configurations.sum(axis=1)
    variance = weights @ (count ** 2) - (weights @ count) ** 2
    residuals['independent_collective_noise'] = close(
        variance / 2, sum(level.measure * p * (1 - p) for level, p in zip(levels, expected)),
        'N/sqrt(two cells) connected variance uses linear multiplicity')
    # Exact canonical N=2 ground ensemble: one low mode and one of three boundary modes.
    energy = configurations @ expanded
    selected = (count == 2) & (energy == energy[count == 2].min())
    canonical = selected / selected.sum()
    close(canonical @ count ** 2 - (canonical @ count) ** 2, 0, 'exact-N variance')
    require(variance > 0.6, 'fixed-mean N must not be relabelled exact N')
    close(canonical @ configurations, weights @ configurations, 'same one-body density, different ensemble')
    canonical_cross = canonical @ (configurations[:, 1] * configurations[:, 2])
    grand_cross = weights @ (configurations[:, 1] * configurations[:, 2])
    require(abs(canonical_cross - grand_cross) > .1, 'population pair correlations distinguish ensembles')
    mutations = {
        'spin_two': abs(2 * filling.weighted_number - 1),
        'squared_number_weights': abs(sum(l.measure ** 2 * p for l, p in zip(levels, expected)) - 1),
        'squared_fluctuation_weights': abs(sum(l.measure ** 2 * p * (1-p) for l, p in zip(levels, expected)) - variance / 2),
        'forced_half_fill': abs(sum(l.measure * p for l, p in zip(levels, [1, .5, .5, 0])) - 1),
    }
    require(all(value > 1e-3 for value in mutations.values()), 'normalization/fill mutations must be exposed')
    return filling, residuals, mutations


def frame_oracle(filling):
    q = np.eye(5, dtype=complex)[:, :4]
    root = np.array([[1.4, .12j, .08, 0], [-.12j, 1.1, .07, .03j],
                     [.08, .07, .9, .04], [0, -.03j, .04, 1.2]], complex)
    psi = q @ root  # root is explicitly positive Hermitian, not a unitary frame
    number = np.eye(5)
    probe = np.zeros((5, 5), complex)
    probe[1, 2], probe[2, 1] = 1j, -1j  # exact degenerate off-diagonal population noise
    probe[0, 4] = probe[4, 0] = .4       # initially empty, accessible full-PW complement
    matrices = {'N': number, 'A': probe}
    f = functional(psi, matrices, 'original raw oracle frame')
    d = coefficients_from_orthonormal_populations(f.gram, filling.populations,
                                                 mapping_authority='declared symmetric fixture frame')
    prepared = FactorPreparation(f, d, 'conditional test-only independent quasifree state')
    population = np.array([1., 1 / 3, 1 / 3, 0., 0.])
    configurations = occupations(5)
    probability = np.prod(np.where(configurations, population, 1-population), axis=1)
    # The comparator is a density on Fock space, not another Wick contraction.
    fock_n, fock_a = bilinear(number), bilinear(probe)
    mean = lambda op: probability @ np.diag(op)
    residuals = {'raw_density': close(psi @ d @ psi.conj().T, np.diag(population), 'Gram transport'),
                 'population_noise': close(prepared.ordered('N', 'N'), mean(fock_n @ fock_n)-mean(fock_n)**2, 'Fock number noise'),
                 'degenerate_complement_noise': close(prepared.ordered('A', 'A'), mean(fock_a @ fock_a)-mean(fock_a)**2, 'Fock off-diagonal and complement noise')}
    generator = np.diag([-3., 2., 1., 4., 6.]).astype(complex)
    generator[0, 1] = generator[1, 0] = .3
    rho = psi @ d @ psi.conj().T
    drift = np.linalg.norm(generator @ rho - rho @ generator)
    require(drift > .1, 'noncommuting generator must retain initial drift, not equilibrium')
    gauge = np.eye(4, dtype=complex)
    gauge[1:3, 1:3] = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
    rotated = psi @ gauge
    dg = coefficients_from_orthonormal_populations(rotated.conj().T @ rotated,
            filling.populations, mapping_authority='exact degenerate unitary gauge')
    residuals['block_gauge'] = close(rotated @ dg @ rotated.conj().T, rho, 'block gauge covariance')
    wrong_gram = psi @ np.diag(filling.populations) @ psi.conj().T
    split = np.diag([1., 0., 1., 0.])  # still weighted-neutral, but not block scalar
    mutations = {'wrong_Gram': float(np.linalg.norm(wrong_gram-rho)),
                 'band_label_split': float(np.linalg.norm(q @ (gauge @ split @ gauge.conj().T-split) @ q.conj().T)),
                 'false_H0_equilibrium': float(drift)}
    require(all(value > .1 for value in mutations.values()), 'frame/gauge/equilibrium mutations exposed')
    # Exercise raw instantiation with explicit TEST-ONLY source identities, not material data.
    raw = b'analytic five-row Fock fixture; not an actual WFC'
    binding = SourceBinding('oracle:analytic-frame', len(raw), hashlib.sha256(raw).hexdigest())
    sector = FullGridSector(1, (0, 0, 0), 1, binding, f.factor_identity,
                           f.pw_basis_identity, 5, 'identity fixture transport')
    policy = MaterialPreparationPolicy(*('oracle authority ' + str(i) for i in range(8)))
    other_sector = replace(sector, spin=2, frame_identity='other raw oracle frame')
    full_levels = tuple(ModeLevel('complete fixture ' + str(i), 1 if i < 4 else 2, e, 1.)
                        for i, e in enumerate([-2., 0., 0., 3., -3., 0., 3., 3.]))
    full_filling = fill_zero_temperature(full_levels, 3., policy=replace(POLICY,
        number_authority='three particles in one fixture cell',
        measure_authority='one independent mode per explicit spin/band at single k'))
    specification = MaterialPopulationSpecification(binding, (('fixture', binding),),
                    (('fixture', 3.),), (1, 1, 1), (sector, other_sector),
                    ((0, 4), (4, 8)), full_filling, policy)
    kwargs = dict(source_binding=binding, reference_binding=binding,
                  coverage_authority=policy.coverage_authority, frame_authority='literal analytic raw fixture')
    instance = instantiate_material_sector(specification, sector.key, f, psi, **kwargs)
    close(instance.density_coefficients, d, 'material raw instantiation seam only')
    try:
        instantiate_material_sector(specification, sector.key, replace(f, gram=np.eye(4)), psi, **kwargs)
    except ValueError:
        mutations['wrong_functional_Gram_rejected'] = True
    else:
        raise AssertionError('raw frame/functional Gram mismatch accepted')
    return residuals, mutations


def main():
    filling, residuals, mutations = filling_oracle()
    frame_residuals, frame_mutations = frame_oracle(filling)
    files = (Path(__file__), Path('src/substrate_framework/zero_temperature_preparation.py'),
             Path('src/substrate_framework/material_preparation.py'))
    bindings = []
    for path in files:
        data = path.read_bytes()
        bindings.append({'path': str(path), 'bytes': str(len(data)), 'sha256': hashlib.sha256(data).hexdigest()})
    print(json.dumps({'scope': 'NEW small conditional preparation seam; no actual material result',
                      'residuals': {**residuals, **frame_residuals},
                      'mutation_exposures': {**mutations, **frame_mutations}, 'source_bindings': bindings}, indent=2))


if __name__ == '__main__':
    main()
