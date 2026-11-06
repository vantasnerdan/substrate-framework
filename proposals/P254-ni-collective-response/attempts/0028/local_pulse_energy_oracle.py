"""Independent differential-operator and exact finite-DFT pulse oracle.

No energy-consumer imports and no expected numerical dose fixtures.  The
operator acts on an arbitrary wavefunction, not an assumed zero mean.
Numerical comparisons only connect the supplied Fourier implementation to
an exact roots-of-unity DFT; all physical energy identities are symbolic.
"""
from itertools import product
from math import isfinite

import sympy as sp


def _exact(value):
    return sp.simplify(sp.expand(value))


def _identity(actual, expected):
    defect = _exact(actual-expected)
    return dict(passed=bool(defect == 0), residual=str(defect))


def _operator_checks(energy_function):
    x, a, y = sp.symbols('x a y', real=True)
    j, z, pbar, abar = sp.symbols('J Z pbar abar', real=True, nonzero=True)
    hbar, mass, field_mass, electron_mass = sp.symbols('hbar M Mfield me', positive=True)
    psi = sp.Function('psi')(x, a, y)
    phase = sp.exp(sp.I*j*x/hbar)

    def mechanical(wave):
        return -sp.I*hbar*sp.diff(wave, x)-z*a*wave

    def kinetic(wave):
        return mechanical(mechanical(wave))/(2*mass)

    shifted = _exact(kinetic(phase*psi)/phase-kinetic(psi))
    expected_operator = j*mechanical(psi)/mass+j**2*psi/(2*mass)
    momentum_shift = _exact((-sp.I*hbar*sp.diff(phase*psi, x))/phase
                            -(-sp.I*hbar*sp.diff(psi, x)))
    potential = sp.Function('V')(x, a, y)
    spectators = dict(
        coordinate_potential=lambda wave: potential*wave,
        field_kinetic=lambda wave: -hbar**2*sp.diff(wave, a, 2)/(2*field_mass),
        electronic_kinetic=lambda wave: -hbar**2*sp.diff(wave, y, 2)/(2*electron_mass))
    checks = dict(
        positive_canonical_shift=_identity(momentum_shift, j*psi),
        minimal_coupling_operator_difference=_identity(shifted, expected_operator))
    for name, operator in spectators.items():
        checks[name+'_commutes'] = _identity(operator(phase*psi)/phase, operator(psi))
    # Reading the coefficients of psi and d_x psi yields the expectation for
    # arbitrary first moments.  No pbar=abar=0 specialization licenses this test.
    derived_mean = j*(pbar-z*abar)/mass+j**2/(2*mass)
    actual = energy_function(mass, z, j, pbar, abar)
    expected = dict(mechanical_mean=pbar-z*abar,
                    linear_cross=j*(pbar-z*abar)/mass,
                    quadratic=j**2/(2*mass), total=derived_mean)
    for name, value in expected.items():
        checks['arbitrary_mean_'+name] = _identity(actual[name], value)
    checks['canonical_mean_sensitivity'] = _identity(sp.diff(actual['total'], pbar), j/mass)
    checks['minimal_coupling_sensitivity'] = _identity(sp.diff(actual['total'], abar), -j*z/mass)
    centered = energy_function(mass, z, j, 0, 0)
    checks['centered_specialization'] = _identity(centered['total'], j**2/(2*mass))
    checks['zero_charge_does_not_remove_p_cross'] = _identity(
        energy_function(mass, 0, j, pbar, abar)['total'], j*pbar/mass+j**2/(2*mass))
    # A different amplitude for a different mass gives the same centered dose.
    other_mass = sp.Symbol('M_other', positive=True)
    matched_j = j*sp.sqrt(other_mass/mass)
    checks['matched_dose_mass_rescaling'] = _identity(
        energy_function(other_mass, z, matched_j, 0, 0)['total'], centered['total'])
    baseline = all(c['passed'] for c in checks.values())
    candidates = dict(
        reversed_pulse_sign=energy_function(mass, z, -j, pbar, abar)['total'],
        wrong_minimal_coupling_sign=energy_function(mass, -z, j, pbar, abar)['total'],
        dropped_mechanical_cross=actual['quadratic'],
        dropped_field_cross=j*pbar/mass+actual['quadratic'],
        doubled_quadratic_kinetic=actual['linear_cross']+2*actual['quadratic'],
        missing_Ry_mass_half=energy_function(2*mass, z, j, 0, 0)['total'])
    mutations = {}
    for name, candidate in candidates.items():
        target = centered['total'] if name == 'missing_Ry_mass_half' else derived_mean
        defect = _exact(candidate-target)
        variables = dict(reversed_pulse_sign=(j, pbar), wrong_minimal_coupling_sign=(j, abar),
            dropped_mechanical_cross=(j, pbar), dropped_field_cross=(j, abar),
            doubled_quadratic_kinetic=(j, j), missing_Ry_mass_half=(j, j))[name]
        witness = _exact(sp.diff(defect, *variables))
        mutations[name] = dict(detected=bool(baseline and defect != 0 and witness != 0),
            exact_defect=str(defect), sensitivity=dict(variables=list(map(str, variables)),
                                                       derivative=str(witness)))
    return checks, mutations, dict(
        unitary=str(phase), canonical_convention='[x,p]=i*hbar; active state U psi',
        transformed_kinetic_difference=str(shifted),
        expectation_difference=str(derived_mean),
        hypotheses=['finite kinetic energy and first moments',
                    'canonical pulse independent of field/electronic coordinates',
                    'all non-ionic-momentum Hamiltonian terms commute with this pulse'],
        implication='Quantum fluctuations need not vanish; their unchanged covariance cancels in the initial energy difference.')


def _fourier_checks(fourier_function, mesh, cell):
    """Exact 4^3 roots, non-origin phase, inverse, and all64 Parseval."""
    if tuple(mesh) != (4, 4, 4):
        raise ValueError('this source-bound oracle is for the original 4^3 mesh')
    indices = list(product(range(4), repeat=3))
    n = len(indices)
    root_n = sp.sqrt(n)
    j = sp.Symbol('J', real=True, nonzero=True)
    checks, mutations, evidence = {}, {}, {}
    for label, position in (('registered_cell', tuple(cell)), ('non_origin', (1, 2, 3))):
        # Direct forward DFT of a delta, from exact fourth roots; no production
        # phase/weight intermediate is used by this derivation.
        phases = [sp.I**(-sum(k*r for k, r in zip(step, position))) for step in indices]
        exact_p = [j*phase/root_n for phase in phases]
        norm = sp.expand(sum(sp.conjugate(v)*v for v in exact_p))
        checks[label+'_Parseval_all64'] = _identity(norm, j**2)
        recovered = []
        for position_out in indices:
            value = sp.expand(sum(sp.I**sum(k*r for k, r in zip(step, position_out))*v
                                  for step, v in zip(indices, exact_p))/root_n)
            recovered.append(value)
        checks[label+'_inverse_all64'] = dict(
            passed=all(_exact(value-(j if r == position else 0)) == 0
                       for r, value in zip(indices, recovered)),
            reconstruction=[str(v) for v in recovered])
        produced = fourier_function(list(mesh), list(position), 1.0)
        structure = (len(produced) == n and [tuple(row['mesh_index']) for row in produced] == indices)
        errors, weight_errors, q_errors = [], [], []
        for row, step, phase in zip(produced, indices, phases):
            expected = complex(phase/root_n)
            errors.append(abs(row['coefficient']-expected))
            weight_errors.append(abs(row['quadratic_weight']-1/n))
            q_errors.extend(abs(q-k/4) for q, k in zip(row['q_fractional'], step))
        max_error = max(errors+weight_errors+q_errors, default=float('inf'))
        tolerance = 128*float(sp.Float(2)**-52)
        checks[label+'_implementation_against_exact_DFT'] = dict(
            passed=bool(structure and isfinite(max_error) and max_error <= tolerance),
            coefficient_max_absolute_error=max(errors, default=float('inf')),
            weight_max_absolute_error=max(weight_errors, default=float('inf')),
            q_max_absolute_error=max(q_errors, default=float('inf')),
            tolerance=tolerance, criterion='absolute error for unit impulse; dose itself is never a fixed numeric fixture')
        evidence[label] = dict(cell=list(position), exact_coefficients_per_J=[str(v/j) for v in exact_p],
                               q_count=n, exact_total_weight_per_J_squared=str(norm/j**2))
        if label == 'non_origin':
            # These phase-only defects preserve Parseval; the phase oracle,
            # rather than a norm tally, is what exposes them.
            phase_defects = dict(reversed_Fourier_phase=[sp.conjugate(v) for v in exact_p],
                                 omitted_non_origin_phase=[j/root_n for _ in indices])
            for name, values in phase_defects.items():
                witnesses = [(list(indices[k]), str(_exact(v-exact_p[k])))
                             for k, v in enumerate(values) if _exact(v-exact_p[k]) != 0]
                mutations[name] = dict(detected=bool(witnesses), phase_witnesses=witnesses,
                    Parseval_still_passes=bool(_exact(sum(sp.conjugate(v)*v for v in values)-j**2) == 0))
    # Independent symbolic energy-counting defects.  One Gamma share is not
    # all64; conjugate edges are not an additional set of physical modes.
    mass = sp.Symbol('M', positive=True)
    whole = j**2/(2*mass)
    for name, candidate in dict(Gamma_as_whole_local=whole/n,
                                missing_sqrtN=whole*n,
                                double_count_q_minus_q=2*whole).items():
        defect = _exact(candidate-whole)
        mutations[name] = dict(detected=bool(defect != 0), exact_defect=str(defect))
    baseline = all(c['passed'] for c in checks.values())
    for item in mutations.values():
        item['detected'] = bool(baseline and item['detected'])
    return checks, mutations, evidence


def check_pulse(energy_function, fourier_function, *, mesh, cell, mass_unit):
    """JSON-safe evidence; call only in the reviewed real CLI run."""
    checks, mutations, derivation = _operator_checks(energy_function)
    fchecks, fmutations, fourier = _fourier_checks(fourier_function, mesh, cell)
    checks.update(fchecks)
    mutations.update(fmutations)
    # Starting with p^2/(2 m) Hartree and E_Ry=2 E_Ha gives p^2/m_e_ratio
    # Ry, hence the mass in p^2/(2 M_Ry) is m_e_ratio/2.  Runtime constants
    # come from the authenticated source, not a second decimal fixture bank.
    amu_si = sp.Rational(str(mass_unit['AMU_SI']))
    me_si = sp.Rational(str(mass_unit['ELECTRONMASS_SI']))
    expected = amu_si/(2*me_si)
    actual = float(mass_unit['AMU_RY'])
    error = abs(actual-float(expected))
    bound = 4*float(sp.Float(2)**-52)*abs(float(expected))
    checks['source_Ry_mass_convention'] = dict(passed=bool(isfinite(error) and error <= bound),
        exact_SI_ratio_over_two=str(expected), actual_source_AMU_RY=actual,
        absolute_error=error, floating_roundoff_comparison_bound=bound,
        derivation='E_Ry=2*E_Ha; p^2/(2*m_in_electron_masses) Ha = p^2/(2*M_Ry) Ry implies M_Ry=m_in_electron_masses/2')
    passed = all(c['passed'] for c in checks.values()) and all(m['detected'] for m in mutations.values())
    return dict(status='passed' if passed else 'failed',
        arithmetic='exact SymPy differential operators and fourth-root DFT; explicit floating bridge to implementation',
        checks=checks, mutations=mutations, derivation=derivation, Fourier=fourier,
        scope='Initial input identity and centered Fourier bookkeeping only; no propagated energy/heat/noise oracle.')
