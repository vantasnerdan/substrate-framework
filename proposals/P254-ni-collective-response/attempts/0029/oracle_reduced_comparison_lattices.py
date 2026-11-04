"""Independent exact-value and missing/mismatched-data lattice comparisons."""
from copy import deepcopy
import unittest
import warnings

import numpy as np
import sympy as sp

from reduced_coherence_readout import coherence_curve, normalized_ordered
from reduced_comparison_inputs import CASES
from reduced_comparison_lattices import REPRESENTATIONS, compare_lattices


VALUES = (
    (sp.Rational(1, 5)+sp.I/10, sp.Rational(2, 5), sp.Rational(3, 5), sp.Rational(4, 5)),
    (sp.Rational(9, 10), sp.Rational(91, 100), sp.Rational(23, 25), sp.Rational(93, 100)),
    (sp.Rational(3, 10), sp.Rational(1, 2), sp.Rational(7, 10), sp.Rational(3, 4)),
    (sp.Rational(1, 10), sp.Rational(1, 5), sp.Rational(1, 4), sp.Rational(3, 10)),
)


def covariance(value):
    array = np.zeros(64, complex)
    array[0] = value
    return dict(ordered_ab=array, ordered_ba=array.conj(),
                symmetrized=(array+array.conj())/2, commutator=array-array.conj())


def observations(value):
    rows = []
    for time in (0., 1., 10., 30.):
        variance = 4 if time == 0 else 9
        cross = 4 if time == 0 else 6*complex(value)
        rows.append(dict(time_inverse_ry=time, connected_correlation=dict(
            available=True, reference_time_inverse_ry=0., noise_scope='finite oracle baseline',
            equal_time=covariance(variance), two_time=covariance(cross))))
    return rows


def fixture():
    runs = []
    for case_index, case in enumerate(CASES):
        for index, (bands, features) in enumerate(REPRESENTATIONS):
            runs.append(dict(
                binding=dict(path=f'{case}-{bands}-{features}.json', sha256='0'*64, bytes=1),
                status='complete', specification=dict(case=case, case_identity='own-'+case),
                model_source_identity='original-'+case, coordinate=0,
                configuration=dict(bands=bands, features=features, rtol=2e-8, atol=2e-10),
                seed=dict(atom=1, axis='x', cell=[0, 0, 0], amplitude=1.),
                geometry_arrays=dict(atom_names=np.array(['Ni'])),
                mechanical_arrays=dict(masses=np.array([2., 2., 2.])),
                observations=observations(VALUES[case_index][index])))
    return runs


def evaluate(runs):
    curves = [dict(coherence=coherence_curve(run['observations'])) for run in runs]
    return compare_lattices(runs, curves)


def envelope(report, case, time=30.):
    return next(row for row in report['comparison_groups'][0]['observed_nested_ranges']
                if row['case'] == case and row['time_inverse_ry'] == time)


class LatticeComparisonOracle(unittest.TestCase):
    def test_one_element_variance_arrays_extract_explicit_scalars(self):
        expected = sp.Rational(1, 5)+sp.I/10
        with warnings.catch_warnings():
            warnings.simplefilter('error', DeprecationWarning)
            value = normalized_ordered([6*complex(expected)],
                np.array([[9.+0j]]), np.array([[4.+0j]]))['complex_degree'][0]
        self.assertAlmostEqual(value.real, float(sp.re(expected)))
        self.assertAlmostEqual(value.imag, float(sp.im(expected)))

    def test_exact_complex_range_and_conditional_ordering(self):
        with warnings.catch_warnings():
            warnings.simplefilter('error', DeprecationWarning)
            report = evaluate(fixture())
        self.assertTrue(report['numerical_comparison_complete'])
        self.assertIsNone(report['physical_lattice_ranking'])
        self.assertFalse(report['material_accuracy_certified'])
        row = envelope(report, CASES[0])
        low = sp.Abs(VALUES[0][0])
        self.assertAlmostEqual(row['minimum_magnitude'], float(low))
        self.assertAlmostEqual(row['maximum_magnitude'], float(sp.Rational(4, 5)))
        self.assertAlmostEqual(row['observed_spread'], float(sp.Rational(4, 5)-low))
        contrasts = report['comparison_groups'][0]['descriptive_pairwise_contrasts']
        separated = next(row for row in contrasts if row['left_case'] == CASES[0]
                         and row['right_case'] == CASES[1] and row['time_inverse_ry'] == 30.)
        self.assertEqual(separated['conditional_model_relation'],
                         'right_greater')
        overlapping = next(row for row in contrasts if row['left_case'] == CASES[0]
                           and row['right_case'] == CASES[2] and row['time_inverse_ry'] == 30.)
        self.assertEqual(overlapping['conditional_model_relation'],
                         'overlap_or_touch')

    def test_extreme_representation_mutation_changes_the_range(self):
        runs = fixture()
        before = envelope(evaluate(runs), CASES[0])
        runs[3]['observations'] = observations(sp.Rational(1, 2))
        after = envelope(evaluate(runs), CASES[0])
        expected = sp.Rational(3, 5)-sp.Abs(VALUES[0][0])
        self.assertAlmostEqual(after['observed_spread'], float(expected))
        self.assertAlmostEqual(before['observed_spread']-after['observed_spread'],
                               float(sp.Rational(4, 5)-sp.Rational(3, 5)))

    def test_missing_resolution_is_not_filled_with_another_basis(self):
        runs = fixture()
        del runs[3]
        report = evaluate(runs)
        self.assertFalse(report['numerical_comparison_complete'])
        row = envelope(report, CASES[0])
        self.assertFalse(row['complete_observed_nested_set'])
        self.assertNotIn('minimum_magnitude', row)
        self.assertEqual(row['missing'], [dict(bands=12, features=27, state='absent')])

    def test_nearby_time_does_not_replace_required_thirty(self):
        runs = fixture()
        runs[3]['observations'][-1]['time_inverse_ry'] = 29.999
        report = evaluate(runs)
        self.assertFalse(report['numerical_comparison_complete'])
        self.assertFalse(envelope(report, CASES[0])['complete_observed_nested_set'])
        self.assertTrue(envelope(report, CASES[0], 10.)['complete_observed_nested_set'])

    def test_dose_species_mass_and_tolerances_do_not_mix(self):
        for change in ('dose', 'species', 'mass', 'tolerance'):
            with self.subTest(change=change):
                runs = fixture()
                if change == 'dose':
                    runs[3]['seed']['amplitude'] = 2.
                elif change == 'species':
                    runs[3]['geometry_arrays']['atom_names'][0] = 'H'
                elif change == 'mass':
                    runs[3]['mechanical_arrays']['masses'][0] = 3.
                else:
                    runs[3]['configuration']['rtol'] = 1e-7
                report = evaluate(runs)
                self.assertFalse(report['numerical_comparison_complete'])
                self.assertTrue(all(not group['numerical_comparison_complete']
                                    for group in report['comparison_groups']))

    def test_repeated_resolution_is_ambiguous_not_averaged(self):
        runs = fixture()
        repeat = deepcopy(runs[3])
        repeat['binding']['path'] = 'distinct-repeat.json'
        repeat['observations'] = observations(sp.Rational(1, 2))
        report = evaluate(runs+[repeat])
        self.assertFalse(report['numerical_comparison_complete'])
        row = envelope(report, CASES[0])
        self.assertEqual(row['missing'], [dict(bands=12, features=27, state='ambiguous')])
        self.assertNotIn('observed_spread', row)

    def test_zero_actual_reference_variance_cannot_complete_a_panel(self):
        runs = fixture()
        runs[3]['observations'][0]['connected_correlation']['equal_time'] = covariance(0.)
        report = evaluate(runs)
        self.assertFalse(report['numerical_comparison_complete'])
        self.assertFalse(envelope(report, CASES[0])['complete_observed_nested_set'])


if __name__ == '__main__':
    unittest.main()
