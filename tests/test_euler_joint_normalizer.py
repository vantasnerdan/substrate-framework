"""Check the exact full forms, including complex Bloch diagonal phase rows."""

import pytest
import sympy as sp

from substrate_framework.euler_joint_normalizer import common_null_gram_lift


def test_actual_gram_with_complex_rows_and_scale():
    a, b, c, d, e, f, g, h = sp.symbols("a b c d e f g h", real=True)
    scale = sp.symbols("scale", positive=True)
    energy = sp.Matrix([[a, b + sp.I*c], [b - sp.I*c, d]])
    phase = sp.Matrix([[sp.I*e, f + sp.I*g], [-f + sp.I*g, sp.I*h]])
    columns = common_null_gram_lift(energy, phase, scale)
    # Build the physical input form conventions independently by basis entries.
    h_form, o_form = sp.zeros(8), sp.zeros(8)
    for index in range(2):
        h_form[index, index] = 1
        h_form[2+index, 2+index] = -1
        o_form[4+index, 6+index] = 1
        o_form[6+index, 4+index] = -1
    assert sp.simplify(columns.conjugate().T*h_form*columns) == energy
    assert sp.simplify(columns.conjugate().T*o_form*columns) == phase
    wrong_signature = h_form.copy()
    wrong_signature[2, 2] = 1
    assert sp.simplify(columns.conjugate().T*wrong_signature*columns-energy) != sp.zeros(2)
    assert isinstance(columns, sp.ImmutableMatrix)


def test_zero_targets_retain_a_nonzero_common_null_frame():
    columns = common_null_gram_lift([[0]], [[0]])
    assert columns == sp.Matrix([1, 1, 1, 0])


@pytest.mark.parametrize(
    "energy,phase,scale",
    [([[1, 0]], [[0]], 1), ([[1]], [[1]], 1), ([[sp.I]], [[0]], 1),
     ([[1]], [[0]], 0), ([[1]], [[0]], -1), ([[1]], [[0]], sp.I)],
)
def test_invalid_form_domains(energy, phase, scale):
    with pytest.raises(ValueError):
        common_null_gram_lift(energy, phase, scale)
