"""Independent original-output decoding and semantic mutation cases for A27.

The original comparator tokenizes the source's whitespace-separated decimal
columns, independently of the production fixed-width decoder. Synthetic field
boundary and malformed-record cases are input tests, not material calculations.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile

import numpy as np

from substrate_framework.qe_native_force_ledger import load_native_force_components


HEADERS = (
    ("processed_total", "Forces acting on atoms (cartesian axes, Ry/au):"),
    ("nonlocal", "The non-local contrib.  to forces"),
    ("ion_ion", "The ionic contribution  to forces"),
    ("local", "The local contribution  to forces"),
    ("nlcc", "The core correction contribution to forces"),
    ("hubbard", "The Hubbard contrib.    to forces"),
    ("scf_correction", "The SCF correction term to forces"),
)


def independent_original_rows(lines, atom_types):
    """Decode the actual source block, without the production row/section parser."""
    result = {}
    for name, header in HEADERS:
        matches = [i for i, line in enumerate(lines) if line.strip() == header]
        if len(matches) != 1:
            raise ValueError("this original-data oracle requires one complete SCF evaluation")
        cursor = matches[0] + 1
        while not lines[cursor].strip():
            cursor += 1
        records = []
        for atom, kind in enumerate(atom_types, start=1):
            words = lines[cursor].split()
            if words[:6] != ["atom", str(atom), "type", str(kind), "force", "="] or len(words) != 9:
                raise ValueError("original independent atom/type/column interpretation failed")
            records.append((cursor + 1, tuple(words[6:])))
            cursor += 1
        result[name] = records
    return result


def compare_original(ledger, expected, require):
    for name, records in expected.items():
        block = ledger.blocks[name]
        decimals = tuple(row for _, row in records)
        values = np.array([[float(value) for value in row] for row in decimals])
        require(name + "_all_original_decimals", block.decimal_values == decimals)
        require(name + "_all_original_coordinates", np.array_equal(block.values_ry_per_bohr, values))
        require(name + "_original_row_identity", block.row_lines == tuple(line for line, _ in records))


def fixture_text(values=("-0.00000000", "0.12345678", "-1.23456789")):
    row = "     atom    1 type  1   force = " + "".join(value.rjust(14) for value in values)
    blocks = ["     " + header + "\n" + row for _, header in HEADERS]
    return "\n".join(blocks) + "\n     Total force = 0.0     Total SCF correction = 0.0\n"


def _bound_fixture(directory, name, text):
    path = Path(directory) / (name + ".stdout")
    data = text.encode("utf-8")
    path.write_bytes(data)
    return dict(path=str(path), bytes=str(len(data)), sha256=hashlib.sha256(data).hexdigest())


def exercise_boundaries(require):
    with tempfile.TemporaryDirectory(prefix="p254-native-force-input-") as directory:
        values = ("-9999.12345678", "99999.12345678", "-9999.12345678")
        wide = _bound_fixture(directory, "abutting-F14.8", fixture_text(values))
        parsed = load_native_force_components(wide, atom_type_indices=(1,), evaluation_index=0)
        require("valid_abutting_signed_F14.8", parsed.blocks["ion_ion"].decimal_values == (values,))
        zero = _bound_fixture(directory, "signed-zero", fixture_text())
        parsed = load_native_force_components(zero, atom_type_indices=(1,), evaluation_index=0)
        require("original_negative_zero_preserved", parsed.blocks["nonlocal"].decimal_values[0][0] == "-0.00000000"
                and bool(np.signbit(parsed.blocks["nonlocal"].values_ry_per_bohr[0, 0])))
        two = _bound_fixture(directory, "explicit-two-evaluations", fixture_text() + fixture_text(values))
        parsed = load_native_force_components(two, atom_type_indices=(1,), evaluation_index=1)
        require("explicit_evaluation_selects_second", parsed.evaluation_count == 2
                and parsed.blocks["local"].decimal_values == (values,))
        _semantic_mutations(directory, require)


def _semantic_mutations(directory, require):
    original = fixture_text()
    lines = original.splitlines()
    mutations = {
        "wrong_units": original.replace("Ry/au", "Ha/au", 1),
        "wrong_atom_order": original.replace("atom    1", "atom    2", 1),
        "wrong_atom_type": original.replace("type  1", "type  2", 1),
        "missing_component": "\n".join(lines[:4] + lines[6:]) + "\n",
        "missing_terminal": "\n".join(lines[:-1]) + "\n",
        "extra_force_row": "\n".join(lines[:2] + [lines[1]] + lines[2:]) + "\n",
        "duplicate_component": "\n".join(lines[:4] + lines[2:4] + lines[4:]) + "\n",
        "fortran_overflow_field": original.replace("   -0.00000000", "**************", 1),
        "wrong_decimal_precision": original.replace("0.12345678", " 0.1234567", 1),
    }
    for name, text in mutations.items():
        # Each semantic mutant has a NEW valid count/digest. Integrity refusal
        # cannot masquerade as detection of its malformed force record.
        binding = _bound_fixture(directory, name, text)
        try:
            load_native_force_components(binding, atom_type_indices=(1,), evaluation_index=0)
        except ValueError:
            refused = True
        else:
            refused = False
        require(name + "_refused", refused)
    binding = _bound_fixture(directory, "index-contract", original)
    for label, selected in (("negative", -1), ("out_of_range", 1), ("boolean", True)):
        try:
            load_native_force_components(binding, atom_type_indices=(1,), evaluation_index=selected)
        except ValueError:
            refused = True
        else:
            refused = False
        require(label + "_evaluation_refused", refused)
