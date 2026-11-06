"""Read one explicit QE 7.5 printed force evaluation from authenticated bytes.

This input ledger preserves rounded decimal strings, units, atom order and
processing stages. It neither reconstructs XML forces nor defines a physical
preparation. All seven ordinary-periodic component blocks are required; extra
unassigned force rows are rejected instead of silently losing optional forces.
"""
from __future__ import annotations

from dataclasses import dataclass
from operator import index
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .qe_native_force_format import (
    BLOCK_HEADERS, PRINT_FORMAT, REQUIRED_BLOCKS, TOTAL_HEADER, UNITS,
)
from .qe_saved_native import read_bound_bytes

_DECIMAL = r"([+-]?\d+\.\d{8})"
_ROW = re.compile(r"\s*atom\s+(\d+)\s+type\s+(\d+)\s+force = (.{42})\s*")


@dataclass(frozen=True)
class PrintedForceBlock:
    name: str
    header: str
    header_line: int
    row_lines: tuple[int, ...]
    decimal_values: tuple[tuple[str, str, str], ...]
    values_ry_per_bohr: np.ndarray
    processing_stage: str
    print_format: str = PRINT_FORMAT
    units: str = UNITS

    def as_dict(self):
        return {
            "name": self.name, "header": self.header,
            "header_line": str(self.header_line),
            "row_lines": [str(line) for line in self.row_lines],
            "decimal_values": [list(row) for row in self.decimal_values],
            "values_ry_per_bohr": self.values_ry_per_bohr.tolist(),
            "processing_stage": self.processing_stage,
            "print_format": self.print_format, "units": self.units,
            "rounded_output_not_exact_zeros": True,
        }


@dataclass(frozen=True)
class NativeForceLedger:
    source: Mapping[str, str]
    atom_type_indices: tuple[int, ...]
    evaluation_index: int
    evaluation_count: int
    summary_line: int
    summary_literal: str
    blocks: Mapping[str, PrintedForceBlock]

    def as_dict(self):
        return {
            "source": dict(self.source),
            "atom_type_indices": [str(value) for value in self.atom_type_indices],
            "evaluation_index": str(self.evaluation_index),
            "evaluation_count": str(self.evaluation_count),
            "summary_line": str(self.summary_line),
            "summary_literal": self.summary_literal,
            "blocks": {name: block.as_dict() for name, block in self.blocks.items()},
            "scope": "printed component evidence; not full-precision raw/XML decomposition",
        }


def _positive_types(values):
    raw = tuple(values)
    if not raw or any(isinstance(value, (bool, np.bool_)) for value in raw):
        raise ValueError("nonempty original positive atom-type indices required")
    try:
        types = tuple(index(value) for value in raw)
    except TypeError as error:
        raise ValueError("atom-type indices must be integers") from error
    if any(value <= 0 for value in types):
        raise ValueError("atom-type indices must be positive")
    return types


def _evaluation_span(lines, selection):
    starts = [i for i, line in enumerate(lines)
              if line.strip().startswith("Forces acting on atoms")]
    if isinstance(selection, (bool, np.bool_)):
        raise ValueError("explicit nonnegative force-evaluation index required")
    try:
        selected = index(selection)
    except TypeError as error:
        raise ValueError("force-evaluation index must be an integer") from error
    if not 0 <= selected < len(starts):
        raise ValueError("force-evaluation index outside original output")
    start = starts[selected]
    if lines[start].strip() != TOTAL_HEADER:
        raise ValueError("unsupported force axes/units/header")
    stop = starts[selected + 1] if selected + 1 < len(starts) else len(lines)
    for end in range(start + 1, stop):
        if lines[end].strip().startswith("Total force ="):
            return start, end, selected, len(starts)
    raise ValueError("force evaluation lacks its terminal force summary")


def _read_block(lines, header_index, stop, atom_types):
    header = lines[header_index].strip()
    name, stage = BLOCK_HEADERS[header]
    cursor = header_index + 1
    while cursor < stop and not lines[cursor].strip():
        cursor += 1
    decimals, row_lines = [], []
    for atom, expected_type in enumerate(atom_types, start=1):
        match = _ROW.fullmatch(lines[cursor]) if cursor < stop else None
        if match is None:
            raise ValueError(f"{name}: missing or malformed atom {atom} force row")
        if int(match[1]) != atom or int(match[2]) != expected_type:
            raise ValueError(f"{name}: atom/type order differs at atom {atom}")
        fields = tuple(match[3][offset:offset + 14].strip() for offset in (0, 14, 28))
        if any(re.fullmatch(_DECIMAL, field) is None for field in fields):
            raise ValueError(f"{name}: invalid F14.8 field at atom {atom}")
        decimals.append(fields)
        row_lines.append(cursor + 1)
        cursor += 1
    values = np.array(decimals, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"{name}: nonfinite printed force")
    values.flags.writeable = False
    block = PrintedForceBlock(name, header, header_index + 1, tuple(row_lines),
                              tuple(decimals), values, stage)
    return block, cursor


def load_native_force_components(stdout_binding, *, atom_type_indices,
                                 evaluation_index) -> NativeForceLedger:
    """Consume bound stdout; atom types come from the caller's bound geometry.

    The explicit evaluation index prevents a last-iteration guess. No optional
    force block is synthesized and no printed SCF correction is discarded.
    """
    atom_types = _positive_types(atom_type_indices)
    data, source = read_bound_bytes(stdout_binding, max_bytes=64 << 20)
    lines = data.decode("utf-8", errors="strict").splitlines()
    start, stop, selected, count = _evaluation_span(lines, evaluation_index)
    blocks = {}
    cursor = start
    while cursor < stop:
        text = lines[cursor].strip()
        if text in BLOCK_HEADERS:
            block, cursor = _read_block(lines, cursor, stop, atom_types)
            if block.name in blocks:
                raise ValueError("duplicate force component: " + block.name)
            blocks[block.name] = block
        else:
            if text.startswith("atom") and "force" in text:
                raise ValueError("unassigned or extra force row at line " + str(cursor + 1))
            cursor += 1
    if blocks.keys() != REQUIRED_BLOCKS:
        raise ValueError("missing required force components: " + ",".join(sorted(REQUIRED_BLOCKS - blocks.keys())))
    return NativeForceLedger(source, atom_types, selected, count, stop + 1,
                             lines[stop].strip(), MappingProxyType(blocks))
