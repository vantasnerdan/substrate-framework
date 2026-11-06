"""QE 7.5 forces.f90 print vocabulary and stage metadata, not a force model.

Source: PW/src/forces.f90:264-399,532; force_us.f90:190-197.
The total is printed before if_pos masking; forcescc is printed after masking.
No component printed as zero is asserted to vanish at full precision.
"""
from types import MappingProxyType

TOTAL_HEADER = "Forces acting on atoms (cartesian axes, Ry/au):"
PRINT_FORMAT = "3F14.8"
UNITS = "Ry/bohr"

# name, processing stage; do not infer these stages from the numerical values.
BLOCK_HEADERS = MappingProxyType({
    TOTAL_HEADER: (
        "processed_total",
        "assembled and translation/symmetry/external processed; before if_pos masking",
    ),
    "The non-local contrib.  to forces": (
        "nonlocal",
        "nonlocal component already internally symmetrized by force_us",
    ),
    "The ionic contribution  to forces": (
        "ion_ion",
        "periodic ionic component before total-force translation/symmetry processing",
    ),
    "The local contribution  to forces": (
        "local",
        "local electronic component before total-force translation/symmetry processing",
    ),
    "The core correction contribution to forces": (
        "nlcc",
        "core-correction component, not an inferred exact zero",
    ),
    "The Hubbard contrib.    to forces": (
        "hubbard",
        "Hubbard component, not an inferred exact zero",
    ),
    "The SCF correction term to forces": (
        "scf_correction",
        "incomplete-SCF correction after if_pos masking; not a bare preparation force",
    ),
})
REQUIRED_BLOCKS = frozenset(value[0] for value in BLOCK_HEADERS.values())
