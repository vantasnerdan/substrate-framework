"""Compatibility re-exports for the P253/0095 exact helper API.

New code imports :mod:`substrate_framework.euler_p2_principal` directly.
"""

from substrate_framework.euler_p2_principal import (  # noqa: F401
    bas_velocity_derivative,
    beltrami_hessian_symbols,
    cao_odd_border_beta,
    charged_column_slow_block,
    column_bas_matrix,
    column_metric,
    hf_detuned_return_generator,
    hf_general_first_harmonic,
    hf_growth_bracket,
    hf_hill_potential,
    hf_resonant_slow_matrix,
    maxwell_transverse_block,
    metric_bent_c1,
)
