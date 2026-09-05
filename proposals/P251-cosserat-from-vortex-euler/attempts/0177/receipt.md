# Additive observed-response API receipt

The two canonical response functions derive the actual full-pressure
first-cell forcing/preparation and Eulerian stress/current rows. Existing
acoustic functions and their definitions are unchanged. The PDE evolution
of supplied cell functions remains the explicit caller hypothesis, not
an algebraic result returned by these functions.

First execution:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python -m pytest tests/test_euler_acoustic_response.py -q

Native first-pytest.stdout reports4 passed in2.04s, exit0. Ruff passes
on the new module/test surface. Independent Cartesian integrals expose
each pressure return and the physical/material current subtraction;
the oblique one-wave comparison checks both transverse polarizations,
including the correct microscopic pressure and full preparation.
The negative-SV limiting case is retained rather than replaced by a
positive coefficient. Invalid stationary-pressure and domain inputs reject.

SHA256 pins:

- euler_acoustic.py: 3a231339c8e3c9e50702e889438b972cc6d3f479d8f6ed4d5835130765de164c
- test_euler_acoustic_response.py: a590432d8729a38440aa58dc6a92628c198893d079aacd3b3b907d1a1cb87cec
- first-pytest.stdout: 89a1d9bbe39cb723c3552a68bf43551ff9a1f947f030a0e90a0e6ad15488d7d1

The source0170 analytic proof and independent0173 acceptance remain
unchanged. This is the reusable implementation receipt for C-CST-013,
not a parent-completion claim. Validation follows the user's blast radius;
no full scientific suite was repeated.
