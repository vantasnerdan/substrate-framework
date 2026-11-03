"""Seed adaptive solves without evaluating a finite-domain baseline outside it.

SciPy1.10's select_initial_step probes fun(t0 + h0*direction) before bounding h0
by t_bound (installed scipy/integrate/_ivp/common.py104-121, rk.py94-100).
Use its scale-based first guess, bounded by the actual interval/max_step, through
solve_ivp's public first_step option. No Euler probe, extrapolation or tolerance
change; the existing adaptive solver still controls every accepted step.
Source: https://github.com/scipy/scipy/blob/v1.10.1/scipy/integrate/_ivp/common.py
"""
import numpy as np


def bounded_first_step(fun, interval, initial, *, rtol, atol, max_step=np.inf):
    length = abs(float(interval[1])-float(interval[0]))
    if length == 0:
        return None
    y = np.asarray(initial)
    if not y.size:
        return min(length, max_step)
    scale = atol+np.abs(y)*rtol
    d0 = np.linalg.norm(y/scale)/np.sqrt(y.size)
    d1 = np.linalg.norm(np.asarray(fun(interval[0], y))/scale)/np.sqrt(y.size)
    guess = 1e-6 if d0 < 1e-5 or d1 < 1e-5 else .01*d0/d1
    return min(length, max_step, float(guess))
