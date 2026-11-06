"""Direct force origins against symbolic pairing and original full coupled probe."""
import json
import sys
from types import SimpleNamespace
import numpy as np
import sympy as sp

import saved_force_origins as force
from mixed_coordinate_response import coordinate_rhs
from run_joint_response_oracle import FiniteSystem
from substrate_framework.joint_response_observables import total_force_probe, readout
from substrate_framework.streamed_material_state import file_binding


def symbolic():
    n = sp.symbols('N', positive=True)
    c = sp.Matrix(2, 2, sp.symbols('c0:4'))
    t = sp.Matrix(2, 1, sp.symbols('t0:2'))
    u, g, v = (sp.Matrix(sp.symbols(name+'0:2')) for name in ('u', 'g', 'v'))
    a = sp.Matrix([sp.symbols('a')])
    probe_pairing = ((-c.H*v).H*(sp.sqrt(n)*u)+(-sp.sqrt(n)*t.H*v).H*a-sp.sqrt(n)*v.H*g)[0]
    physical = (sp.sqrt(n)*v.H*(-c*u-g-t*a))[0]
    assert sp.expand(probe_pairing-physical) == 0
    return 'Exact Hilbert-pairing identity for arbitrary complex C/T/coordinates/gradient; global field and sqrtN matter conversion explicit.'


def coupled():
    system = FiniteSystem()
    result = []
    for gamma in (True, False):
        r = system.response(gamma)
        stage = system.stage(.073)
        for direction in ('plus', 'minus', 'mixed'):
            zero = r.zero_tangent()
            x = r.tangent(tuple(system.random(v.shape) for v in zero.ket),
                tuple(system.random(v.shape) for v in zero.bra),
                system.random(zero.u.shape), system.random(zero.p.shape),
                system.random(zero.coordinate.shape), system.random(zero.momentum.shape))
            _, gradient, ell = r.moments(stage, x)
            physical_u, physical_p = x.u/r.sqrt_n, x.p/r.sqrt_n
            mixed = r.mixed[:, :r.nfield]
            components = force.force_components(r.curvature, mixed, physical_u,
                x.coordinate, gradient/r.sqrt_n, ncells=r.ncells)
            for axis in range(r.ndof):
                probe = total_force_probe(r, np.eye(r.ndof)[axis], identity='finite direct force origins')
                observed = readout(r, stage, x, probe)
                for name in ('electronic', 'ionic', 'field', 'total'):
                    np.testing.assert_allclose(components[name][axis], observed[name], rtol=3e-12, atol=3e-12)
            small = SimpleNamespace(u=physical_u, p=physical_p,
                coordinate=x.coordinate, momentum=x.momentum)
            rhs = coordinate_rhs(reference=SimpleNamespace(masses=r.masses,ncells=r.ncells),
                coefficients=SimpleNamespace(delta_curvature=system.dc,mixed=system.mixed,contact=system.contact),
                field=SimpleNamespace(mass=r.field_mass), charges=[system.z], polarizations=r.polarizations,
                curvature=r.curvature-system.dc, gamma=gamma, small=small,
                gradient=gradient/r.sqrt_n, field_gradient=ell/r.sqrt_n)
            np.testing.assert_allclose(components['total'], r.sqrt_n*rhs['p'], rtol=3e-12, atol=3e-12)
            wrong = force.force_components(r.curvature, mixed, physical_u,
                x.coordinate, gradient.conj()/r.sqrt_n, ncells=r.ncells)
            if np.allclose(wrong['electronic'], components['electronic'], rtol=1e-8, atol=1e-10):
                raise AssertionError('conjugated direct gradient mutation was insensitive')
            if gamma:
                wrong = force.force_components(r.curvature, r.mixed, physical_u,
                    x.coordinate/r.sqrt_n, gradient/r.sqrt_n, ncells=r.ncells)
                if np.allclose(wrong['field'], components['field'], rtol=1e-8, atol=1e-10):
                    raise AssertionError('global-field rescaling mutation was insensitive')
            result.append(dict(gamma=gamma, direction=direction, actual_coupled_origins='matched', physical_p_rhs='matched'))
    return result


if __name__ == '__main__':
    names = ('saved_force_origins', 'mixed_coordinate_response', 'run_joint_response_oracle',
        'substrate_framework.joint_bloch_response', 'substrate_framework.joint_response_observables')
    sources = [file_binding(__file__)]+[file_binding(sys.modules[name].__file__) for name in names]
    result = dict(schema='P254-saved-direct-force-origins-oracle-v1', symbolic=symbolic(), coupled=coupled())
    if sources != [file_binding(source['path']) for source in sources]:
        raise ValueError('force oracle source changed during execution')
    result.update(status='passed',source_bindings=sources,
        scope='Direct force representation/normalization/origin algebra, not material accuracy or trajectory convergence.')
    print(json.dumps(result,indent=2))
