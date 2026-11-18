"""Expose the conditional ionic response formula against its Fourier kernel."""
import argparse
import json
from math import asinh, atan, cosh, log, pi, sqrt
from pathlib import Path
import time
import numpy as np
from scipy.integrate import quad
from scipy.constants import elementary_charge, epsilon_0
import sympy as sp
from substrate_framework.ionic_diffusion import diffusive_ionic_shift_eV, static_ionic_shift_eV


def factor_from_kernel(s):
    # Integrate q=kappa*exp(t), after the explicit angular rational integral.
    limit = max(40.0, log(max(s, 1.0))+35.0)
    if s == 0:
        function = lambda t: 1/(pi*cosh(t))
    else:
        function = lambda t: atan(s/cosh(t))/(pi*s)
    return quad(function, -limit, limit, epsabs=1e-13, epsrel=1e-11, limit=200)[0]


def verify(material):
    q, k, alpha, a, b = sp.symbols('q k alpha a b', positive=True)
    # d/dalpha atan(alpha*q/(q²+k²))/q, independently factorized at ab=k².
    rational = (q*q+a*b)/((q*q+a*a)*(q*q+b*b))
    partial = a/(a+b)/(q*q+a*a)+b/(a+b)/(q*q+b*b)
    assert sp.cancel(rational-partial) == 0
    integrals = [sp.integrate(1/(q*q+r*r), (q, 0, sp.oo)) for r in (a,b)]
    derivative = sp.simplify(a/(a+b)*integrals[0]+b/(a+b)*integrals[1])
    assert sp.simplify(derivative-sp.pi/(a+b)) == 0
    assert sp.simplify(sp.diff(sp.pi*sp.asinh(alpha/(2*k)),alpha)
                       -sp.pi/sp.sqrt(alpha*alpha+4*k*k)) == 0
    assert sp.limit(sp.pi*sp.asinh(alpha/(2*k)), alpha, 0, dir='+') == 0
    A, B, mu = sp.symbols('A B mu', positive=True)
    angular = sp.integrate(A/(A*A+B*B*mu*mu),(mu,-1,1))
    assert sp.simplify(angular-2*sp.atan(B/A)/B) == 0
    # Fourier elimination before angle integration; no fitted output enters this identity.
    w, D = sp.symbols('w D', real=True)
    kernel = 1/q**2 - (D*q*q-sp.I*w)/(q*q*(D*(q*q+k*k)-sp.I*w))
    assert sp.cancel(kernel-k*k/(q*q*(q*q+k*k-sp.I*w/D))) == 0
    samples = []
    for s in (0.0,1e-4,0.1,1.0,10.0,1e3,1e6):
        spectral = factor_from_kernel(s)
        candidate = 1.0 if s == 0 else asinh(s)/s
        relative = abs(spectral-candidate)/candidate
        assert relative < 1e-8
        samples.append(dict(speed_ratio=s, spectral=spectral, closed_form=candidate, relative_error=relative))
    assert 1e-6 < factor_from_kernel(1e6) < 1e-4
    assert abs(1.0-factor_from_kernel(1e6))/factor_from_kernel(1e6) > 1000
    # Physical units: static Coulomb constant and lab/cm speed are checked separately.
    c = elementary_charge/(4*pi*epsilon_0)
    assert abs(static_ionic_shift_eV(1e10)-c*1e10) < 1e-12
    flux_errors = []
    for point in material['PdD_parameters']:
        x = point['mobile_x']
        # Source Darken current under an electric-potential gradient gives sigma.
        sigma = (elementary_charge * point['metal_density_m3'] *
                 material['tracer_diffusivity_m2_s'] * x * (1-x) /
                 material['thermal_energy_eV'])
        relaxation = point['collective_diffusivity_m2_s'] * point['kappa_per_m']**2
        flux_errors.append(abs(relaxation-sigma/epsilon_0)/(sigma/epsilon_0))
    assert max(flux_errors) < 1e-12
    actual = []
    for row in material['PdD_energy_scan']:
        assert row['E_lab_keV'] == 2*row['E_cm_keV']
        expected_speed = sqrt(2*row['E_lab_keV']*1000*elementary_charge/material['deuteron_mass_kg'])
        assert abs(row['speed_m_s']-expected_speed)/expected_speed < 1e-12
        i = int(np.argmax(row['Ue_by_cell_eV']))
        p = material['PdD_parameters'][i]
        s = row['speed_m_s']/(2*p['collective_diffusivity_m2_s']*p['kappa_per_m'])
        independent = c*p['kappa_per_m']*factor_from_kernel(s)
        reported = row['Ue_by_cell_eV'][i]
        error = abs(independent-reported)/reported
        assert error < 1e-8
        actual.append(dict(E_cm_keV=row['E_cm_keV'], independent_max_cell_eV=independent,
                           reported_eV=reported, relative_error=error))
    assert material['source_derivative_step_relative_change'] <= 1e-5
    return dict(status='passed', symbolic=['Fourier elimination','angular integral',
                'radial partial fractions','integrated derivative and zero anchor'],
                kernel_quadrature=samples, material_crosscheck=actual,
                source_Darken_Maxwell_identity_max_relative_error=max(flux_errors),
                static_at_fast_speed_mutation_rejected=True,
                scope='conditional continuum response; no atomic applicability theorem')


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--material',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    start=time.process_time()
    result=verify(json.loads(Path(args.material).read_text()))
    result['process_cpu_s']=time.process_time()-start
    Path(args.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result))
