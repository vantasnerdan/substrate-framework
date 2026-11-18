"""Evaluate a declared linear ionic model against the actual frozen run-end state."""
import argparse
import hashlib
import importlib.util
import json
from math import exp, pi, sqrt
from pathlib import Path
import time
import numpy as np
from scipy.constants import Boltzmann, elementary_charge, epsilon_0, physical_constants
from substrate_framework.ionic_diffusion import (darken_charge_parameters,
    diffusive_ionic_shift_eV, static_ionic_shift_eV)
from surface_io import manifest, save


def evaluate(profiles_path, transport_path):
    raw = Path(profiles_path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != 'e2e91077595de602b815ae86f45c1a675d6a553153272d0e1295b72e2a5a9943':
        raise ValueError('Not the actual frozen run-end profile source')
    source = importlib.util.spec_from_file_location('ionic_source_transport', transport_path)
    transport = importlib.util.module_from_spec(source)
    source.loader.exec_module(transport)
    record = json.loads(raw)['records'][6]
    case = record['case']
    if case['condition'] != 'PdD' or case['beam_time_s'] != 7200 or case['E_cm_keV'] != 2.5:
        raise ValueError('Wrong run-end reference state')
    temperature = case['T_K']
    isotherm = transport.Isotherm(transport.PD, temperature)
    x = np.asarray(record['state']['x'])
    z = np.asarray(record['state']['zc'])
    nmetal = np.asarray([layer[6] for layer in record['layers']])
    step = 1e-3 * np.minimum(x, isotherm.ths - x)
    slope1 = (isotherm.mu(x+step)-isotherm.mu(x-step))/(2*step)
    slope2 = (isotherm.mu(x+step/2)-isotherm.mu(x-step/2))/step
    derivative_change = float(np.max(np.abs(slope2-slope1)/np.abs(slope2)))
    if np.any(slope2 <= 0) or derivative_change > 1e-5:
        raise ValueError('Source compressibility is unstable or numerically unresolved')
    tracer = transport.PD['D0'] * exp(-transport.PD['Ea']/isotherm.kT)
    parameters = [darken_charge_parameters(float(n), float(occupancy), isotherm.ths,
                  tracer, isotherm.kT, float(slope))
                  for n, occupancy, slope in zip(nmetal, x, slope2)]
    rows = [dict(z_m=float(depth), mobile_x=float(occupancy),
                 mobile_D_density_m3=float(n*occupancy), metal_density_m3=float(n),
                 mu_slope_eV=float(slope), kappa_per_m=kappa,
                 collective_diffusivity_m2_s=collective)
            for depth, occupancy, n, slope, (kappa, collective)
            in zip(z, x, nmetal, slope2, parameters)]
    mass = physical_constants['deuteron mass'][0]
    scan = []
    for energy in (0.25, 0.5, 1.0, 2.5, 6.5):
        speed = sqrt(4 * energy * 1000 * elementary_charge / mass)
        shifts = np.asarray([diffusive_ionic_shift_eV(kappa, collective, speed)
                             for kappa, collective in parameters])
        scan.append(dict(E_cm_keV=energy, E_lab_keV=2*energy, speed_m_s=speed,
                         Ue_front_eV=float(shifts[0]),
                         Ue_max_top100nm_eV=float(shifts[z <= 100e-9].max()),
                         Ue_max_whole_foil_eV=float(shifts.max()),
                         Ue_by_cell_eV=shifts.tolist()))
    # Measured hydration content; the cited cubic reference fixes local formula-unit volume.
    oxide_temperature = 300.0
    oxide_lattice = 4.206e-10
    oxide_protons = 0.03 / oxide_lattice**3
    oxide_kappa = sqrt(elementary_charge**2 * oxide_protons /
                        (epsilon_0 * Boltzmann * oxide_temperature))
    oxide_static = static_ionic_shift_eV(oxide_kappa)
    maximum = max(row['Ue_max_whole_foil_eV'] for row in scan)
    return dict(schema='coherence-pairs-conditional-ionic-v1',
        source_profiles=manifest(profiles_path), source_transport=manifest(transport_path),
        source_driver=manifest(__file__),
        source_api=manifest(diffusive_ionic_shift_eV.__code__.co_filename),
        source_case=case, mobile_reference='actual7200s state x(z); trapped xt excluded; v not equilibrated',
        model='infinite homogeneous local linear drift-diffusion; unit ionic charge; no electronic background',
        formula='U_static*asinh(s)/s; s=v/(2*D_collective*kappa)',
        source_derivative_step_relative_change=derivative_change,
        tracer_diffusivity_m2_s=tracer, thermal_energy_eV=isotherm.kT,
        deuteron_mass_kg=mass, PdD_parameters=rows, PdD_energy_scan=scan,
        oxide=dict(material='BZY10 BaZr0.9Y0.1O3-delta', isotope='measured H; static model bound is mass-independent',
                   temperature_K=oxide_temperature, H_per_formula=0.03,
                   cited_cubic_reference_lattice_A=4.206, local_H_density_m3=oxide_protons,
                   kappa_per_m=oxide_kappa, ideal_model_static_upper_eV=oxide_static,
                   RT_mobility_status='not recovered; no finite-speed mobility or screening number invented'),
        trigger=dict(threshold_eV=100.0, PdD_max_finite_speed_model_eV=maximum,
                     oxide_ideal_static_model_upper_eV=oxide_static,
                     reached=maximum >= 100 or oxide_static >= 100,
                     further_ionic_Native_launched=False),
        literature=[dict(url='https://arxiv.org/pdf/1802.02777v2',
                        location='SectionII,Eqs1-5', use='linear drift-diffusion/Poisson equations; explicit vibration-timescale exclusion'),
                    dict(url='https://arxiv.org/pdf/1106.1924v1',
                        location='Experimental,Eq2-3', inputs='3mol%H;cubic reference a4.206A;RT QENS immobile-background assumption, not D'),
                    dict(doi='10.1140/epjb/s10051-025-00941-2',
                        use='identified inertial/phonon alternative; not parameterized or executed')],
        physical_atomic_response_qualified=False,
        applicability='Formal continuum estimate only: excludes correlations, strong fields, atomic-q and hopping/inertial response. Below-trigger is not a physical exclusion of every ionic mechanism.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profiles', required=True)
    parser.add_argument('--transport', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = time.process_time()
    result = evaluate(args.profiles, args.transport)
    result['process_cpu_s'] = time.process_time()-started
    save(args.output, result)
    print(json.dumps(dict(trigger=result['trigger'],
                         scan=[{k:v for k,v in row.items() if k != 'Ue_by_cell_eV'}
                               for row in result['PdD_energy_scan']], process_cpu_s=result['process_cpu_s'])))


if __name__ == '__main__':
    main()
