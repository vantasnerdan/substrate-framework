"""Dirac hopping -> plaquette flux -> selected 25-state colour-spin cell.

Finite mechanism calculation, not QCD or a fitted spectrum. The canonical
1+3+bar3 SU3 link truncation is used for actual vacuum loop matrix elements;
classical unitary links are a separate Dirac/commutator coefficient consumer.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import expm

from rich_composite_local_flux_binding_20261017_model import link_operators
from rich_nucleon_hierarchy_20261019 import COLOUR, SIGMA, matter_sector

HERE = Path(__file__).resolve().parent
Wave = dict[tuple[int, ...], complex]


def number(value: complex) -> dict[str, float]:
    return {'real': float(np.real(value)), 'imag': float(np.imag(value))}


def inner(left: Wave, right: Wave) -> complex:
    return sum(value.conjugate() * right.get(key, 0j) for key, value in left.items())


def dirac_map() -> dict[str, Any]:
    """Exact finite two-direction Pauli identity and oriented-loop comparison."""
    a, q, mass, theta, phi = 1., .7, 3., .31, .43
    ux, uy = expm(1j*theta*COLOUR[0]), expm(1j*phi*COLOUR[1])
    # Covariant shifts on a 2x2 periodic lattice; tensor order site, colour.
    tx = np.kron(SIGMA[0], np.eye(2))
    ty = np.kron(np.eye(2), SIGMA[0])
    sx, sy = np.kron(tx, ux), np.kron(ty, uy)
    px, py = (sx-sx.conj().T)/(2j*a), (sy-sy.conj().T)/(2j*a)
    comm = px@py-py@px
    b = comm/(1j*q)
    pauli = np.kron(SIGMA[0], px)+np.kron(SIGMA[1], py)
    spinless = np.kron(np.eye(2), px@px+py@py)
    magnetic = -q*np.kron(SIGMA[2], b)
    # One local closed loop transforms by conjugation at its base vertex.
    def loop_b(x: Any, y: Any) -> Any:
        w = x@y@x.conj().T@y.conj().T
        raw = 1j*(w-w.conj().T)/(2*q*a*a)
        return raw-np.trace(raw)*np.eye(3)/3
    local_b = loop_b(ux, uy)
    coefficients = [float((2*np.trace(t@local_b)).real) for t in COLOUR]
    v = expm(1j*(.27*COLOUR[2]-.18*COLOUR[7]))
    covariance = np.linalg.norm(loop_b(v@ux@v.conj().T, v@uy@v.conj().T)
                                -v@local_b@v.conj().T)
    # Low-momentum uniform spatial mode: shifts have eigenvalue +1.
    low_px, low_py = (ux-ux.conj().T)/(2j*a), (uy-uy.conj().T)/(2j*a)
    low_b = (low_px@low_py-low_py@low_px)/(1j*q)
    scale_rows = []
    for scale in (1., .5, .25):
        x, y = expm(1j*scale*theta*COLOUR[0]), expm(1j*scale*phi*COLOUR[1])
        p1, p2 = (x-x.conj().T)/(2j*a), (y-y.conj().T)/(2j*a)
        discrete_b = (p1@p2-p2@p1)/(1j*q)
        scale_rows.append({'angle_scale': scale, 'absolute_loop_commutator_difference':
                           float(np.linalg.norm(loop_b(x,y)-discrete_b)),
                           'relative_difference': float(np.linalg.norm(loop_b(x,y)-discrete_b)
                                                        /np.linalg.norm(discrete_b))})
    # A fixed fundamental colour with off-diagonal physical spin matrix element:
    # choose x-polarized states |+x> -> |-x> for the z-directed field.
    spin_plus, spin_minus = np.array([1,1])/np.sqrt(2), np.array([1,-1])/np.sqrt(2)
    initial, final = np.kron(spin_plus, [1,0,0]), np.kron(spin_minus, [1,0,0])
    pauli_entry = np.vdot(final, -q*np.kron(SIGMA[2],low_b)@initial)/(2*mass)
    loop_entry = np.vdot(final, -q*np.kron(SIGMA[2],local_b)@initial)/(2*mass)
    commuting = loop_b(ux, expm(1j*phi*COLOUR[0]))
    zero = loop_b(np.eye(3),np.eye(3))
    return {'inputs': {'spacing_a':a, 'gauge_charge_q':q, 'Dirac_mass_m':mass,
                       'theta_T0':theta, 'phi_T1':phi},
            'finite_site_count':4, 'positive_Dirac_block_dimension':24,
            'Pauli_square_identity_residual':float(np.linalg.norm(pauli@pauli-spinless-magnetic)),
            'finite_covariant_commutator_norm':float(np.linalg.norm(comm)),
            'closed_loop_B_adjoint_coefficients':coefficients,
            'loop_B_norm':float(np.linalg.norm(local_b)),
            'loop_base_covariance_residual':float(covariance),
            'commutator_Pauli_transition':number(pauli_entry),
            'local_loop_Pauli_transition':number(loop_entry),
            'spinless_transition':number(np.vdot(final,np.kron(np.eye(2),low_px@low_px+low_py@low_py)@initial)),
            'flat_loop_B_norm':float(np.linalg.norm(zero)),
            'commuting_loop_B_norm':float(np.linalg.norm(commuting)),
            'loop_vs_central_difference':scale_rows,
            'meaning':'Finite central-difference B is two-hop and not exactly a local Wilson B at finite angles. '
                      'Their common smooth-link leading coefficient is -q/(2m). The loop-local Pauli action '
                      'used below is the corresponding Wilson/clover effective interaction, not an exact '
                      'identity for the finite central-difference hopping Hamiltonian.'}


def loop_terms(generator: Any, edges: tuple[int, int, int, int], a: float, q: float
               ) -> list[tuple[complex, tuple[tuple[int, int, int, bool], ...]]]:
    """B^c=i/(q a^2) tr[T^c(W-Wdag)]; W=U1 U2 U3dag U4dag."""
    terms = []
    for row, col in zip(*np.nonzero(generator)):
        for i,j,k in itertools.product(range(3), repeat=3):
            coefficient = complex(generator[row,col])
            indices = [(col,i,False), (i,j,False), (k,j,True), (row,k,True)]
            ops = tuple((edge, int(x), int(y), adj) for edge,(x,y,adj) in zip(edges,indices))
            terms.append((1j*coefficient/(q*a*a),ops))
            terms.append((-1j*coefficient.conjugate()/(q*a*a),
                          tuple((edge,x,y,not adj) for edge,x,y,adj in ops)))
    return terms


def act(terms: Any, wave: Wave, transitions: Any) -> Wave:
    result: dict[tuple[int, ...],complex] = defaultdict(complex)
    for coefficient, ops in terms:
        for state, amplitude in wave.items():
            changed, weight = list(state), coefficient*amplitude
            for edge,x,y,adj in ops:
                hop = transitions[x,y,adj].get(changed[edge])
                if hop is None:
                    break
                changed[edge], value = hop
                weight *= value
            else:
                result[tuple(changed)] += weight
    return {key:value for key,value in result.items() if abs(value)>1e-15}


def quantum_flux() -> dict[str, Any]:
    """Actual shared nine-link, three-plaquette source Gram and electric gap."""
    a, q, mass, electric = 1., .7, 3., .27
    transport,left,right = link_operators()
    transitions = {}
    for key, matrix in transport.items():
        entries = {}
        for col in range(19):
            rows = np.flatnonzero(matrix[:,col])
            if len(rows)>1:
                raise ValueError('Expected at most one canonical transition per link matrix entry')
            if len(rows):
                entries[col] = (int(rows[0]),complex(matrix[rows[0],col]))
        transitions[key] = entries
    # Root-positive x,y,z links 0,1,2 are shared. Six outer edges are distinct.
    # Planes yz,zx,xy correspond to physical axial-vector components x,y,z.
    planes = [(1,5,6,2),(2,7,8,0),(0,3,4,1)]
    vacuum = {(0,)*9:1.+0j}
    terms = [[loop_terms(t,plane,a,q) for t in COLOUR] for plane in planes]
    sources = [act(terms[k][c],vacuum,transitions) for c in range(8) for k in range(3)]
    gram = np.array([[inner(x,y) for y in sources] for x in sources])
    beta2 = 1/(3*q*q*a**4)
    gap = 16*electric/3
    # All nonvacuum link blocks are fundamental or antifundamental: C_F=4/3.
    casimir = sum((x@x for x in left),np.zeros((19,19),complex))
    casimir_expected = np.diag([0.]+[4/3]*18)
    energy_residual = max(np.sqrt(sum(abs(value)**2*(electric*sum(4/3 for n in state if n)-gap)**2
                                     for state,value in wave.items())) for wave in sources)
    # Compute all same-plane excited-to-excited B matrix elements. Cross-plane
    # triple products vanish because each plane has two unique outer links;
    # a single U/Udag on such an edge cannot return that link to its vacuum.
    excited_max = 0.
    for c in range(8):
        for b in range(8):
            image = act(terms[2][c],sources[3*b+2],transitions)
            excited_max = max(excited_max, max(abs(inner(sources[3*d+2],image))/beta2 for d in range(8)))
    coupling = -q*np.sqrt(beta2)/(2*mass)
    kappa = coupling**2/gap
    data, levels = matter_sector(3)
    n = next(level for level in levels if level['S']==.5)
    delta = next(level for level in levels if level['S']==1.5)
    # Pauli -> loop source -> electric inverse -> return, on the actual parent G.
    def shift(level: Any) -> float:
        return -q*q/(4*mass*mass)*beta2/gap*sum(
            float(np.vdot(op@level['psi'],op@level['psi']).real) for op in data['mixed'])
    en, ed = shift(n), shift(delta)
    # All three root edges have coordinates directed outward. In the canonical
    # convention their root charges are -L, including a link used as Udag
    # on the returning side of an oriented loop.
    root_charge = []
    for generator in range(8):
        components = {}
        for edge, matrix in [(0,-left[generator]),(1,-left[generator]),(2,-left[generator])]:
            components[edge] = matrix
        root_charge.append(components)
    def root_act(wave: Wave, components: Any) -> Wave:
        out: dict[tuple[int,...],complex] = defaultdict(complex)
        for state,value in wave.items():
            for edge,matrix in components.items():
                for row in np.flatnonzero(matrix[:,state[edge]]):
                    changed=list(state); changed[edge]=int(row)
                    out[tuple(changed)] += matrix[row,state[edge]]*value
        return dict(out)
    # At root each source includes outgoing links + inward dagger: both use -L
    # because their coordinate edges are directed outward from the root.
    covariance_residual = 0.
    for c,t in enumerate(COLOUR):
        for b,u in enumerate(COLOUR):
            image=root_act(sources[3*b+2],root_charge[c])
            # [t,u]=i f_cbd T_d and root action on B^b|0> is -i f_cbd B^d|0>.
            expected: dict[tuple[int,...],complex] = defaultdict(complex)
            for d,v in enumerate(COLOUR):
                coefficient=-2*np.trace((t@u-u@t)@v)
                for state,value in sources[3*d+2].items():
                    expected[state] += coefficient*value
            keys=set(image)|set(expected)
            covariance_residual=max(covariance_residual,np.sqrt(sum(abs(image.get(s,0)-expected.get(s,0))**2 for s in keys)))
    # Full locality consumer: six nonroot vertices must remain Gauss singlets.
    edge_vertices = [(0,1),(0,2),(0,3),(1,4),(2,4),(2,5),(3,5),(3,6),(1,6)]
    nonroot_gauss = 0.
    for vertex in range(1,7):
        for c in range(8):
            components = {edge: (-left[c] if start==vertex else -right[c])
                          for edge,(start,end) in enumerate(edge_vertices)
                          if start==vertex or end==vertex}
            for wave in sources:
                image = root_act(wave,components)
                nonroot_gauss = max(nonroot_gauss,float(np.sqrt(inner(image,image).real)))
    return {'inputs':{'spacing_a':a,'gauge_charge_q':q,'Dirac_mass_m':mass,'electric_Casimir_coefficient':electric},
            'link_representation':'canonical 1+3+bar3 Peter-Weyl, dimension19 per edge',
            'shared_link_count':9,'oriented_plaquettes_xyz':planes,
            'source_order':'colour a outer, axial component k inner; same as parent mixed G',
            'source_support_sizes':[len(x) for x in sources],
            'source_Gram_diagonal_predicted':beta2,
            'source_Gram_eigenvalue_min':float(np.linalg.eigvalsh(gram)[0]),
            'source_Gram_eigenvalue_max':float(np.linalg.eigvalsh(gram)[-1]),
            'source_Gram_prediction_residual':float(np.linalg.norm(gram-beta2*np.eye(24))),
            'root_adjoint_source_residual':float(covariance_residual),
            'nonroot_six_vertex_Gauss_max_norm':nonroot_gauss,
            'link_electric_Casimir_residual':float(np.linalg.norm(casimir-casimir_expected)),
            'loop_source_electric_gap':gap,'electric_eigen_residual':float(energy_residual),
            'same_plane_projected_B_excited_block_max':float(excited_max),
            'cross_plane_excited_block':'zero by unique outer-edge vacuum selection, including shared root edges',
            'selected_loop_space_dimension':25,
            'parent_coupling_g_signed':coupling,'parent_gap_omega':gap,
            'elimination_g2_over_omega':kappa,
            'nucleon_G2':n['lambda_G2'],'Delta_G2':delta['lambda_G2'],
            'nucleon_second_order_shift':en,'Delta_second_order_shift':ed,
            'Delta_minus_nucleon_second_order':ed-en,
            'matching':'g=-q beta/(2m)=-1/(2 sqrt(3) m a^2), omega=16 kE/3; '
                       'g^2/omega=1/(64 m^2 a^4 kE). Electric-positive loop states supply '
                       'the parent vacuum+24 source multiplet, not an added oscillator.',
            'scope':'The electric Hamiltonian is an assumed positive canonical gauge-link dynamics. '
                    'The 25-state subspace is a selected source projection; B also accesses higher '
                    'link sectors outside it. Exact source Gram/gap and order-Pauli^2 elimination '
                    'do not assert a nonperturbative invariant 25-state truncation.'}


def main() -> None:
    result = {'status':'executed finite Dirac and canonical quantum-loop supplier',
              'runtime':{'python':platform.python_version(),'numpy':np.__version__},
              'dirac_hopping_map':dirac_map(),'quantum_link_flux':quantum_flux(),
              'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'dependency_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                  for name in ['rich_nucleon_hierarchy_20261019.py',
                               'rich_composite_local_flux_binding_20261017_model.py']},
              'remaining_inputs':['canonical SU3 link representation and positive electric law',
                                  'Dirac kinetic normalization 1/(2a) and mass m',
                                  'Wilson/clover local Pauli effective action or smooth-link expansion',
                                  'spatial transport of composite loop states and Gauss-compatible multi-cell projection'],
              'carrier_connection':'Replace each independent local adjoint-spin cell by normalized '
                  'B_x,ak[U]|electric vacuum>, retain the existing canonical U quark hopping and '
                  'Peierls factor. At order Pauli^2 the same CAR G gives the computed nucleon/Delta '
                  'shifts. Existing field propagation h must instead be computed from shared-loop '
                  'electric/Wilson matrix elements; no value of h is supplied here.'}
    target=HERE/'rich_colour_spin_root_20261019_result.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
