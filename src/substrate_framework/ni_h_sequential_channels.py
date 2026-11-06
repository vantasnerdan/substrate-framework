"""Conditional ordinary Ni58 -> Cu59 -> Ni56 + alpha, never a 2p amplitude.

Two separately supplied proton encounter models compete with evaluated Cu59
decay. The finite network retains Ni59->Co59 and Ni56->Co56->Fe56. It is closed
only by an explicit caller assumption about other reactions and ground-state
preparation; this is not a complete material reaction network. Photons are not
removed from the capture vertex and Q is never automatically heat.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Mapping

from .ni_h_measured_channels import (
    ProtonEncounterModel, exact_nonnegative, rate_from_matched_flux,
)
from .ni_h_nuclear_ledger import (
    H1, HE4, Nuclide, Reaction, decay_paths, daughter_network, mass_difference,
    reaction_energy,
)

NI58, CU59, NI59, CO59, NI56, CO56, FE56 = (
    Nuclide(58, 28), Nuclide(59, 29), Nuclide(59, 28), Nuclide(59, 27),
    Nuclide(56, 28), Nuclide(56, 27), Nuclide(56, 26),
)
SPECIES = (NI58, CU59, NI59, CO59, NI56, CO56, FE56)
SPECIES_NAMES = ('Ni58', 'Cu59', 'Ni59', 'Co59', 'Ni56', 'Co56', 'Fe56')
EVENT_NAMES = ('capture', 'palpha', 'Cu59_decay', 'Ni59_decay', 'Ni56_decay', 'Co56_decay')
DECAY_EDGES = ((CU59, NI59), (NI59, CO59), (NI56, CO56), (CO56, FE56))
CAPTURE = Reaction('Ni58(p,gamma)Cu59', 'ordinary binary capture',
                   (NI58, H1), (CU59,), ('gamma cascade',))
PALPHA = Reaction('Cu59(p,alpha)Ni56', 'ordinary second binary encounter',
                  (CU59, H1), (NI56, HE4))
NET = Reaction('Ni58+2p->Ni56+alpha+gamma', 'conditional two-encounter sequence',
               (NI58, H1, H1), (NI56, HE4), ('first-capture gamma cascade',))


def sequence_energy(masses, decays) -> dict:
    """Exact evaluated neutral-atom endpoints; no photon/particle partition."""
    q1, q2, qnet = (reaction_energy(reaction, masses) for reaction in (CAPTURE, PALPHA, NET))
    if q1.q_atomic_kev + q2.q_atomic_kev != qnet.q_atomic_kev:
        raise ValueError('Intermediate mass did not cancel')
    paths = daughter_network((CU59, NI56, HE4), decays)
    decay_energies = []
    for parent, branches in sorted(paths.items()):
        for branch in branches:
            q = mass_difference((parent,), (branch.daughter,) + branch.nuclear_ejectiles, masses)
            decay_energies.append({'parent': parent, 'path': branch, 'energy': q})
    return {
        'capture': q1, 'second_encounter': q2, 'net_before_weak_decays': qnet,
        'decay_endpoints': decay_energies,
        'successful_terminal_atomic_endpoint': mass_difference((NI58, H1, H1), (FE56, HE4), masses),
        'intermediate_decay_terminal_atomic_endpoint': mass_difference((NI58, H1), (CO59,), masses),
        'endpoint_convention': ('AME neutral-atom ground-state mass differences in keV. '
            'Ionized proton/alpha/nuclear endpoints require explicit electronic binding and charge-state corrections.'),
        'step_energy_equations': [
            'E_gamma_cascade + T_Cu_recoil + E_Cu_exc = Q1_atomic + E1_CM + delta_electronic_1',
            'T_alpha + T_Ni_recoil + E_Ni_exc = Q2_atomic + E2_CM + delta_electronic_2',
            'Each equation uses its OWN CM frame. Lab totals require supplied momentum transformations and intermediate transport; E1_CM+E2_CM is not automatically one global incident energy.',
            'Weak Q_EC endpoints separately feed neutrinos, EC relaxation or positrons/annihilation, daughter recoil and cascades; subtract 2*m_e*c^2 for each beta+ kinetic endpoint.',
        ],
        'capture_photon_spectrum': None, 'exclusive_excited_state_branches': None,
        'particle_and_recoil_spectra': None, 'weak_decay_spectra': None,
        'escaping_energy': None, 'deposited_heat': None,
        'not_a_two_body_partition': 'The successful net route retains capture photons; do not partition net Q solely between alpha and Ni56.',
    }


@dataclass(frozen=True)
class SequentialPreparation:
    cascade_ground_state_justification: str
    incoherent_binary_factorization_justification: str
    competing_reactions_closure_justification: str
    decay_environment_justification: str

    def require_explicit(self) -> None:
        if not all(isinstance(value, str) and value.strip() for value in (
            self.cascade_ground_state_justification,
            self.incoherent_binary_factorization_justification,
            self.competing_reactions_closure_justification,
            self.decay_environment_justification,
        )):
            raise ValueError('Explicit preparation, factorization, closure and decay environment required')


def decay_half_lives_seconds(decays) -> dict[Nuclide, Fraction]:
    """Source central half-lives, not new fitted widths; Julian year for ky.

    NUBASE B+=100 is inclusive EC/positron topology, not a measured split.
    Applying these evaluated lifetimes to a new charge state requires caller
    justification; no fully stripped-ion EC lifetime is inferred here.
    """
    factors = {'s': Fraction(1), 'd': Fraction(86400), 'ky': Fraction(1000 * 36525 * 86400, 100)}
    result = {}
    for parent, daughter in DECAY_EDGES:
        record = decays[parent]
        branches = decay_paths(record)
        if len(branches) != 1 or branches[0].daughter != daughter or branches[0].source_branch_raw != 'B+=100':
            raise ValueError(f'Evaluated decay topology changed for {parent}')
        if record.half_life_unit not in factors:
            raise ValueError('Unsupported evaluated half-life unit')
        half = exact_nonnegative(record.half_life_raw) * factors[record.half_life_unit]
        if half == 0:
            raise ValueError('Zero half-life')
        result[parent] = half
    return result


@dataclass(frozen=True)
class ConditionalNetwork:
    capture_per_s: Fraction
    palpha_per_s: Fraction
    decay_half_lives_s: Mapping[Nuclide, Fraction]
    preparation: SequentialPreparation
    first_encounter: ProtonEncounterModel
    second_encounter: ProtonEncounterModel

    def generator(self):
        """Column generator plus integrated event-count rows; no particle heat."""
        import numpy as np
        self.preparation.require_explicit()
        hazards = [float(self.capture_per_s), float(self.palpha_per_s)]
        hazards.extend(math.log(2) / float(self.decay_half_lives_s[parent]) for parent, _ in DECAY_EDGES)
        if not all(math.isfinite(rate) and rate >= 0 for rate in hazards):
            raise ValueError('Invalid conditional hazard')
        matrix = np.zeros((len(SPECIES) + len(EVENT_NAMES),) * 2)
        edges = ((NI58, CU59), (CU59, NI56)) + DECAY_EDGES
        for event, ((parent, daughter), rate) in enumerate(zip(edges, hazards)):
            source, target = SPECIES.index(parent), SPECIES.index(daughter)
            matrix[source, source] -= rate
            matrix[target, source] += rate
            matrix[len(SPECIES) + event, source] += rate
        return matrix

    def evolve(self, initial_population: Mapping[str, str], duration_s: str) -> dict:
        """Constant external proton reservoirs; time in seconds, counts not pairs.

        Matrix exponential handles coincident hazards without dividing by their
        difference. No clipping, renormalization or mutation of negative outputs.
        The reservoir is external: protons consumed are counted, not replenished
        by an invented material supply. Unknown competing reactions remain an
        explicit closure assumption, not a physical prediction of their absence.
        """
        import numpy as np
        from scipy.linalg import expm
        if set(initial_population) != set(SPECIES_NAMES):
            raise ValueError('Supply every initial heavy-species population explicitly')
        initial = np.zeros(len(SPECIES) + len(EVENT_NAMES))
        initial[:len(SPECIES)] = [float(exact_nonnegative(initial_population[name])) for name in SPECIES_NAMES]
        time = float(exact_nonnegative(duration_s))
        if not np.isfinite(initial).all() or not math.isfinite(time):
            raise ValueError('Non-finite initial population or duration')
        final = expm(self.generator() * time) @ initial
        counts = dict(zip(EVENT_NAMES, final[len(SPECIES):].tolist()))
        return {
            'conditional': True, 'duration_s': duration_s,
            'populations': dict(zip(SPECIES_NAMES, final[:len(SPECIES)].tolist())),
            'events': counts, 'alpha_particles_produced': counts['palpha'],
            'protons_consumed': counts['capture'] + counts['palpha'],
            'completed_sequences_from_new_capture': counts['palpha'] if exact_nonnegative(initial_population['Cu59']) == 0 else None,
            'multiplicity_note': 'palpha events can also use initially supplied Cu59; a capture followed by decay consumes only its first proton.',
            'photon_cascade_events': counts['capture'],
            'photon_count': None, 'emitted_energy': None, 'escaping_energy': None, 'deposited_heat': None,
            'scope': 'Finite conditional seven-heavy-species network, not complete N-node material evolution',
        }


def build_conditional_network(tables, decays, first: ProtonEncounterModel,
                              second: ProtonEncounterModel,
                              preparation: SequentialPreparation) -> ConditionalNetwork:
    preparation.require_explicit()
    first_rows = [row for row in tables['capture'] if row['target_a'] == 58]
    first_rate = rate_from_matched_flux(first_rows, first)
    second_rate = rate_from_matched_flux(tables['cu59_palpha'], second)
    return ConditionalNetwork(first_rate.per_target_per_second, second_rate.per_target_per_second,
                              decay_half_lives_seconds(decays), preparation, first, second)


ENERGY_CARRIERS = ('capture_photons', 'alpha_kinetic', 'nuclear_recoil',
                   'decay_photons_and_electrons_positrons', 'neutrinos', 'remaining_excitation')

def account_energy(available_kev: str, produced_kev: Mapping[str, str], *,
                   transport: Mapping[str, Mapping[str, str]] | None = None,
                   transport_provenance: str | None = None) -> dict:
    """Check a SUPPLIED common-frame energy budget; do not infer one from Q.

    Each carrier must remain separate. Transport entries specify deposited,
    escaped and retained energy in keV; retained energy is not deposited heat.
    This validates accounting, not the physical truth of a user's transport model.
    """
    if set(produced_kev) != set(ENERGY_CARRIERS):
        raise ValueError('Supply every radiation/particle/recoil/excitation energy separately')
    produced = {key: exact_nonnegative(value) for key, value in produced_kev.items()}
    available = exact_nonnegative(available_kev)
    if sum(produced.values(), Fraction(0)) != available:
        raise ValueError('Produced energy does not match supplied common-frame budget')
    result = {'produced_kev': produced, 'deposited_heat_kev': None,
              'escaped_kev': None, 'retained_kev': None}
    if transport is None:
        if transport_provenance is not None:
            raise ValueError('Provenance alone is not a transport model')
        return result
    if not transport_provenance or not transport_provenance.strip():
        raise ValueError('Transport/deposition provenance required; Q is not heat')
    if set(transport) != set(ENERGY_CARRIERS):
        raise ValueError('Transport must cover every carrier')
    totals = {key: Fraction(0) for key in ('deposited_kev', 'escaped_kev', 'retained_kev')}
    for carrier, energy in produced.items():
        if set(transport[carrier]) != set(totals):
            raise ValueError('Require deposited, escaped and retained energy separately')
        partition = {key: exact_nonnegative(value) for key, value in transport[carrier].items()}
        if sum(partition.values(), Fraction(0)) != energy:
            raise ValueError('Carrier transport energy does not close')
        for key, value in partition.items():
            totals[key] += value
    result.update(deposited_heat_kev=totals['deposited_kev'], escaped_kev=totals['escaped_kev'],
                  retained_kev=totals['retained_kev'], transport_provenance=transport_provenance)
    return result
