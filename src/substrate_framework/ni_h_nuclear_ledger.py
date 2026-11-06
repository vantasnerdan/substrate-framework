"""Source-backed Ni/H reaction energetics, not reaction or deposited-heat rates.

AME2020 neutral-atom mass excesses are parsed as exact decimal rationals in keV.
A,Z conservation cancels the A*u term and electron *rest masses*, not electronic
binding energies. Q values therefore use neutral ground-state atomic endpoints;
bare-ion/material endpoints need their separate electronic-energy correction.
NUBASE2020 B+ denotes EC PLUS positron emission, not positron emission alone.

References: Huang et al., CPC45 030002 (2021); Wang et al., CPC45 030003
(2021); Kondev et al., CPC45 030001 (2021), especially section 2.5. The
attempt0028 source registration binds these parsers to original official bytes.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Mapping


@dataclass(frozen=True, order=True)
class Nuclide:
    a: int
    z: int

    def __post_init__(self) -> None:
        if type(self.a) is not int or type(self.z) is not int:
            raise TypeError("A and Z must be integers")
        if self.a < 1 or not 0 <= self.z <= self.a:
            raise ValueError("Require A >= 1 and 0 <= Z <= A")


H1 = Nuclide(1, 1)
NEUTRON = Nuclide(1, 0)
HE4 = Nuclide(4, 2)
NI_ISOTOPES = (58, 60, 62, 64)


def exact_decimal(field: str) -> Fraction | None:
    """AME '#' replaces the decimal point; absent quantities remain absent."""
    value = field.strip()
    if not value or '*' in value:
        return None
    return Fraction(value.replace('#', '.'))


@dataclass(frozen=True)
class MassRecord:
    nuclide: Nuclide
    element: str
    excess_kev: Fraction | None
    uncertainty_kev: Fraction | None
    excess_raw: str
    uncertainty_raw: str
    extrapolated: bool
    line_number: int
    raw_line: str


@dataclass(frozen=True)
class DecayRecord:
    nuclide: Nuclide
    half_life_raw: str
    half_life_unit: str
    half_life_uncertainty_raw: str
    spin_parity_raw: str
    branches_raw: str
    line_number: int
    raw_line: str

    @property
    def stable(self) -> bool:
        return self.half_life_raw == 'stbl'


def read_ame2020(path: str | Path) -> dict[Nuclide, MassRecord]:
    """Read the unrounded mass_1.mas20 fixed format documented in its header."""
    records = {}
    lines = Path(path).read_text(encoding='ascii').splitlines()
    if not any('f14.6,f12.6' in line for line in lines[:36]):
        raise ValueError('Not the supported AME2020 unrounded format')
    for number, line in enumerate(lines[36:], start=37):
        if not line.strip():
            continue
        n, z, a = int(line[4:9]), int(line[9:14]), int(line[14:19])
        if n + z != a:
            raise ValueError(f'Invalid AME A,N,Z at line {number}')
        key = Nuclide(a, z)
        mass, uncertainty = line[28:42], line[42:54]
        if key in records:
            raise ValueError(f'Duplicate AME isotope at line {number}')
        records[key] = MassRecord(
            key, line[20:23].strip(), exact_decimal(mass),
            exact_decimal(uncertainty), mass, uncertainty,
            '#' in mass or '#' in uncertainty, number, line,
        )
    return records


def read_nubase2020(path: str | Path) -> dict[Nuclide, DecayRecord]:
    """Ground states only; preserve uncertainties/limits and mode qualifiers."""
    records = {}
    lines = Path(path).read_text(encoding='utf-8').splitlines()
    if not lines or 'Nubase2020' not in lines[0]:
        raise ValueError('Not the supported NUBASE2020 format')
    for number, line in enumerate(lines, start=1):
        if not line or line.startswith('#'):
            continue
        if line[7] != '0':
            continue
        key = Nuclide(int(line[:3]), int(line[4:7]))
        if key in records:
            raise ValueError(f'Duplicate NUBASE ground state at line {number}')
        records[key] = DecayRecord(
            key, line[69:78].strip(), line[78:80].strip(),
            line[81:88].strip(), line[88:102].strip(), line[119:209].strip(),
            number, line,
        )
    return records


@dataclass(frozen=True)
class Reaction:
    name: str
    mechanism: str
    reactants: tuple[Nuclide, ...]
    products: tuple[Nuclide, ...]
    radiation: tuple[str, ...] = ()

    def check_conservation(self) -> None:
        for label in ('a', 'z'):
            if sum(getattr(x, label) for x in self.reactants) != sum(
                getattr(x, label) for x in self.products
            ):
                raise ValueError(f'{self.name}: {label.upper()} is not conserved')


def nickel_channels(*, include_binary_capture: bool = False) -> tuple[Reaction, ...]:
    """Kim's three proposed exits for each even Ni isotope; no paired-p state."""
    channels = []
    for a in NI_ISOTOPES:
        entrance = (Nuclide(a, 28), H1, H1)
        for name, products in (
            ('Cu+p', (Nuclide(a + 1, 29), H1)),
            ('Ni+alpha', (Nuclide(a - 2, 28), HE4)),
            ('Zn+n', (Nuclide(a + 1, 30), NEUTRON)),
        ):
            channels.append(Reaction(f'Ni{a}+2p->{name}', 'proposed Ni+2p', entrance, products))
        if include_binary_capture:
            channels.append(Reaction(
                f'Ni{a}+p->Cu{a+1}+gamma', 'separate ordinary binary radiative capture',
                (Nuclide(a, 28), H1), (Nuclide(a + 1, 29),), ('gamma/cascade',),
            ))
    return tuple(channels)


@dataclass(frozen=True)
class EnergyLedger:
    q_atomic_kev: Fraction
    coefficients: tuple[tuple[Nuclide, int], ...]
    uncertainty_l1_kev: Fraction | None
    extrapolated_inputs: tuple[Nuclide, ...]


def mass_difference(
    reactants: tuple[Nuclide, ...], products: tuple[Nuclide, ...],
    masses: Mapping[Nuclide, MassRecord],
) -> EnergyLedger:
    """Atomic mass difference; weak decays require explicit lepton bookkeeping.

    The reported L1 uncertainty is sum |coefficient|*tabulated uncertainty
    after cancellation. It upper-bounds a standard uncertainty for arbitrary
    correlations consistent with those marginal uncertainties, not a hard
    confidence interval. No covariance matrix or independence is assumed.
    """
    if sum(x.a for x in reactants) != sum(x.a for x in products):
        raise ValueError('A must balance before using mass excess differences')
    coefficients = Counter(reactants)
    coefficients.subtract(products)
    terms = tuple(sorted((key, value) for key, value in coefficients.items() if value))
    q = Fraction(0)
    uncertainty: Fraction | None = Fraction(0)
    estimated = []
    for key, coefficient in terms:
        record = masses[key]
        if record.excess_kev is None:
            raise ValueError(f'Unavailable AME mass for {key}')
        q += coefficient * record.excess_kev
        if record.uncertainty_kev is None:
            uncertainty = None
        elif uncertainty is not None:
            uncertainty += abs(coefficient) * record.uncertainty_kev
        if record.extrapolated:
            estimated.append(key)
    return EnergyLedger(q, terms, uncertainty, tuple(estimated))


def reaction_energy(reaction: Reaction, masses: Mapping[Nuclide, MassRecord]) -> EnergyLedger:
    reaction.check_conservation()
    return mass_difference(reaction.reactants, reaction.products, masses)


def _rational(value: Fraction | int | str) -> Fraction:
    if not isinstance(value, (Fraction, int, str)) or isinstance(value, bool):
        raise TypeError('Use an exact rational, integer, or decimal string, not float')
    return Fraction(value)


def two_body_kinetic_energies(
    available_kev: Fraction | int | str,
    first_rest_energy_kev: Fraction | int | str,
    second_rest_energy_kev: Fraction | int | str,
) -> tuple[Fraction, Fraction]:
    """Exact relativistic CM kinematics for declared two-body endpoints.

    Inputs are actual product rest energies and Q+entrance CM kinetic energy
    minus final excitation, with consistent electronic corrections. An AME
    atomic Q must not be combined with bare nuclear rest masses uncorrected.
    No angular distribution, reaction probability, or transport is supplied.
    A zero-mass photon is allowed, but both rest energies cannot be zero.
    """
    q, m1, m2 = map(_rational, (available_kev, first_rest_energy_kev, second_rest_energy_kev))
    if q < 0 or m1 < 0 or m2 < 0 or m1 + m2 <= 0:
        raise ValueError('Require nonnegative available energy and valid product rest energies')
    total = m1 + m2 + q
    return q * (q + 2 * m2) / (2 * total), q * (q + 2 * m1) / (2 * total)


@dataclass(frozen=True)
class DecayPath:
    mode: str
    source_branch_raw: str
    daughter: Nuclide
    nuclear_ejectiles: tuple[Nuclide, ...]
    leptons: tuple[str, ...]
    energy_convention: str
    inclusive_parent_mode: str | None = None


def decay_paths(record: DecayRecord) -> tuple[DecayPath, ...]:
    """Expand relevant NUBASE modes without inventing gamma feeding or fractions.

    B+p is a SUBSET of B+, not an independent additional branch. B+ daughter
    entries are inclusive before delayed-p emission; no network yields or
    exclusive probabilities are inferred. '?' modes remain unobserved proposals.
    """
    a, z = record.nuclide.a, record.nuclide.z
    paths = []
    for raw in record.branches_raw.split(';'):
        raw = raw.strip()
        if not raw or raw.startswith('IS='):
            continue
        mode = raw.split('=', 1)[0].split('<', 1)[0].split('>', 1)[0].split('~', 1)[0].strip().split()[0]
        if mode == 'B+':
            paths.append(DecayPath(mode, raw, Nuclide(a, z - 1), (),
                ('nu_e', 'EC consumed electron OR emitted positron; split unresolved'),
                'Q_EC=atomic parent-daughter; positron kinetic endpoint subtracts 2*m_e*c^2'))
        elif mode == 'B+p':
            paths.append(DecayPath(mode, raw, Nuclide(a - 1, z - 2), (H1,),
                ('nu_e', 'EC consumed electron OR emitted positron; split unresolved'),
                'atomic parent-daughter-H difference; positron branch subtracts 2*m_e*c^2', 'B+'))
        elif mode == 'B-':
            paths.append(DecayPath(mode, raw, Nuclide(a, z + 1), (), ('electron', 'anti_nu_e'),
                'Q_beta_minus=atomic parent-daughter'))
        elif mode == '2B+':
            paths.append(DecayPath(mode, raw, Nuclide(a, z - 2), (),
                ('2 nu_e', 'two EC/positron operations; mixture unresolved'),
                'atomic mass difference; subtract 2*m_e*c^2 for each emitted positron'))
        elif mode in ('EC', 'e+'):
            # Retained in raw data, but not extra branches on top of inclusive B+.
            if not any(x.mode == 'B+' for x in paths):
                raise ValueError(f'{mode} without preceding inclusive B+ requires explicit treatment')
        else:
            raise ValueError(f'Unsupported decay mode {raw!r} for {record.nuclide}')
    return tuple(paths)


def daughter_network(
    products: tuple[Nuclide, ...], decays: Mapping[Nuclide, DecayRecord],
) -> dict[Nuclide, tuple[DecayPath, ...]]:
    """Reachable product/decay species, including neutron decay and delayed protons.

    This is a topology, not a yield or time evolution. Stable entries may retain
    unobserved '?' modes; the source qualifier must travel with each edge.
    """
    network = {}
    pending = list(products)
    while pending:
        key = pending.pop()
        if key in network:
            continue
        paths = decay_paths(decays[key])
        network[key] = paths
        for path in paths:
            pending.append(path.daughter)
            pending.extend(path.nuclear_ejectiles)
    return network
