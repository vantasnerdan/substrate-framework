"""Independent finite-model scalar, exact differentiation and Fourier oracles.

No propagation or material model is represented by these finite fixtures.
The scalar below is written literally, not in terms of the implementation's
mixed-product helper. Finite differences use physical η=λ* perturbations.
"""
from types import SimpleNamespace as NS

import numpy as np
import sympy as sy


FIELDS = ("density", "potential", "u", "p", "shift", "momentum", "gradient",
          "field_gradient_electronic")
KEYS = ("mean", "plus", "minus", "mixed")


def literal_scalar(r, c, f, z, s, electronic, dv):
    """Literal physical scalar, separately encoded for independent differentiation."""
    charge_shift = np.repeat(z, 3)*np.tile(s.shift, len(z))
    mechanical = s.p-charge_shift
    return dict(
        bare_electronic=r.ncells*electronic,
        hartree=r.ncells*dv*sum((s.density-r.nbar)*s.potential)/2,
        electron_ionic=r.ncells*(s.u @ s.gradient),
        ionic_reference_linear=-r.ncells*(s.u @ (r.ionic_force+r.support_force)),
        ionic_C=r.ncells*(s.u @ r.curvature @ s.u)/2,
        ionic_kinetic_minimal=r.ncells*sum(mechanical*mechanical/r.masses)/2,
        delta_C=r.ncells*(s.u @ c.delta_curvature @ s.u)/2,
        electron_field=r.ncells*(s.shift @ s.field_gradient_electronic),
        field_contact=r.ncells*(s.shift @ c.contact @ s.shift)/2,
        ion_field_mixed=r.ncells*(s.u @ c.mixed @ s.shift),
        maxwell=(s.momentum @ s.momentum)/(2*f.mass))


def fixture(*, exact=False, gamma=True):
    """Deterministic all-nonzero, nonstationary inputs; no fitted energy residual."""
    if exact:
        count = iter(range(1, 1000))
        def array(shape):
            return np.array([sy.Rational(next(count) % 13+1, 17)
                             + sy.I*sy.Rational(next(count) % 7+1, 19)
                             for _ in range(int(np.prod(shape)))], dtype=object).reshape(shape)
        scalar = lambda x: sy.Rational(str(x))
    else:
        rng = np.random.default_rng(74121)
        array = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
        scalar = float
    symmetric = lambda a: (a+a.T)/2
    r = NS(ncells=5, masses=np.array([2, 3, 5], dtype=object if exact else float),
           nbar=array((2,)), ionic_force=array((3,)), support_force=array((3,)),
           curvature=symmetric(array((3, 3))))
    c = NS(delta_curvature=symmetric(array((3, 3))),
           contact=symmetric(array((3, 3))), mixed=array((3, 3)))
    f = NS(mass=scalar(7))
    stages = {key: NS(**{name: array((2,) if name in ("density", "potential") else (3,))
                         for name in FIELDS}) for key in KEYS}
    one = {key: array((1,))[0] for key in KEYS}
    if not exact:
        for obj in (r, c):
            for name, value in vars(obj).items():
                if isinstance(value, np.ndarray):
                    setattr(obj, name, value.real)
        for key in ("mean", "mixed"):
            for name in FIELDS:
                setattr(stages[key], name, getattr(stages[key], name).real)
            one[key] = one[key].real
        for name in FIELDS:
            setattr(stages["minus"], name, getattr(stages["plus"], name).conj())
        one["minus"] = one["plus"].conjugate()
    if not gamma:
        for key in ("plus", "minus"):
            stages[key].shift = np.zeros(3)
            stages[key].momentum = np.zeros(3)
    return dict(reference=r, coefficients=c, field=f, charges=np.array([scalar(2)]),
                stages=stages, one_body=one, grid_measure=scalar(0.3),
                curvature_plus=r.curvature.copy(), curvature_minus=r.curvature.copy(),
                gamma=gamma)


def exact_derivatives(args, authoritative):
    """Differentiate each literal term exactly; also bind to original scalar."""
    lam, eta = sy.symbols("lambda eta")
    weights = (1, lam, eta, lam*eta)
    ss = args["stages"]
    s = NS(**{name: sum((w*getattr(ss[key], name) for w, key in zip(weights, KEYS)),
                       np.zeros_like(getattr(ss["mean"], name))) for name in FIELDS})
    electronic = sum(w*args["one_body"][key] for w, key in zip(weights, KEYS))
    literal = literal_scalar(args["reference"], args["coefficients"], args["field"],
                             args["charges"], s, electronic, args["grid_measure"])
    source = authoritative(reference=args["reference"], coefficients=args["coefficients"],
                           field=args["field"], charges=args["charges"], stage=s,
                           one_body=electronic, grid_measure=args["grid_measure"])
    source_defects = {key: sy.simplify(literal[key]-source[key]) for key in literal}
    mixed = {key: sy.diff(value, lam, eta).subs({lam: 0, eta: 0})
             for key, value in literal.items()}
    return mixed, source_defects


def periodic_model(args, n=5, mode=1):
    """C(R-S) with real reciprocal offsite blocks, not transpose-blind data."""
    rng = np.random.default_rng(132)
    kernel = np.zeros((n, 3, 3))
    kernel[1] = rng.normal(size=(3, 3))
    kernel[-1] = kernel[1].T
    kernel[0] = args["reference"].curvature-kernel[1]-kernel[-1]
    angle = 2*np.pi*mode*np.arange(n)/n
    phase = np.exp(1j*angle)
    cp = sum(kernel[d]*np.exp(-2j*np.pi*mode*d/n) for d in range(n))
    cm = sum(kernel[d]*np.exp(+2j*np.pi*mode*d/n) for d in range(n))
    dense = np.block([[kernel[(row-col) % n] for col in range(n)] for row in range(n)])
    return phase, cp, cm, dense


def periodic_scalar(args, phase, dense, lam, eta):
    """Dense real-space energy of a smooth physical finite perturbation family.

    Cubic λ²η/λη² stage terms have zero mixed derivative but give the physical
    central Laplacian a genuine O(h²) remainder. No production formula is used.
    """
    r, c, f = args["reference"], args["coefficients"], args["field"]
    ss = args["stages"]
    n = len(phase)
    total = None
    displacements = []
    for cell in range(n):
        w = (1, lam*phase[cell], eta*phase[cell].conjugate(), lam*eta)
        values = {}
        for name in FIELDS:
            if name in ("shift", "momentum"):
                ws = (1, lam, eta, lam*eta)
                cubic = (lam*lam*eta+lam*eta*eta)*0.13*getattr(ss["mixed"], name)
            else:
                ws = w
                cubic = 0.13*(lam*lam*eta*phase[cell]*getattr(ss["plus"], name)
                              +lam*eta*eta*phase[cell].conjugate()*getattr(ss["minus"], name))
            values[name] = sum(weight*getattr(ss[key], name) for weight, key in zip(ws, KEYS))+cubic
        stage = NS(**values)
        electronic = sum(weight*args["one_body"][key] for weight, key in zip(w, KEYS))
        local = literal_scalar(NS(**{**vars(r), "ncells": 1}), c, f, args["charges"],
                               stage, electronic, args["grid_measure"])
        local["maxwell"] /= n
        displacements.append(stage.u)
        if total is None:
            total = {key: 0j for key in local}
        for key, value in local.items():
            total[key] += value
    u = np.concatenate(displacements)
    total["ionic_C"] = u @ dense @ u/2
    return total


def physical_laplacian(energy, h):
    """∂λ∂η = (∂x²+∂y²)/4 for λ=x+iy, η=x-iy."""
    center = energy(0, 0)
    points = [energy(h, h), energy(-h, -h), energy(1j*h, -1j*h), energy(-1j*h, 1j*h)]
    return {key: (sum(point[key] for point in points)-4*center[key])/(4*h*h)
            for key in center}


def local_impulse_normalization(production):
    """Derive physical amplitudes from the explicit orthonormal DFT matrix.

    xhat=F x; reconstruction x=F† xhat. Per-cell Fourier amplitudes are the
    coefficients multiplying exp(iqR), hence xhat/sqrt(N), not xhat itself.
    Complete-q mixed energy sums to twice the quadratic physical dose because
    a mixed DERIVATIVE has two ordered legs. The energy coefficient is half
    that derivative only when reconstructing a real one-parameter dose.
    """
    n, cell, kick, mass = 7, 2, 0.8, 3.0
    sites = np.arange(n)
    f = np.exp(-2j*np.pi*np.outer(sites, sites)/n)/np.sqrt(n)
    impulse = np.zeros(n)
    impulse[cell] = kick
    collective = f @ impulse
    reconstruction_basis = np.exp(2j*np.pi*np.outer(sites, sites)/n)
    # Solve the explicit reconstruction problem, rather than assume a factor.
    physical = np.linalg.solve(reconstruction_basis, impulse)
    dose = impulse @ impulse/(2*mass)
    mixed_q = n*physical*physical.conj()/mass
    args = fixture(gamma=False)
    args["reference"].ncells = n
    args["reference"].masses[0] = mass
    for stage in args["stages"].values():
        for name in FIELDS:
            setattr(stage, name, np.zeros(getattr(stage, name).shape, dtype=complex))
    args["one_body"] = dict.fromkeys(KEYS, 0.0)
    implementation = []
    for amplitude in physical:
        args["stages"]["plus"].p[0] = amplitude
        args["stages"]["minus"].p[0] = amplitude.conjugate()
        implementation.append(production(**args)["ionic_kinetic_minimal"])
    return dict(ncells=n, cell=cell, kick=kick, mass=mass,
                unitarity_defect=float(np.max(np.abs(f.conj().T@f-np.eye(n)))),
                reconstruction_defect=float(np.max(np.abs(reconstruction_basis@physical-impulse))),
                derived_collective_to_physical_defect=float(np.max(np.abs(physical-collective/np.sqrt(n)))),
                derived_local_amplitude_defect=float(np.max(np.abs(physical-kick*np.exp(-2j*np.pi*sites*cell/n)/n))),
                parseval_defect=float(abs(np.vdot(collective, collective)-np.vdot(impulse, impulse))),
                direct_initial_dose=float(dose), complete_q_half_mixed_sum=float(np.sum(mixed_q).real/2),
                complete_q_dose_defect=float(abs(np.sum(mixed_q)/2-dose)),
                implementation_complete_q_dose_defect=float(abs(sum(implementation)/2-dose)),
                implementation_mode_derivative_defect=float(np.max(np.abs(np.array(implementation)-mixed_q))),
                wrong_collective_units_dose_gap=float(abs(np.sum(n*collective*collective.conj()/mass)/2-dose)),
                missing_energy_half_gap=float(abs(np.sum(mixed_q)-dose)))
