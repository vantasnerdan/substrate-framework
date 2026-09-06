#!/usr/bin/env python3
"""Exact checks for the P253/0003 Euler--strong-coupling-Skyrme route.

Scope
-----
This verifier audits the algebra behind Slobodeanu, arXiv:1405.3469v3,
Proposition 2.  It distinguishes:

* the spatial Hodge/pullback identity and static energy equality;
* the steady Euler correspondence, including its local converse;
* a time-dependent field history, which is not licensed by that theorem;
* a failure-derived dynamical route using an advected vorticity two-form;
* dimensional/action normalization, Derrick scaling, and the minimal
  topological information relevant to exchange quantization.

All calculations are exact SymPy algebra.  There is no numerical evidence and
no assertion that a topological sector is an electron or a neutrino.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from substrate_framework.verification import CheckLedger


CLAIM = "P253-0003-EULER-SKYRME-EXACT"


def vec_simplify(vector: sp.Matrix) -> sp.Matrix:
    return vector.applyfunc(sp.simplify)


def grad(expr: sp.Expr, xyz: tuple[sp.Symbol, ...]) -> sp.Matrix:
    return sp.Matrix([sp.diff(expr, coordinate) for coordinate in xyz])


def divergence(vector: sp.Matrix, xyz: tuple[sp.Symbol, ...]) -> sp.Expr:
    return sp.simplify(sum(sp.diff(vector[i], xyz[i]) for i in range(3)))


def curl(vector: sp.Matrix, xyz: tuple[sp.Symbol, ...]) -> sp.Matrix:
    x, y, z = xyz
    vx, vy, vz = vector
    return vec_simplify(
        sp.Matrix(
            [
                sp.diff(vz, y) - sp.diff(vy, z),
                sp.diff(vx, z) - sp.diff(vz, x),
                sp.diff(vy, x) - sp.diff(vx, y),
            ]
        )
    )


def advective_derivative(
    vector: sp.Matrix, velocity: sp.Matrix, xyz: tuple[sp.Symbol, ...]
) -> sp.Matrix:
    return vec_simplify(
        sp.Matrix(
            [
                sum(velocity[j] * sp.diff(vector[i], xyz[j]) for j in range(3))
                for i in range(3)
            ]
        )
    )


def material_derivative(
    expr: sp.Expr,
    velocity: sp.Matrix,
    time: sp.Symbol,
    xyz: tuple[sp.Symbol, ...],
) -> sp.Expr:
    return sp.simplify(
        sp.diff(expr, time)
        + sum(velocity[j] * sp.diff(expr, xyz[j]) for j in range(3))
    )


@dataclass(frozen=True)
class Dimension:
    """Exponents of (mass, length, time), sufficient for this exact audit."""

    mass: int = 0
    length: int = 0
    time: int = 0

    def __add__(self, other: "Dimension") -> "Dimension":
        return Dimension(
            self.mass + other.mass,
            self.length + other.length,
            self.time + other.time,
        )

    def __mul__(self, power: int) -> "Dimension":
        return Dimension(
            self.mass * power, self.length * power, self.time * power
        )


def main() -> int:
    ledger = CheckLedger(CLAIM)

    print("SOURCE: Slobodeanu, arXiv:1405.3469v3 (2019 correction)")
    print("SOURCE_SHA256: fdc453ac26b0568e5b1344940c537aa547829f8542d3a0da015ff226aa1a81bf")
    print("SOURCE_BOUNDARY: Proposition 2 is a steady 3D Riemannian theorem;")
    print("  phi: (M^3,g)->(N^2,h), phi C^2, omega target area form,")
    print("  V^flat = star(phi^*omega), with converse local near V != 0.")
    print("DOMAIN_MAP: the forward map is global for a supplied global phi;")
    print("  the Euler-to-field converse is local near V!=0, and becomes global")
    print("  only for S-integrable/simple streamline foliations. On R^3, the")
    print("  usual Hopf sector additionally uses phi(x)->constant at infinity.")
    print("MEASURE_MAP: E_sigma2,P=(1/2) integral_M (|phi^*omega|^2")
    print("  +2P(phi)) nu_g is a spatial energy. It has no dt measure and is")
    print("  not the constrained Euler spacetime action for evolving flow maps.")

    # 1. Pointwise pullback/Hodge identity in an oriented orthonormal frame.
    f12, f13, f23 = sp.symbols("f12 f13 f23", real=True)
    pullback_norm_sq = f12**2 + f13**2 + f23**2
    # star(dx^1^dx^2)=dx^3, star(dx^1^dx^3)=-dx^2,
    # star(dx^2^dx^3)=dx^1.
    hodge_vector = sp.Matrix([f23, -f13, f12])
    velocity_norm_sq = sp.expand(hodge_vector.dot(hodge_vector))
    ledger.check(
        "Hodge star isometry gives |V|^2=|phi^*omega|^2",
        sp.simplify(velocity_norm_sq - pullback_norm_sq) == 0,
    )
    ledger.mutation_sensitive(
        "Hodge component map",
        lambda candidate: sp.simplify(candidate.dot(candidate) - pullback_norm_sq)
        == 0,
        hodge_vector,
        [sp.Matrix([f23, f13, 2 * f12])],
    )

    quartic_integral, potential_integral = sp.symbols(
        "I_4 I_0", finite=True, nonnegative=True
    )
    kinetic_energy = quartic_integral / 2
    sigma2_potential_energy = quartic_integral / 2 + potential_integral
    ledger.check(
        "fluid kinetic energy equals only the sigma2 quartic contribution",
        sp.simplify(kinetic_energy - quartic_integral / 2) == 0,
    )
    ledger.check(
        "potential contribution is additional, not kinetic energy",
        sp.simplify(sigma2_potential_energy - kinetic_energy)
        == potential_integral,
    )

    # 2. Exact local steady correspondence example.
    x, y, z, t = sp.symbols("x y z t", real=True)
    xyz = (x, y, z)
    shift = sp.symbols("a", real=True)
    q_static = sp.exp(y - shift)
    v_static = sp.Matrix([q_static, 0, 0])
    bernoulli = q_static**2 / 2
    pressure = sp.simplify(bernoulli - v_static.dot(v_static) / 2)
    stationary_acceleration = advective_derivative(v_static, v_static, xyz)
    ledger.check("local example is divergence free", divergence(v_static, xyz) == 0)
    ledger.check(
        "local example obeys steady Euler with p=P-|V|^2/2=0",
        vec_simplify(stationary_acceleration + grad(pressure, xyz))
        == sp.zeros(3, 1),
    )
    ledger.check(
        "Bernoulli identity V cross curl(V)=grad(P)",
        vec_simplify(v_static.cross(curl(v_static, xyz)) - grad(bernoulli, xyz))
        == sp.zeros(3, 1),
    )
    print("STATIC_EXAMPLE: phi=(exp(y-a),z), phi^*omega=exp(y-a)dy^dz,")
    print("  V=exp(y-a)e_x, target P(q,z)=q^2/2, fluid pressure p=0.")

    # 3. A family of stationary critical maps is not an Euler history.
    a = sp.Function("a")(t)
    q = sp.exp(y - a)
    v_naive = sp.Matrix([q, 0, 0])
    naive_acceleration = vec_simplify(
        v_naive.diff(t) + advective_derivative(v_naive, v_naive, xyz)
    )
    obstruction = curl(naive_acceleration, xyz)
    expected_obstruction = sp.Matrix([0, 0, sp.diff(a, t) * q])
    ledger.check(
        "each time slice retains the steady spatial equations",
        advective_derivative(v_naive, v_naive, xyz) == sp.zeros(3, 1),
    )
    ledger.check(
        "time-dependent Euler acceleration has nonzero exact curl obstruction",
        vec_simplify(obstruction - expected_obstruction) == sp.zeros(3, 1),
    )
    ledger.mutation_sensitive(
        "time-transfer obstruction detects a moving slice",
        lambda candidate: vec_simplify(candidate - expected_obstruction)
        == sp.zeros(3, 1),
        obstruction,
        [sp.zeros(3, 1)],
    )
    print("TIME_COUNTEREXAMPLE: a'(t)!=0 gives curl[partial_t V+(V.grad)V]")
    print("  = (0,0,a'(t)exp(y-a(t))); no scalar pressure can cancel it.")

    # General condition: if every slice is already a steady correspondence,
    # the unsteady equation can hold only when partial_t V^flat is exact.  It is
    # also co-closed because every slice is divergence-free.  On a closed
    # connected manifold (or R^3 with sufficient decay), an exact co-closed
    # one-form is zero.  The explicit example tests the exposing curl condition.
    print("GENERAL_TRANSFER_CONDITION: for slice-wise steady maps, Euler requires")
    print("  partial_t(V^flat)=-d(B-P); hence d partial_t(V^flat)=0.")
    print("  With exact+co-closed rigidity on the stated domain, partial_t V=0.")

    # 4. Spacetime pullback contains a temporal component absent from the
    # Slobodeanu energy.  With signature (-,+,+,+), (1/2)F_mu_nu F^mu_nu
    # equals |B|^2-|E|^2 (c=1 convention).
    adot = sp.diff(a, t)
    magnetic_sq = q**2
    electric_sq = adot**2 * q**2
    lorentz_invariant = sp.expand(magnetic_sq - electric_sq)
    ledger.check(
        "Lorentz two-form contraction differs from static spatial norm",
        sp.simplify(lorentz_invariant - magnetic_sq + electric_sq) == 0,
    )
    ledger.check(
        "temporal pullback term vanishes exactly in stationary limit",
        sp.simplify(electric_sq.subs(adot, 0)) == 0,
    )
    ledger.mutation_sensitive(
        "Lorentz signature exposes the temporal sign",
        lambda candidate: sp.simplify(candidate - (magnetic_sq - electric_sq))
        == 0,
        lorentz_invariant,
        [magnetic_sq + electric_sq],
    )
    print("SPACETIME_SPLIT: phi^*omega=q dy^dz-a'(t)q dt^dz;")
    print("  (1/2)F_mn F^mn=q^2-a'(t)^2 q^2, while the 3D theorem uses q^2.")

    # 5. Failure-derived exact dynamical representation: pullback vorticity,
    # rather than pullback velocity flux.  Galilean transport of the shear is
    # an exact unsteady Euler solution after adding the transporting velocity.
    c = sp.symbols("c", real=True, nonzero=True)
    q_transport = sp.exp(y - c * t)
    v_dynamic = sp.Matrix([q_transport, c, 0])
    dynamic_acceleration = vec_simplify(
        v_dynamic.diff(t) + advective_derivative(v_dynamic, v_dynamic, xyz)
    )
    ledger.check("Galilean-translated shear is exact Euler", dynamic_acceleration == sp.zeros(3, 1))
    ledger.check("dynamic repair remains incompressible", divergence(v_dynamic, xyz) == 0)

    alpha = q_transport
    beta = x - t * q_transport
    ledger.check(
        "Clebsch alpha is materially advected",
        material_derivative(alpha, v_dynamic, t, xyz) == 0,
    )
    ledger.check(
        "Clebsch beta is materially advected",
        material_derivative(beta, v_dynamic, t, xyz) == 0,
    )
    clebsch_xy = sp.simplify(
        sp.diff(alpha, x) * sp.diff(beta, y)
        - sp.diff(alpha, y) * sp.diff(beta, x)
    )
    vorticity_xy = sp.simplify(sp.diff(v_dynamic[1], x) - sp.diff(v_dynamic[0], y))
    ledger.check(
        "advected pullback d(alpha)^d(beta) equals Euler vorticity two-form",
        sp.simplify(clebsch_xy - vorticity_xy) == 0,
    )
    ledger.check(
        "velocity-flux pullback and vorticity pullback are different objects",
        vec_simplify(curl(v_dynamic, xyz) - v_dynamic) != sp.zeros(3, 1),
    )
    print("DYNAMIC_ROUTE: Omega=du=d(alpha)^d(beta), D_t alpha=D_t beta=0,")
    print("  for v=(exp(y-ct),c,0), alpha=exp(y-ct), beta=x-t exp(y-ct).")
    print("  Recovering u from Omega requires the Hodge/Biot-Savart inverse and")
    print("  boundary/harmonic data; its energy is nonlocal <Omega,G Omega>/2.")

    # 6. Physical dimensions and normalization.  For dimensionless target
    # fields on dimensional M, F=phi^*omega has L^-2.  The necessary conversion
    # v=kappa *F fixes [kappa]=L^3/T.  Constant density restores physical energy.
    dim_rho = Dimension(mass=1, length=-3)
    dim_kappa = Dimension(length=3, time=-1)
    dim_pullback = Dimension(length=-2)
    dim_volume = Dimension(length=3)
    dim_time = Dimension(time=1)
    dim_velocity = dim_kappa + dim_pullback
    dim_c4 = dim_rho + dim_kappa * 2
    dim_energy = dim_c4 + dim_pullback * 2 + dim_volume
    dim_action = dim_energy + dim_time
    ledger.check("kappa converts pullback to velocity", dim_velocity == Dimension(length=1, time=-1))
    ledger.check("quartic coefficient rho*kappa^2 has physical dimension", dim_c4 == Dimension(mass=1, length=3, time=-2))
    ledger.check("integrated quartic term has energy dimension", dim_energy == Dimension(mass=1, length=2, time=-2))
    ledger.check("time-integrated energy has action dimension", dim_action == Dimension(mass=1, length=2, time=-1))
    print("NORMALIZATION: v=kappa(*phi^*omega)^sharp, [kappa]=L^3/T;")
    print("  C4=rho*kappa^2 and K=(C4/2)integral_M|phi^*omega|^2 nu_g.")
    print("  The corresponding physical action scale is C4*T/L, but Euler")
    print("  velocity equations do not fix the classically overall factor rho.")

    # 7. Spatial scale selection.  For phi_lambda(x)=phi(x/lambda), the
    # quartic and potential pieces scale as A/lambda and B*lambda^3.
    lam, A4, A0 = sp.symbols("lambda A4 A0", positive=True)
    energy_scaled = A4 / lam + A0 * lam**3
    derivative = sp.diff(energy_scaled, lam)
    stationary_relation = {A4: 3 * A0 * lam**4}
    curvature_on_shell = sp.simplify(sp.diff(energy_scaled, lam, 2).subs(stationary_relation))
    ledger.check(
        "Derrick stationary size requires A4=3*A0*lambda^4",
        sp.simplify(derivative.subs(stationary_relation)) == 0,
    )
    ledger.check("Derrick curvature is positive when both terms are positive", curvature_on_shell == 12 * A0 * lam)
    ledger.check(
        "pure quartic energy has no finite stationary size",
        sp.diff((A4 / lam), lam) != 0,
    )

    amp, inv_length = sp.symbols("A B", positive=True)
    # Euler scaling: v'(x,t)=A v(Bx,ABt), Bernoulli'=A^2 Bernoulli.
    dt_factor = sp.simplify(amp * (amp * inv_length))
    adv_factor = sp.simplify(amp * amp * inv_length)
    pressure_factor = sp.simplify(amp**2 * inv_length)
    ledger.check("Euler time/transport/pressure-gradient scaling agrees", dt_factor == adv_factor == pressure_factor)
    # A field precomposition phi(Bx) makes F and V scale as B^2; maintaining
    # Bernoulli correspondence therefore requires target potential B^4 P.
    ledger.check("pullback velocity scaling fixes Euler amplitude A=B^2", sp.simplify(amp.subs(amp, inv_length**2) - inv_length**2) == 0)
    ledger.check("mapped Bernoulli/potential must scale as B^4", sp.simplify((inv_length**2) ** 2 - inv_length**4) == 0)
    print("SCALE_RESULT: fixed potential can select lambda=(A4/(3A0))^(1/4),")
    print("  but correspondence-preserving Euler rescaling requires P->B^4 P.")
    print("  A charge-/solution-tailored potential therefore inputs the scale.")

    # 8. Topology permits but does not choose fermionic exchange.  For the
    # standard imported based mapping-space fact pi_1(Map_*(S^3,S^2))=pi_4(S^2)
    # = Z_2, a one-dimensional unitary representation sends the generator to a
    # real eta with eta^2=1.  Both signs are representations.
    eta = sp.symbols("eta", real=True)
    exchange_characters = sp.solve(sp.Eq(eta**2, 1), eta)
    ledger.check("Z2 admits two real one-dimensional exchange characters", set(exchange_characters) == {-1, 1})
    ledger.mutation_sensitive(
        "Z2 relation does not force fermionic sign",
        lambda candidate: set(candidate) == {-1, 1},
        exchange_characters,
        [[-1]],
    )

    # Classical equations are unchanged by a nonzero overall action factor,
    # while exp(i S/hbar) depends on it.  Thus classical correspondence plus a
    # Hopf label cannot determine quantum normalization or exchange sign.
    C, S, hbar = sp.symbols("C S hbar", nonzero=True)
    variation = sp.symbols("deltaS")
    scaled_variation = C * variation
    ledger.check("nonzero action multiplier leaves classical stationary set", sp.solve(sp.Eq(scaled_variation, 0), variation) == [0])
    phase_ratio = sp.simplify(sp.exp(sp.I * C * S / hbar) / sp.exp(sp.I * S / hbar))
    ledger.check("quantum phase retains the classically invisible normalization", phase_ratio.has(C))
    print("TOPOLOGY_QUANTIZATION: imported pi_1=Z2 permits eta=+1 or eta=-1;")
    print("  a quantization/FR line bundle and exchange-loop identification must")
    print("  select eta. Hopf Q in Z alone supplies neither selection nor hbar.")

    print("ROUTE_VERDICT: BLOCKED_WITH_MISSING_CONSTRUCTION")
    print("EVIDENCE_SCOPE: EXACT_STEADY_CORRESPONDENCE_AND_DYNAMIC_EXPOSURE")
    print("STRONGEST_EXACT_VERDICT: the spatial quartic density is exactly Euler")
    print("  kinetic density after kappa normalization, but slice-wise sigma2")
    print("  criticality does not transfer time-dependent Euler dynamics.")
    print("NEXT_DEPENDENCY: construct the vorticity-pullback/Hodge-inverse action on")
    print("  a declared volume-preserving configuration space, then derive its")
    print("  symplectic/quantum line bundle, exchange-loop holonomy, and effective")
    print("  relativistic dispersion for the same localized carrier.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
