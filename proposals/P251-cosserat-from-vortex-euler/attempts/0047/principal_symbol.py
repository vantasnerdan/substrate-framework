"""Exact frozen-vorticity carrier symbol; not a constant Beltrami solution."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0047-principal-carrier-symbol")
    carrier, lam, vort = sp.symbols("K lambda W", positive=True)
    px, py, pz = sp.symbols("px py pz", real=True)
    omega = sp.Matrix([0, 0, vort])
    hessian = sp.zeros(2)
    kks = sp.zeros(2)
    for sign in (-1, 1):
        wave = sp.Matrix([px, py, sign * carrier + pz])
        polarization = sp.Matrix([1, -sp.I * sign, 0]) / 2
        amplitudes = (polarization, sp.I * sign * polarization)
        generators = [
            -sp.I * wave.cross(amplitude) / carrier for amplitude in amplitudes
        ]
        forces = []
        for generator in generators:
            force = generator.cross(omega)
            forces.append(force.applyfunc(sp.simplify))
        for i in range(2):
            ledger.check(
                f"carrier {sign}: generator {i} is shifted divergence free",
                sp.simplify(wave.dot(generators[i])) == 0,
            )
            for j in range(2):
                # P is orthogonal and curl P=curl: reduce before rational expansion.
                entry = (
                    sp.conjugate(forces[i]).dot(forces[j])
                    - sp.conjugate(wave.dot(forces[i]))
                    * wave.dot(forces[j])
                    / wave.dot(wave)
                    - sp.conjugate(forces[i]).dot(sp.I * wave.cross(forces[j])) / lam
                )
                hessian[i, j] += sp.factor(entry)
                kks[i, j] += omega.dot(sp.conjugate(generators[i]).cross(generators[j]))
    hessian = hessian.applyfunc(sp.factor)
    kks = kks.applyfunc(sp.factor)
    zero = {px: 0, py: 0, pz: 0}
    ledger.check(
        "principal full Hessian recovers compact-cage leading value",
        (hessian.subs(zero) - vort**2 * (1 + carrier / lam) * sp.eye(2)).applyfunc(
            sp.simplify
        )
        == sp.zeros(2),
    )
    ledger.check(
        "principal KKS recovers positive carrier pairing",
        kks.subs(zero) == vort * sp.Matrix([[0, 1], [-1, 0]]),
    )
    ledger.check(
        "H is Hermitian",
        (hessian - hessian.conjugate().T).applyfunc(sp.simplify) == sp.zeros(2),
    )
    ledger.check(
        "Omega is anti-Hermitian",
        (kks + kks.conjugate().T).applyfunc(sp.simplify) == sp.zeros(2),
    )
    symplectic = sp.Matrix([[0, 1], [-1, 0]])
    ratio = pz / carrier
    transform = (sp.eye(2) - sp.I * ratio * symplectic) / (1 - ratio**2)
    normalized_h = (transform.conjugate().T * hessian * transform).applyfunc(
        sp.simplify
    )
    normalized_o = (transform.conjugate().T * kks * transform).applyfunc(sp.simplify)
    ledger.check(
        "exact Darboux map retains varying KKS",
        (normalized_o - vort * symplectic).applyfunc(sp.simplify) == sp.zeros(2),
    )
    a0 = vort**2 * (1 + carrier / lam)
    ax = sp.simplify(sp.diff(normalized_h[0, 0], px, 2).subs(zero) / 2)
    az = sp.simplify(sp.diff(normalized_h[0, 0], pz, 2).subs(zero) / 2)
    drift_h = sp.simplify(sp.diff(normalized_h[0, 1], pz).subs(zero) / sp.I)
    ledger.check(
        "transverse normalized Hessian curvature is negative",
        ax == -(vort**2) / (2 * carrier**2),
    )
    ledger.check("axial diagonal normalized curvature is zero", az == 0)
    ledger.check("chiral drift survives complete action", drift_h == vort**2 / lam)
    inertia0 = sp.simplify(vort**2 / a0)
    inertia2_transverse = sp.simplify(-(vort**2) * ax / a0**2)
    c_transverse = sp.simplify(ax - a0 * inertia2_transverse / inertia0)
    c_axial = sp.simplify(-inertia0 * (drift_h / vort) ** 2)
    ledger.check(
        "complete transverse angle gradient remains negative",
        c_transverse == -(vort**2) / carrier**2,
    )
    ledger.check(
        "parity-averaged axial action gradient is negative",
        c_axial == -(vort**2) / (lam * (carrier + lam)),
    )
    qdot, qz, q = sp.symbols("qdot qz q", real=True)
    speed = vort / lam
    plus = inertia0 * (qdot + speed * qz) ** 2 / 2 - a0 * q**2 / 2
    minus = inertia0 * (qdot - speed * qz) ** 2 / 2 - a0 * q**2 / 2
    average = sp.expand((plus + minus) / 2)
    ledger.check(
        "action averaging retains negative axial stiffness",
        sp.simplify(average - (inertia0 * qdot**2 - a0 * q**2 - c_axial * qz**2) / 2)
        == 0,
    )
    ledger.check(
        "frequency-square averaging would give the opposite drift sign",
        sp.simplify(speed**2 + c_axial / inertia0) == 0 and speed.is_positive,
    )
    print("H_FULL =", hessian)
    print("OMEGA_FULL =", kks)
    print("DARBOUX_H =", normalized_h)
    print("DARBOUX_OMEGA =", normalized_o)
    print("A0 =", a0, "; transverse A2 =", ax, "; drift coefficient =", drift_h)
    print(
        "I0 =",
        inertia0,
        "; C_transverse =",
        c_transverse,
        "; C_axial_after_parity =",
        c_axial,
    )
    print(
        "Scope: frozen local-vorticity principal symbol only. No constant nonzero-lambda Beltrami background is asserted."
    )
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
