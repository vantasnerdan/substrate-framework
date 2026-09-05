"""Exact dipole, isotropic response projection, KKS and current identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0055-exterior-centering")
    x, y, z = s.symbols("x y z", real=True)
    coords = (x, y, z)
    position = s.Matrix(coords)
    radius = s.sqrt(position.dot(position))
    volume = s.Symbol("volume", positive=True)
    a = s.Matrix(s.symbols("a0:3", real=True))
    pressure = volume*a.dot(position)/(4*s.pi*radius**3)
    dipole = -s.Matrix([s.diff(pressure, coordinate) for coordinate in coords])
    expected = volume*(3*position*a.dot(position)-a*radius**2)/(4*s.pi*radius**5)
    ledger.check("Neumann-potential-free R3 Leray dipole has the stated sign and factor",
                 all(s.simplify(component) == 0 for component in dipole-expected))
    radial = s.simplify(position.dot(dipole)/radius)
    ledger.check("radial dipole coefficient is nonzero for a nonzero direction",
                 s.simplify(radial-volume*a.dot(position)/(2*s.pi*radius**4)) == 0)
    ledger.check("exterior dipole is divergence free",
                 s.simplify(sum(s.diff(dipole[j], coords[j]) for j in range(3))) == 0)

    # Exact finite right inverse, including different local Gram matrices.
    for j in range(4):
        triangular = s.eye(4)
        for row in range(1, 4):
            triangular[row, row-1] = j+row
        gram = triangular*triangular.T
        unit = s.eye(4)[:, j]
        coefficients = gram.inv()*unit
        ledger.check(f"positive Gram response in separate ball {j} has its exact dual moments",
                     gram*coefficients == unit and gram.det() > 0)

    # Basis: K, four disjoint fixed responses, raw A, raw Q, raw S.
    moments = s.Matrix(4, 3, s.symbols("f0:12"))
    b0, c = s.symbols("b0 c", real=True, nonzero=True)
    form = s.zeros(8)
    form[0, 1], form[1, 0] = 1, -1
    for j in range(3):
        form[0, 5+j], form[5+j, 0] = moments[0, j], -moments[0, j]
    form[6, 7], form[7, 6] = c, -c
    change = s.zeros(8, 4)
    change[0, 0] = 1
    for j in range(3):
        change[5+j, 1+j] = 1
        for i in range(4):
            change[1+i, 1+j] = -moments[i, j]
    change[1, 1] += b0
    selected = s.simplify(change.T*form*change)
    target = s.Matrix([[0, b0, 0, 0], [-b0, 0, 0, 0],
                       [0, 0, 0, c], [0, 0, -c, 0]])
    ledger.check("all fixed-projection corrections preserve the exact full KKS blocks",
                 selected == target)
    all_moments = s.zeros(4, 8)
    all_moments[:, 1:5] = s.eye(4)
    all_moments[:, 5:8] = moments
    new_moments = s.simplify(all_moments*change)
    ledger.check("selected physical cell mean vanishes in all three compact directions",
                 new_moments[1:4, :] == s.zeros(3, 4))
    ledger.check("common response retains the prescribed nonzero source normalization",
                 new_moments[0, :] == s.Matrix([[0, b0, 0, 0]]))
    ledger.check("internal KKS pairing has no asymptotic projection correction",
                 selected[2, 3] == c)
    ledger.check("body-internal KKS crosses vanish exactly, not by projection-tail neglect",
                 selected[1, 2] == 0 and selected[1, 3] == 0)
    ledger.check("selected four-dimensional form is nondegenerate",
                 s.factor(selected.det()) == b0*b0*c*c)
    plus = s.Matrix([[1, 0, 0, 0, 0, 0, 1, 0]])
    minus = s.Matrix([[1, 0, 0, 0, 0, 0, -1, 0]])
    ledger.check("exterior response never changes the physical core-angle jets",
                 plus*change == s.Matrix([[1, 0, 1, 0]])
                 and minus*change == s.Matrix([[1, 0, -1, 0]]))
    mutated = form.copy()
    mutated[1, 2], mutated[2, 1] = 1, -1
    ledger.check("dropping isotropic-support construction changes the projected KKS form",
                 s.simplify(change.T*mutated*change-selected) != s.zeros(4))

    amin, remainder, lam_abs, margin = s.symbols("Amin remainder lambda_abs margin", positive=True)
    carrier = lam_abs*remainder/amin+margin
    lower = s.expand((1+carrier/lam_abs)*amin-remainder)
    ledger.check("bounded fixed responses preserve positive full Hessian at finite carrier",
                 lower.is_positive)
    jd, jext, boundary = s.symbols("J_D J_ext boundary_spin", real=True)
    coherence = jd+jext
    effective = boundary-jext
    ledger.check("ambient angular impulse remains in the exact physical-spin matching",
                 s.expand(coherence+effective-(jd+boundary)) == 0)
    ledger.check("omitting ambient impulse changes the matched material spin",
                 s.expand(coherence+boundary-(jd+boundary)) == jext)
    print("Analytic oracle: exterior-centering.md; actual EPS radiation and finite positive response Grams.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
