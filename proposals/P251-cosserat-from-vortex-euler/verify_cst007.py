"""C-CST-007 verifier: EPS existence import and Beltrami regression.

The Enciso--Peralta-Salas existence statements are imported from hash-pinned
primary sources; this script does not re-prove them. The explicit periodic
field

    u = (sin(lambda z), cos(lambda z), 0)

is a separate exact regression for the Beltrami-to-steady-Euler implication.
It has constant speed and planar field lines, so it is not presented as a
localized or knotted vortex tube. The vortex-tube existence role belongs to the
audited EPS theorem statements, not to this elementary field.
"""

import hashlib
import sys
from pathlib import Path

import sympy as sp

from substrate_framework.verification import CheckLedger

x, y, z = sp.symbols("x y z", real=True)
lam = sp.Symbol("lam", positive=True)
u = sp.Matrix([sp.sin(lam * z), sp.cos(lam * z), 0])


def curl(field: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [
            sp.diff(field[2], y) - sp.diff(field[1], z),
            sp.diff(field[0], z) - sp.diff(field[2], x),
            sp.diff(field[1], x) - sp.diff(field[0], y),
        ]
    )


def check_beltrami_regression(ledger: CheckLedger) -> None:
    divergence = sp.simplify(
        sum(sp.diff(u[index], (x, y, z)[index]) for index in range(3))
    )
    vorticity = curl(u)
    helicity = sp.simplify(u.dot(vorticity))
    norm2 = sp.simplify(u.dot(u))
    ledger.check("periodic plane wave is incompressible", divergence == 0, "div u = 0")
    ledger.check(
        "periodic plane wave is a Beltrami eigenfield",
        all(sp.simplify(vorticity[index] - lam * u[index]) == 0 for index in range(3)),
        "curl u = lambda u",
    )
    ledger.check(
        "helicity density is lambda times speed squared",
        sp.simplify(helicity - lam * norm2) == 0,
        "u.curl(u) = lambda |u|^2",
    )

    convective = sp.Matrix(
        [
            sum(u[jj] * sp.diff(u[ii], (x, y, z)[jj]) for jj in range(3))
            for ii in range(3)
        ]
    )
    ledger.check(
        "periodic plane wave solves stationary constant-pressure Euler",
        all(sp.simplify(value) == 0 for value in convective),
        "(u.grad)u = 0",
    )
    ledger.check(
        "regression field is explicitly nonlocalized and is not used as the EPS tube",
        sp.simplify(norm2 - 1) == 0 and u[2] == 0,
        "|u| is constant and dz/ds=0 along its field lines",
    )


def check_source_integrity_and_scope(ledger: CheckLedger) -> None:
    base = Path(__file__).parent / "sources"
    expected = {
        "1003.3122.pdf": "a179105fd823baecc71b10b37cab5a3f",
        "1210.6271.pdf": "6349631cfdfe0d71a4673340f1056f29",
        "1505.01605.pdf": "af8e5f917bd5dc1f52642426f6de48f1",
    }
    actual = {
        name: hashlib.md5((base / name).read_bytes()).hexdigest()
        for name in expected
        if (base / name).is_file()
    }
    ledger.check(
        "all three EPS primary PDFs match their recorded digests",
        actual == expected,
        f"matched={sorted(actual)}",
    )

    tube_text = (base / "1210.6271.pdf.txt").read_text(encoding="utf-8")
    torus_text = (base / "1505.01605.pdf.txt").read_text(encoding="utf-8")
    ledger.check(
        "archived EPS text exposes the cited vortex-tube theorem statements",
        "Theorem 1.1." in tube_text
        and "vortex tubes" in tube_text
        and "Theorem 1.1." in torus_text
        and "vortex tubes" in torus_text
        and "curl u" in torus_text,
        "statement-location audit only; the theorem is a declared import",
    )


def check_mutations(ledger: CheckLedger) -> None:
    u_bad = sp.Matrix([sp.sin(lam * z), sp.cos(lam * z), sp.Symbol("eps") * x])
    curl_bad = curl(u_bad)
    ledger.check(
        "M1 non-Beltrami perturbation rejected",
        not all(
            sp.simplify(curl_bad[index] - lam * u_bad[index]) == 0 for index in range(3)
        ),
        "curl u_bad != lambda u_bad",
    )

    vorticity = curl(u)
    ledger.check(
        "M2 wrong-sign helicity rejected",
        sp.simplify(u.dot(vorticity) + lam * u.dot(u)) != 0,
        "the exact sign is positive for lambda > 0",
    )

    source = (Path(__file__).parent / "sources" / "1210.6271.pdf").read_bytes()
    tampered = bytearray(source)
    tampered[0] ^= 1
    ledger.check(
        "M3 one-byte source mutation changes the digest",
        hashlib.md5(tampered).hexdigest() != hashlib.md5(source).hexdigest(),
        "digest sensitivity is exercised on the archived bytes",
    )


def main() -> int:
    ledger = CheckLedger("C-CST-007")
    check_beltrami_regression(ledger)
    check_source_integrity_and_scope(ledger)
    check_mutations(ledger)
    return int(ledger.finish())


if __name__ == "__main__":
    sys.exit(main())
