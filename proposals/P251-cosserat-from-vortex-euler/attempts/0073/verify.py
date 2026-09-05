"""Exact parcel tensor bookkeeping; not a microscopic Euler model."""

import sympy as s

from substrate_framework.verification import CheckLedger

ledger = CheckLedger("P251/0073")
mass, radius = s.symbols("m a", positive=True)
velocity = s.Matrix(s.symbols("v0:3", real=True))
spin_rate = s.Matrix(s.symbols("w0:3", real=True))
strain = s.Matrix([[s.Symbol("e0"), s.Symbol("e1"), s.Symbol("e2")],
                   [s.Symbol("e1"), s.Symbol("e3"), s.Symbol("e4")],
                   [s.Symbol("e2"), s.Symbol("e4"), s.Symbol("e5")]])
positions = [sign*radius*s.eye(3)[:, axis] for axis in range(3) for sign in (-1, 1)]
internal = [spin_rate.cross(r)+strain*r for r in positions]
q = sum((mass*w*r.T/6 for r, w in zip(positions, internal, strict=True)), s.zeros(3))
shape_rate = sum((mass*(r*w.T+w*r.T)/6
                  for r, w in zip(positions, internal, strict=True)), s.zeros(3))
angular = sum((mass*r.cross(w)/6 for r, w in zip(positions, internal, strict=True)),
              s.zeros(3, 1))
antisymmetric = s.Matrix(3, 3, lambda i, j:
    -sum(s.LeviCivita(i, j, k)*angular[k] for k in range(3))/2)
ledger.check("exact physical spin separates from symmetric shape rate",
             s.simplify(angular-2*mass*radius**2*spin_rate/3) == s.zeros(3, 1))
ledger.check("complete momentum dipole includes shape and half the spin",
             s.simplify(q-shape_rate/2-antisymmetric) == s.zeros(3))
k = s.Matrix(s.symbols("k0:3", real=True))
first = sum((mass*(velocity+w)*(1-s.I*k.dot(r))/6
             for r, w in zip(positions, internal, strict=True)), s.zeros(3, 1))
target = mass*velocity-s.I*shape_rate*k/2+s.I*k.cross(angular)/2
ledger.check("Eulerian first Fourier moment equals centre plus shape and spin currents",
             s.simplify(first-target) == s.zeros(3, 1))
ledger.check("using a full rather than half curl spin changes the actual momentum",
             s.simplify(first-(mass*velocity-s.I*shape_rate*k/2+s.I*k.cross(angular)))
             != s.zeros(3, 1))
energy = sum(mass*(velocity+w).dot(velocity+w)/12 for w in internal)
split = mass*velocity.dot(velocity)/2+sum(mass*w.dot(w)/12 for w in internal)
ledger.check("material centring gives the exact kinetic split",
             s.expand(energy-split) == 0)
raise SystemExit(ledger.finish())
