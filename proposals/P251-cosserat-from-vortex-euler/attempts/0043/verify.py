"""Exact affine Biot--Savart energy; no solver or fitted elastic input."""

import sympy as s

from substrate_framework.euler_affine import affine_vorticity_energy
from substrate_framework.homogenization import sphere_fourth_moment_isotropic

checks = []


def check(name, condition):
    passed = bool(condition)
    checks.append(passed)
    print(f"{'PASS' if passed else 'FAIL'} {name}")


t = s.symbols("t", real=True)
e1, e2, e12, e13, e23 = s.symbols("e1 e2 e12 e13 e23", real=True)
strain = s.Matrix([[e1, e12, e13], [e12, e2, e23],
                   [e13, e23, -e1-e2]])
n = s.Matrix(s.symbols("n1:4", real=True))
sn = (n.T*strain*n)[0]
un = (n.T*strain**2*n)[0]
tr2 = s.trace(strain**2)
# C=F^T F=exp(2t E), C^-1=exp(-2t E). The helicity polarization
# covariance is (Id-n n^T)/2; n.n=1 is imposed before expansion.
numerator = 1-t*sn+t*t*(tr2-un)
denominator = 1-2*t*sn+2*t*t*un
ratio = s.series(numerator/denominator, t, 0, 3).removeO().expand()
check("full covector and vorticity transformation expansion",
      s.expand(ratio-(1+t*sn+t*t*(2*sn**2+tr2-3*un))) == 0)
moment4 = sphere_fourth_moment_isotropic()
average_sn2 = sum(strain[i, j]*strain[k, last]*moment4[i, j, k, last]
                  for i in range(3) for j in range(3)
                  for k in range(3) for last in range(3))
average_sn = s.trace(strain)/3
average_un = s.trace(strain**2)/3
coefficient = s.simplify(2*average_sn2+tr2-3*average_un)
check("isotropic first variation is only pressure", average_sn == 0)
check("full five-dimensional shear coefficient", s.expand(coefficient-4*tr2/15) == 0)
check("positive quadratic form on traceless strains",
      s.expand(tr2-2*(e1**2+e1*e2+e2**2+e12**2+e13**2+e23**2)) == 0)

# Independent exact diagonal affine image of each circular Fourier pair.
f = s.diag(s.exp(t), s.exp(-t), 1)
polarization = s.Matrix([0, 1, s.I])/s.sqrt(2)
wave = s.Matrix([1, 0, 0])
new_wave = f.inv().T*wave
new_omega = f*polarization
check("affine pushforward remains solenoidal", (new_wave.T*new_omega)[0] == 0)
new_velocity = s.I*new_wave.cross(new_omega)/(new_wave.dot(new_wave))
check("curl inverse returns the pushed vorticity",
      s.simplify(s.I*new_wave.cross(new_velocity)-new_omega) == s.zeros(3, 1))
energy_ratio = s.simplify((s.conjugate(new_velocity).T*new_velocity)[0])
check("importable affine energy equals independent curl inversion",
      s.simplify(2*affine_vorticity_energy(f, {tuple(wave): polarization}, 1)-energy_ratio) == 0)
check("direct Fourier energy agrees with helicity covariance",
      s.simplify(energy_ratio-(1+s.exp(2*t))/2) == 0)
check("fixed-wavevector mutation changes energy",
      s.simplify(energy_ratio-(1+s.exp(-2*t))/2) != 0)
angle = s.symbols("angle", real=True)
rotation = s.Matrix([[s.cos(angle), -s.sin(angle), 0],
                     [s.sin(angle), s.cos(angle), 0], [0, 0, 1]])
rot_wave, rot_omega = rotation*wave, rotation*polarization
check("common rigid rotation has exactly zero energy increment",
      s.simplify((s.conjugate(rot_omega).T*rot_omega)[0]/rot_wave.dot(rot_wave)-1) == 0)

rho, a, b, energy0, mu = s.symbols("rho a b energy0 mu", positive=True)
actual_energy = rho*(a*a+b*b)/2
matched_mu = s.solve(s.Eq(mu*tr2, energy0*coefficient), mu)[0]
check("shear modulus derived by coefficient matching", matched_mu == 4*energy0/15)
check("same actual two-mode tube gives nonzero shear",
      s.simplify(matched_mu.subs(energy0, actual_energy)-2*rho*(a*a+b*b)/15) == 0)
check("structure-free amplitude limit removes this shear",
      s.limit(s.limit(matched_mu.subs(energy0, actual_energy), a, 0), b, 0) == 0)
print(f"{sum(checks)}/{len(checks)} checks passed")
raise SystemExit(0 if all(checks) else 1)
