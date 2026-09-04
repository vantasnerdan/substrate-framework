"""Exact one-action field map, helix length variation and modal-order check."""
import sympy as s

from substrate_framework.homogenization import sphere_fourth_moment_isotropic
from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0034-longwave")
    I, K = s.symbols("I K", positive=True)  # noqa: E741 -- physical inertia symbol
    pd, bd, sd, psi, beta = s.symbols("phi_dot beta_dot psi_dot psi beta", real=True)
    kinetic = s.Rational(11, 20)*I*(pd**2+bd**2)-s.Rational(9, 10)*I*pd*bd
    mapped = s.expand(kinetic.subs(pd, (2*sd+9*bd)/11))
    Jpsi, Jbeta = 2*I/55, 4*I/11
    Kpsi = 4*K/121
    ledger.check("covariant angle map diagonalizes the complete kinetic form",
                 s.simplify(mapped-(Jpsi*sd**2+Jbeta*bd**2)/2) == 0)
    ledger.check("same angle map fixes the elastic coefficient",
                 s.simplify(K*((2*psi+9*beta)/11-beta)**2/2
                            -Kpsi*(psi-beta)**2/2) == 0)
    common = s.Symbol("common_angle")
    fieldmap = (11*s.Symbol("phi")-9*beta)/2
    ledger.check("collective microrotation has unit common-rotation weight",
                 s.expand(fieldmap.subs({s.Symbol("phi"): s.Symbol("phi")+common,
                                        beta: beta+common}, simultaneous=True)-fieldmap) == common)
    T, radius, eps = s.symbols("T radius epsilon", positive=True)
    ps, bs, ss = s.symbols("phi_s beta_s psi_s", real=True)
    line_energy = 3*T*(s.sqrt(1+radius**2*(eps*ps)**2)
                       +s.sqrt(1+radius**2*(eps*bs)**2)-2)
    second = s.diff(line_energy, eps, 2).subs(eps, 0)/2
    ledger.check("phase-gradient energy follows from actual helical line lengths",
                 s.simplify(second-3*T*radius**2*(ps**2+bs**2)/2) == 0)
    mapped_gradient = s.expand(second.subs(ps, (2*ss+9*bs)/11))
    Cpsi = 12*T*radius**2/121
    ledger.check("derived spin gradient coefficient is Cpsi=12*T*r^2/121",
                 s.simplify(mapped_gradient.coeff(ss, 2)-Cpsi/2) == 0)
    ledger.check("macro-frame derivatives remain in the complete gradient energy",
                 s.simplify(mapped_gradient-Cpsi*ss**2/2) != 0)

    # Density counts bundle axes; original filament density is 6*L_cell.
    Lv = s.Symbol("L_cell", positive=True)
    kp = s.Matrix(3, 3, lambda i, j: s.Symbol(f"kappa{i}{j}"))
    fourth = sphere_fourth_moment_isotropic()
    mean = sum(kp[i, j]*kp[ell, n]*fourth[i, j, ell, n]
               for i in range(3) for j in range(3)
               for ell in range(3) for n in range(3))
    sym = (kp+kp.T)/2
    target = Lv*Cpsi*s.trace(kp)**2/30+Lv*Cpsi*sum(v**2 for v in sym)/15
    ledger.check("isotropic fourth moment fixes all couple coefficients",
                 s.simplify(Lv*Cpsi*mean/2-target) == 0)
    alpha, j = Lv*Kpsi/12, Lv*Jpsi/3
    ledger.check("locking and inertia come from one action and one field map",
                 s.simplify(4*alpha/j-Kpsi/Jpsi) == 0)

    k, rho, mu, al, jj, c, g, A, B = s.symbols("k rho mu alpha j c g A B", positive=True)
    z = s.Symbol("omega_squared")
    determinant = ((rho+g*k**2)*z-(mu+al)*k**2)*(jj*z-4*al-c*k**2)-4*al**2*k**2
    acoustic = s.expand(determinant.subs(z, A*k**2)).coeff(k, 2)
    optical = s.expand(determinant.subs(z, 4*al/jj+B*k**2)).coeff(k, 2)
    a_coeff = s.solve(acoustic, A)[0]
    b_coeff = s.solve(optical, B)[0]
    ledger.check("macro gradient inertia leaves acoustic omega^2 through k^2 unchanged",
                 s.simplify(a_coeff-mu/rho) == 0)
    ledger.check("macro gradient inertia leaves optical omega^2 through k^2 unchanged",
                 s.simplify(b_coeff-c/jj-al/rho) == 0)
    ledger.check("finite-k characteristic equation retains the gradient-inertia correction",
                 s.simplify(determinant-determinant.subs(g, 0)) != 0)
    print("Jpsi =", Jpsi, "; Kpsi =", Kpsi, "; Cpsi =", Cpsi)
    print("alpha =", alpha, "; j =", j, "; optical k^2 coefficient =", b_coeff)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())
