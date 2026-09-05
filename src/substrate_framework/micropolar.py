"""Conditional, unpromoted isotropic micropolar action (P251/0056).

This module encodes the action and its Fourier operator, not the microscopic
Euler construction or an accepted claim. Coefficients must come from that
construction when the caller makes a microscopic claim. Incompressible Euler
uses div U=0 and a pressure multiplier; the optional Lame coefficient describes
only the formal compressible extension. Rotation itself need not be solenoidal.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class MicropolarCoefficients:
    """W=lam tr(e)²/2+mu||e||²+alpha|curl U-2Phi|²/2+W_curvature.

    W_curvature=c_tr tr(G)²+c_s||sym G||²+c_a||skew G||².
    Positivity is a caller hypothesis: mu,alpha,c_s,c_a>0,
    3*c_tr+c_s>0 (and the appropriate bulk condition if compressible).
    Signs are deliberately not enforced, permitting counterexample probes.
    """

    lame: sp.Expr
    shear: sp.Expr
    locking: sp.Expr
    trace_curvature: sp.Expr
    symmetric_curvature: sp.Expr
    skew_curvature: sp.Expr

    @property
    def transverse_curvature(self):
        return sp.sympify(self.symmetric_curvature)+self.skew_curvature

    @property
    def longitudinal_curvature(self):
        return 2*(sp.sympify(self.symmetric_curvature)+self.trace_curvature)


def _matrix(value, shape):
    value = sp.Matrix(value)
    if value.shape != shape:
        raise ValueError(f"expected matrix shape {shape}")
    return value


def isotropic_micropolar_energy(displacement_gradient, rotation, rotation_gradient,
                                coefficients):
    """Return the local quadratic energy, with h_ij=partial_j U_i.

    Real fields and coefficients are intended. All gradients and rotations
    remain independent arguments so stress, couple stress and local torque
    follow by direct differentiation. Boundary terms have not been discarded.
    """
    h = _matrix(displacement_gradient, (3, 3))
    phi = _matrix(rotation, (3, 1))
    g = _matrix(rotation_gradient, (3, 3))
    c = coefficients
    strain = (h+h.T)/2
    symmetric, skew = (g+g.T)/2, (g-g.T)/2
    curl = sp.Matrix([h[2, 1]-h[1, 2], h[0, 2]-h[2, 0], h[1, 0]-h[0, 1]])
    relative = curl-2*phi
    return (sp.sympify(c.lame)*sp.trace(h)**2/2
            +c.shear*sum(value**2 for value in strain)
            +c.locking*relative.dot(relative)/2
            +c.trace_curvature*sp.trace(g)**2
            +c.symmetric_curvature*sum(value**2 for value in symmetric)
            +c.skew_curvature*sum(value**2 for value in skew))


def micropolar_fourier_stiffness(wave_vector, coefficients):
    """Hermitian 6x6 K for exp(i k.x), ordered (U_x,U_y,U_z,Phi_x,Phi_y,Phi_z).

    The equation is (K-omega²*diag(rho I,j I))*amplitude=0. Both off-diagonal
    blocks are -2*alpha*curl(k), and curl(k)=i*[k cross] is Hermitian for
    real k. The Euler displacement sector is restricted to k.U=0 by pressure;
    this leaves longitudinal spin, whose curvature is 2*(c_s+c_tr).
    """
    k = _matrix(wave_vector, (3, 1))
    c = coefficients
    square, longitudinal = k.dot(k), k*k.T
    cross = sp.Matrix([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    curl = sp.I*cross
    u = ((c.shear+c.locking)*square*sp.eye(3)
         +(c.lame+c.shear-c.locking)*longitudinal)
    phi = ((c.transverse_curvature*square+4*c.locking)*sp.eye(3)
           +(2*c.trace_curvature+c.symmetric_curvature-c.skew_curvature)*longitudinal)
    coupling = -2*c.locking*curl
    return sp.ImmutableMatrix(sp.BlockMatrix([[u, coupling], [coupling, phi]]))


def uniform_phase_average(expression, phases):
    """Exact product-uniform phase integral, one normalized circle per symbol.

    The input may be any integrable SymPy expression. Invariance under frame
    shifts additionally requires periodicity in every phase. Independence is
    encoded by the product measure; correlated phases must instead be pulled
    back to their actual common variables before calling this function.
    Unresolved symbolic integrals remain explicit rather than becoming numbers.
    """
    phases = tuple(phases)
    if any(not isinstance(phase, sp.Symbol) for phase in phases) or len(set(phases)) != len(phases):
        raise ValueError("phases must be distinct symbols")
    value = sp.sympify(expression)
    for phase in phases:
        value = sp.integrate(value, (phase, 0, 2*sp.pi))/(2*sp.pi)
    return sp.simplify(value)
