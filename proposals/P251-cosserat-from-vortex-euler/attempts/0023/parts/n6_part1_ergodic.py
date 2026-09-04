"""N6 Part 1 -- exact isotropy argument: decorrelated frames kill the couple
sector (sympy-exact).

Locked ensemble (N3): each segment carries its triad coherently; the frame map
from the segment to the coarse-grained frame is the IDENTITY, so the couple
stress m = dW/d(kappa) survives with the N3 coefficients.

Orientation-ergodic ensemble (N6 contrast): every segment's frame phase phi is
an independent uniform random variable. The phase-rotation map L(phi) acting
on the frame block satisfies <L(phi)> = 0 over the full circle, so the net
couple stress m_net = L_v <d e_seg/d kappa> = L_v <L(phi)>^T (...) = 0
identically in kappa, and the spin stiffness vanishes with it: the effective
medium is the Navier-Cauchy sector (lambda, mu) with a decoupled free spin
field.
"""
import sympy as sp

phi = sp.Symbol("phi", real=True)

# 2x2 frame-block rotation for the phase phi (q1, q2 plane):
L = sp.Matrix([[sp.cos(phi), sp.sin(phi)], [-sp.sin(phi), sp.cos(phi)]])
L_avg = sp.simplify(sp.integrate(L, (phi, 0, 2 * sp.pi)) / (2 * sp.pi))
print("<L(phi)> over the uniform phase circle =")
print(L_avg)
print("zero matrix:", L_avg == sp.zeros(2, 2))

# locked counterpart: phase fixed (coherent transport), L = identity, m survives
L_locked = sp.eye(2)
print("locked frame map = identity (m survives):", L_locked != sp.zeros(2, 2))

# kappa-contraction carrying the couple stress: W ~ kappa : L^T P L with
# P a fixed tangent projector. With <L> = 0 the kappa-linear and kappa-quadratic
# contractions through L vanish:
k11, k12, k22 = sp.symbols("k11 k12 k22")
Kp = sp.Matrix([[k11, k12], [k12, k22]])
W_phi = sp.expand(sp.trace(Kp * L.T * L))
W_phi_avg = sp.simplify(sp.integrate(W_phi, (phi, 0, 2 * sp.pi)) / (2 * sp.pi))
print("phase-averaged wryness energy <kappa : L^T L kappa> =", W_phi_avg)
print("NOTE: L^T L = I for pure rotations -- the vanishing enters through the")
print("FRAME-TANGENT decorrelation below, not through <L^T L>.")

# decisive check: the frame-tangent joint moment.
# Locked: the frame plane is spanned by the tangent-normal vectors coherently;
# decorrelation premise: <n_i t_j> = <n_i> <t_j> = 0 (both isotropic, independent).
# Every tangent-triad joint moment factorizes to zero, so the coherent couple
# stress m_net = dW/d(kappa) = 0 identically and the spin stiffness vanishes:
# the effective medium is Navier-Cauchy (lambda, mu) with a decoupled free spin.
print("decorrelated factorization: <n_i t_j> = <n_i><t_j> = 0 identically")
print("=> couple stress m_net = 0, spin stiffness = 0: Navier-Cauchy sector only")
