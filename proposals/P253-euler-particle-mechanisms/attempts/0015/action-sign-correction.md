# P253/0015 bounded KKS action-sign correction

Date: 2026-09-06

The 0015 critical-loop calculation originally used a KKS form with the sign
opposite to the repository convention. It was internally consistent, but the
shared convention is

\[
 i_{X_H}\Omega=dH,\qquad -d\Theta=\Omega.
\]

For the 0015 tangent
\(\delta_\chi q_j=\{q_j,\chi_j\}\), this requires

\[
 \Omega_q(\delta_\chi q,\delta_\psi q)
 =-\sum_j\int q_j\{\chi_j,\psi_j\}\,dx.
\]

The local action is therefore

\[
 \int\Theta(\dot q)\,dt-H\,dt
 =-\int_\Sigma\Omega-\int H\,dt,
\]

and in a linear symplectic chart its quadratic kinetic term is
\(-\frac12\Omega(q,\dot q)\). Its first variation is
\(-\Omega(\eta,\dot q)-dH(\eta)\), which vanishes at
\(\dot q=X_H\).

Both the KKS form and symplectic-area term were flipped together. The
published pair remains an exact critical loop. The raw-action Hessian remains
strongly indefinite because the two high-time-frequency test families still
give opposite signs; only which family carries which sign is exchanged.
No source scope, autonomous-maximizer verdict, Jacobi equation, monodromy
identity, or remaining dependency changes.
