# The actual multi-coordinate observation chart and its second-order action

Let z_dot=B(t)z on a finite actual phase family with constant invertible
skew Omega and B^T Omega+Omega B=0. The Hamiltonian is H=-Omega B.
Let q=C(t)z be the registered n physical positions in a2n-dimensional
family. No position is replaced by a convenient canonical coordinate.
Define D=C_dot+CB and T=[C;D]. Where det T is nonzero, y=(q,v)=Tz
is the actual configuration/rate chart, with v=q_dot and

    T_dot=[C_dot; C_ddot+C_dot B+C B_dot],
    E=T^-1, E_dot=-E T_dot E,
    Omega_y=E^T Omega E,
    H_y=E^T H E+sym(E^T Omega E_dot),
    B_y=(T_dot+TB)E.

These are precisely the complete moving_phase_pullback identities at a
full-rank embedding; the ambient residual is zero. Both Omega_y_dot and
the symmetric connection in H_y remain. The physical configuration
Poisson bracket is -C Omega^-1 C^T. Its vanishing at one time does not
prove vanishing at other times.

## Ordinary mechanical action, when the physical positions license it

If the lower-right block of Omega_y is zero THROUGHOUT the window,
write Omega_y=[[A,M],[-M^T,0]]. The upper row of B_y is [0,I] because
v is the actual derivative of q. Differentiating symplectic preservation
gives Omega_y_dot+B_y^T Omega_y+Omega_y B_y=0. Its lower-right block
then gives M=M^T. The full Hamiltonian necessarily has the form

    H_y=[[K,M_dot/2],[M_dot/2,M]],  K=K^T, A=-A^T.

Consequently the first-order action, after an explicit total derivative
and elimination of v, is exactly

    L(q,q_dot)=q_dot^T M q_dot/2-q^T A q_dot/2-q^T K q/2.

Its canonical momentum and Euler-Lagrange equation are

    p_can=M q_dot+Aq/2,
    M q_ddot+(M_dot+A)q_dot+(K+A_dot/2)q=0.

The actual measured momentum row P(t)z is instead PE y; its difference
from [A/2,M]y is returned explicitly. Equality, positivity of M, and
time-independent constitutive coefficients are separate physical facts.

The implementation returns an exact ordinary-action condition built from
the lower-right symplectic block, its time derivative, M-M^T, and the
displayed Hamiltonian block identities. Undecidable symbolic equalities
remain conditions, not a false proof. A nonzero lower-right block is
retained in the full phase action; it cannot be discarded to manufacture
an ordinary q/q_dot Lagrangian. In particular configurations may commute
at t=0 while the time derivative of that bracket is nonzero: the ordinary
action condition then exposes the missing all-time license even at t=0.

## Oracle and parent role

Direct variation of a two-coordinate moving configuration map is the
strongest exact test. Independent tests compare the scalar limit, a
noncommuting configuration, an initially commuting but later noncommuting
map, and measured versus canonical momentum. The computation is an exact
finite phase identity.0158 still supplies actual Euler columns, physical
moment conditions and estimates; this identity alone supplies none of them.
