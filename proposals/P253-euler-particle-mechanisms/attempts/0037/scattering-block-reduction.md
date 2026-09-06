# One-way column scattering and the radial-complement scale

The supercritical speed has a second consequence beyond the scalar limit:
all linear column radiation travels backward relative to the solitary wave.
This section derives the exact free estimate and isolates what must still be
proved for the variable solitary operator.

## 1. Exact branchwise local-energy estimate

Let `omega_n(k)=k c_n(k)` be a positive-frequency axisymmetric column branch
in the physical metric. In the frame moving with the exact solitary speed
`c>c0`, put

    nu_n(k)=omega_n(k)-c k.

The exact group-speed bound in `0030` gives

    nu_n'(k)=v_g,n(k)-c <=-(c-c0)<0.                 (1)

For a spectrally smooth packet on this branch,

    u(z,t)=integral exp[i k z-i nu_n(k)t] a(k) dk,

Plancherel in `t` and the change of variable `k -> nu_n(k)` give, at every
fixed `z`,

    integral_R |u(z,t)|^2 dt
      =2pi integral |a(k)|^2/|nu_n'(k)| dk
      <=2pi ||a||_2^2/(c-c0).                       (2)

Consequently a spatial window of length `L` satisfies

    integral_R ||1_window u(t)||^2 dt
      <=C L ||u(0)||^2/(c-c0).                      (3)

This is an exact one-pass estimate. It is not a dispersive pointwise guess.
It also displays the critical loss `(c-c0)^-1=O(mu^-1)`.

At `k=0`, the finite-interval radial form has a simple critical eigenvector
and a positive gap on its radial orthogonal complement. Continuity in `k`
and the large-`k` form bound give a number `delta_2>0` such that every
noncritical radial branch obeys

    c-v_g,n(k)>=delta_2                              (4)

for sufficiently small `mu`, after the stationary kernel is separated.
On that complement, (3) improves to `C L/delta_2`.

The branch decomposition must be implemented as the measurable spectral
resolution of the positive column oscillator; a formal sum of radial modes
is not used to claim (3) on the full energy space. Stationary kernel data are
transported by `-c partial_z` and satisfy the same one-way estimate directly.

## 2. Why the radial complement is perturbative

The exact `0027` profile bounds imply

    ||partial_r psi_s||_infinity=O(mu),
    ||partial_z psi_s||_infinity=O(mu/L_mu),          (5)

and the same scale for the smooth label coefficients on their fixed radial
support. Since

    mu L_mu=f0(R)^2 log(L_mu)/L_mu ->0,              (6)

a noncritical packet crossing the length-`L_mu` solitary region accumulates
only `O(mu L_mu)` radial-strain and cross-channel coupling. Terms carrying an
axial derivative of a coefficient are smaller. Thus the complement is not a
second slow system.

The critical branch has relative speed `O(mu)`, remains in the wave for time
`O(L_mu/mu)`, and accumulates an order-one effect. That effect is exactly the
KdV block in `threshold-reduction.md`; treating it perturbatively would lose
the required bound.

## 3. Block theorem still required

Let `Pi_mu` denote the transported critical spectral projection of the column
oscillator and `Q_mu=1-Pi_mu`. The full proof earns the desired linear result
by constructing a bounded near-identity graph

    q_Q=G_mu q_P+r_Q,   ||G_mu||=O(mu L_mu),          (7)

for the actual solitary generator, such that:

1. the `q_P` generator converges in the graph norm to the KdV linearization
   `-sigma partial_X L_*` with the translation pole and the physical
   Casimir/momentum rows retained;
2. `r_Q` has the one-pass estimate obtained from (2)--(4), uniformly for all
   time;
3. the off-diagonal Duhamel products are `O(mu L_mu)` rather than merely
   `O(1)`; and
4. the Bessel exterior and `k=0` stationary kernel remain in the graph.

The required cancellation in item 3 is structural: the source entering the
complement is an axial derivative after the mixed-Casimir row is imposed.
Without proving that row, estimate (3) alone carries the critical
`mu^-1` loss and does not close.

This reduces Route A to a concrete operator identity and estimate: construct
the exact `Pi/Q` symplectic projection and prove the derivative-form
off-diagonal factor in the physical metric. If it holds, (2), (6), and scalar
KdV coercivity give the linear all-time scattering bound. The nonlinear
bootstrap then requires the corresponding quadratic local-energy estimate
and conditional classical existence. If the derivative cancellation fails,
its surviving coefficient becomes Route C's explicit Evans coupling rather
than a generic obstruction.
