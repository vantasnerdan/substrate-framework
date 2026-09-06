# Exact Schwinger--Hopf algebra and its Euler realization test

## 1. Two positive canonical modes give the algebra exactly

Assume first that a physical Euler carrier has supplied two positive
canonical pairs `(q_a,p_a)`, `a=1,2`, with identical KKS normalization

    Omega = B sum_a dq_a wedge dp_a,       B>0.          (1)

Define

    z_a=sqrt(B/2)(q_a+i p_a),
    {z_a,conj(z_b)}=-i delta_ab,            (2)

using the physical bracket

    {F,G}=B^-1 sum_a(F_qa G_pa-F_pa G_qa). (3)

Put

    I=(|z_1|^2+|z_2|^2)/2,
    S_i=z^dagger sigma_i z/2.               (4)

Direct differentiation gives

    {S_i,S_j}=epsilon_ijk S_k,
    {I,S_i}=0,
    S_x^2+S_y^2+S_z^2=I^2.                 (5)

For `I=I_0>0`, the level set is `S3`.  The common phase
`z -> exp(i alpha)z` is a free `U(1)` action, and its quotient is
`CP1=S2`.  The reduced KKS sphere has radius `I_0` in (5) and area

    integral_S2 Omega_red=4 pi I_0.         (6)

Thus the Schwinger representation supplies the compact `su(2)` algebra that
one noncompact plane lacked in `0046`.  It also exposes the load-bearing
premise: two independent complex modes are required.  The two real
quadratures of one oscillator provide only one `z`, not `z in C2`.

If the autonomous quadratic Hamiltonian is

    H_2=omega_1 |z_1|^2+omega_2 |z_2|^2,   (7)

then `I` is conserved even when the frequencies differ.  Exact degeneracy
`omega_1=omega_2` makes (7) proportional to `I` and leaves every Stokes vector
fixed under free evolution.  A Hermitian mixing Hamiltonian
`H_h=z^dagger h z` produces `zdot=-i h z`, preserves `I`, and rotates `S`.
But each off-diagonal entry of `h` is a physical cross-mode interaction which
must be derived from Euler.  Writing down `U in U(2)` is only a coordinate or
control prescription until that interaction and its complement are built.

## 2. Route A: the reviewed single-carrier inputs do not yet contain a doublet

The exact `0030/0035` column operator decomposes into simple radial
eigenvalues for each nonzero axial wave number.  A real field pairs the
`+k` and `-k` coefficients as complex conjugates.  Its axial cosine and sine
profiles are the canonical quadratures of the same oscillator.  Counting
them as `z_1,z_2` would double the physical degrees of freedom.  Distinct
radial indices give two genuine modes, but simplicity supplies no exact equal
frequency and no Euler operation that mixes them while preserving a closed
four-real-dimensional subspace.

The ring bending pair has the same issue: cosine and sine around one azimuthal
harmonic form one traveling/standing-wave oscillator.  The internal Kelvin
modes appearing in `0048` are distinct candidates, but the full fixed-domain
Riesz projection, all-sector complement, and nonlinear invariant carrier are
active `0052` obligations.  They cannot be used here to assert a physical
doublet.

**Route A verdict:** the exact algebraic sufficiency criterion is established;
the current same-carrier transfer is blocked by a second independent positive
mode, a symmetry-protected degeneracy or controlled splitting, and a physical
off-diagonal mixing interaction with full Euler leakage control.

## 3. Route B: two disjoint copies give a product, not a one-carrier analyzer

Two identical noninteracting compact Euler cells would provide two copied KKS
planes and hence the algebra (1)--(6).  At strict support separation their
direct-sum Hamiltonian has the required degeneracy.  Yet the two complex
amplitudes then label two spatial carriers.  The quotient is a collective
relative-amplitude sphere, not an internal state of one carrier.

More decisively, disjoint noninteracting cells supply only independent phase
rotations.  Off-diagonal `U(2)` mixing requires a coupling which transfers
action between the cells.  Turning that coupling on must retain a single
Euler field, the pressure response, support/collision control, and leakage
into deformation modes.  No accepted compact-pair theorem constructs this
operation.  Whole-law exchange or simultaneous rotation of both copies is an
ensemble symmetry and does not perform the mixing.

**Route B verdict:** the direct-product Schwinger sphere is conditionally
established as classical symplectic reduction, while its particle-internal
and analyzer interpretations are blocked by the missing same-field coherent
coupling, invariant two-copy sector, and actual exchange path.

## 4. Route C: prepared response columns form one controlled oscillator block

The accepted C-CST optical construction contains two exact real
polarizations with a nonzero cross KKS entry.  Those columns are a canonical
pair and therefore one `z`.  Its three Cartesian observation targets provide
additional prepared source columns, but the accepted theorem is a
finite-window source-to-history inverse.  It permits the source to depend on
the requested history and accuracy.  It does not identify an autonomous
four-real-dimensional invariant bundle or a degenerate second oscillator.

The full `theta,G,S,current` map is nevertheless useful: once a physical
doublet exists, two independent classical current rows could provide analyzer
ports.  At present it is simultaneous classical tomography.  A prepared
linear combination of histories is not an Euler-generated noncommuting
operation, and its `o(|K|^2)` fixed-window remainder is not an all-time leakage
bound.

**Route C verdict:** an exact classical one-oscillator KKS block and
finite-window analyzer readout are established at the accepted scope.  The
two-mode autonomous Schwinger transfer is blocked by the missing second mode,
invariant projection, physical `U(2)` mixing, and complement estimate.

## 5. Quantization and measurement remain separate even after CP1

Suppose a future carrier earns (1)--(6).  Geometric quantization would require

    (1/(2 pi hbar)) integral Omega_red=2 I_0/hbar=N in Z. (8)

This is the same integral class exposed in `0043/0047`.  Euler similarity
rescales the physical KKS action continuously, so neither `I_0` nor `hbar` is
selected by (8).  The `N=1` Hilbert space has dimension two only after the
action class and quantization rule are supplied.  The classical Stokes
functions do not become measurement operators merely because their Poisson
brackets are `su(2)`.

Composing the classical analyzer with `0049` also remains conditional.  A
first-exit probability requires a measure over unresolved detector state and
physical port fluxes.  Exponential clocks reproduce normalized intensities
only after their rates are postulated.  Reversible Euler cannot implement a
many-to-one reset as an autonomous map on finite invariant phase volume; a
scattering construction must retain radiation, while external
repreparation is an added operation.

Finally, no step changes the Galilean Euler boost law or elliptic pressure.
A finite group-speed branch can type a controlled signal band, but exact
Lorentz covariance, action selection, and detector probabilities remain
independent obligations on the same carrier.

## 6. Strongest result and continuation

The Schwinger--Hopf construction is exact and does evade the narrow
finite-CCR obstruction: two supplied positive canonical modes reduce to a
compact classical `CP1` sphere with Stokes `su(2)`.  Current Euler evidence
does not supply the physical doublet.  Every inspected positive pair is
either one oscillator counted twice, a prepared finite-window response, or
two disconnected carriers without coherent mixing.

The next positive construction is therefore specific: on the persistent
carrier sought by P2, exhibit two independent positive Riesz modes; derive
their physical KKS normalization and exact or controlled frequency matrix;
compute an Euler interaction that spans two noncommuting Stokes rotations;
and bound the complement in the same energy/graph norm.  Only then should the
detector, integral action class, exchange path, and finite-speed band be
joined.  This route remains active and is not a P4 no-go.
