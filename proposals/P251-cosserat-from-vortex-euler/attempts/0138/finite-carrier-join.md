# Finite-carrier physical second jet and controlled optical window

This is a construction attachment, not a change to accepted C-CST-008/009/010.
The actual input is the pressure-matched material surface branch of 0135;
its smooth/global realization is not supplied by this calculation. The
carrier calculation below is an **axial-column** calculation. Section 6 is
an explicitly conditional orientation assembly, not an assertion that
transverse modulation of a column solves three-dimensional Euler.

## 1. Actual Euler input and order of limits

Inside a Rankine core of radius a, the base velocity is Ωr eθ; outside it
is Ωa²/r eθ. The pressure-matched azimuthal mode m has pressure
P(q,r) exp(i mθ+i qz−iωt). With σ=ω−mΩ, its interior velocity is

    vr = i(σP′−2ΩmP/r)/(4Ω²−σ²),
    vθ = (2ΩP′−σmP/r)/(4Ω²−σ²),   vz=qP/σ.

The material displacement is η=iv/σ, not a relabeling of the laboratory
frequency. Put x=|q|a, s=σ/Ω and λ=x sqrt(4−s²)/|s|. The actual interface
condition is

    (4−s²) x Km′(x)/Km(x) + s² λ Jm′(λ)/Jm(λ) − 2ms = 0.

The 0135 branch has s=−1−x²/[2(m²−1)]+o(x²). For m=2 its actual material
frequency squared is

    f(q)=σ(q)²=Ω²+B0 q²+o(q²),    B0=Ω²a²/3>0.

The action and the spin observation are the actual 0135 rows. In particular,
with ηr=A_r and ηθ=iB_r,

    β=4πρΩLz ∫ A_r B_r r dr>0,   M_single=−β/(σ c²)>0.

Here c is the physical material-shape angle row, normalized by the actual
rotation generator m. The passive nonnegative material tag is
μχ(r)χz(z)[1+ε b(r)cos(mθ)]. Define

    Q=∫χ b r^(m+1)dr,   T(q)=∫χ b rP(q,r)dr,
    G(q)=∫χ r^m(P′+mP/r)dr,
    Z0=∫χz dz,          Zc(q)=∫χz cos(qz)dz.

Its angle row and measured mechanical spin/rate are

    c = Zc G/[ε Z0 Q σ(2Ω+σ)],
    j_tag = ρμmπ ε² Z0 Q T(2Ω+σ)/(σG).

Both displacement contributions to mechanical angular momentum are
included in this spin row; pressure torque is retained by the core input.
At q=0 this tag gives the wrong sign when Q≠0. Thus q0>0 is fixed first,
with Q≠0 and the tag chosen there. The independent macro wave number d is
then expanded with |d|≪q0. No positive zero-carrier rotor or nonvanishing
inertia in the q0→0 limit is asserted.

For fixed q0 let

    f0=f(q0)=ν0²,    A=f′(q0),    B=f″(q0)/2.

The actual two carrier branches obey
f±=f(±q0+d)=f0±Ad+Bd²+O(d³). Their fields and observation rows, not merely
their frequencies, are Taylor expanded about these fixed carriers.

## 2. A fixed physical tag and coherent normalization

The single-carrier equality j_tag=M_single would require

    T/Q=R(q),
    R=−β Z0 σ²(2Ω+σ)/(ρμmπ Zc²G).

`radial-jet-rank.md` proves the four actual radial rows Q,T,T′,T″ are
independent at sufficiently small fixed q0≠0. Four radial bumps therefore
prescribe Q≠0 and all three derivatives of T/Q to any smooth target.
Scaling their common amplitude makes |εb|<1 without changing T/Q.
The construction is finite at each fixed q0; its conditioning is not
uniform as q0→0. The proof uses the actual Bessel radial equation.

There is a necessary normalization repair for a coherent standing field.
Its physical branch amplitudes are one half of the angle-normalized
traveling fields. The quadratic action has weights 1/4, whereas the
linear physical angle and spin have weights 1/2. Thus the correct fixed
tag target for this standing construction is

    (T/Q)^(j)(q0) = (R/2)^(j)(q0),   j=0,1,2.

It yields j_tag=M_single/2 through the carrier second jet. Reusing the
single-carrier target gives a measured spin twice the coherent canonical
momentum. Time-reversal population averaging does not remove that factor:
it weights each action and each physical spin by the same population
weight. This repair changes the tag preparation, not the observed angle
or a supplied inertia. The radial rank lemma already allows this target.

The diagonal carrier action used here is the Fourier/Bloch action with
the two distinct axial carriers orthogonal under the full-fluid spatial
average (a compatible periodic cell or the corresponding spectral
decomposition). A finite tag is an observation, not that action average.
The analogous orthogonality in a newly constructed bounded or periodic
array has to come from its actual cell action.

## 3. The retained carrier coordinate is a material moment

Write the two angle-normalized carrier coordinates as

    θ+=q+ip,   θ−=q−ip.

The actual even axial tag observes q. The coordinate p can also be
reconstructed physically, without changing that angle. Choose χz even,
put Z(k)=∫χz(ζ)exp(ikζ)dζ and h(k)=Z′(k)/Z(k), and observe the corresponding
central material-shape moment with one extra axial factor ζ. Its normalized
row is

    r = −i h_even q + h_odd p,
    h_even=[h(q0+d)+h(−q0+d)]/2,
    h_odd =[h(q0+d)−h(−q0+d)]/2.

Consequently p=(r+i h_even q)/h_odd. At d=0, h_even=0 and
h_odd=h(q0)≠0 for an ordinary finite window and sufficiently small q0>0;
h(q0)=−q0 Var(ζ)+o(q0). The extra signed moment weight ζ is not a negative
mass density. All coordinates remain actual moments of the same
nonnegative material tag. Fourier reality exchanges the two carriers.

The pressure-mode, displacement and angle-row derivatives are included
by using the physical angle-normalized branch fields at q0±d. Their axial
derivatives need not coincide with derivatives of an arbitrarily held
canonical amplitude. An independent time-reversed history population,
when needed for the two-phase rotor preparation, does not change these
spatial identities.

## 4. Exact second jet, complete action and physical spin

Define effective masses M(q)=M_single(q)/2 and write
M(q0±d)=M0±M1d+M2d²+O(d³), where M2=M″(q0)/2. In the state (q,p),

    Mmat = [[M0+M2d², iM1d], [−iM1d, M0+M2d²]],
    Fmat = [[f0+Bd², iAd], [−iAd, f0+Bd²]],
    Hmat = Mmat Fmat.

These follow by the same change of variables T=[[1,i],[1,−i]] in the
actual carrier action: Mmat=T†diag(M_single,+,M_single,−)T/4. Thus
L=(xdot†Mmat xdot−x†Hmat x)/2, with the usual real Fourier pairing.
The mass is positive for sufficiently small |d| at fixed positive M0.
In particular the diagonal stiffness second coefficient is
H2=M0B+M1A+M2f0 and the off-diagonal first coefficient is
H1=M0A+M1f0. Freezing the mass before this multiplication loses real
observation/preparation terms.

The half-target tag gives the actual spin

    S=(M0+M2d²) qdot+iM1d pdot+O(d³),

which equals the q canonical momentum of this full action. It is not
merely M0 qdot at second spatial order. Its derivative is the physical
tag torque row inherited from the pressure-matched branch. The remaining
shape/current moments are not set to zero by this equality.

The actual equations are

    qtt+(f0+Bd²)q+iAd p=O(d³),
    ptt+(f0+Bd²)p−iAd q=O(d³).

For standing preparation p(0)=pdot(0)=0, qtt(0)=−(f0+Bd²)q(0), but
q^(4)(0)=[(f0+Bd²)²+A²d²]q(0). This exposes why averaging squared
frequencies is not a closed physical equation. With
G(t)=sin(ν0t)/ν0, elimination of the retained physical odd moment gives

    qtt+(f0+Bd²)q − A²d² ∫_0^t G(t−s)q(s)ds = O(d³).

Nonzero initial p adds its explicitly propagated homogeneous source.
Equivalently the exact two-branch model has the fourth-order operator
[(∂t²+f0+Bd²)²−A²d²]. Neither time reversal nor carrier reversal removes
its second-jet memory. An autonomous second-order identity requires the
actual group derivative A=0 (or retention of the internal state).

## 5. Positive optical equation on a controlled finite time window

The absence of an autonomous identity does not prevent a quantitative
positive optical approximation. Compare the actual standing second jet
to q_PDE solving qtt+(f0+Bd²)q=0 with the same initial q and qdot. Let

    R0=sqrt(|q(0)|²+|qdot(0)/ν0|²),
    ε_group=A²/(f0 B),             0≤t≤T=C/ν0,

with B>0. The free optical history has magnitude at most R0. At second
spatial order their difference is

    δq = A²d² (G*G*q_free).

Using |G|≤1/ν0 and |G′|≤1 gives

    max(|δq|, |δqdot|/ν0)
      ≤ A²d² R0 T²/(2ν0²)
      = (C²/2) ε_group (B d²/f0) R0.

More directly, the omitted physical memory torque is bounded by

    sup |j A²d² (G*q_free)| ≤ C ε_group [j B d² R0],

where j is the positive standing inertia (M0 per column, or the correctly
reconstructed density below). The bracket is the natural positive
gradient-torque scale. A small inertia does not hide a large relative
error, because it cancels in this ratio.

For the m=2 small-carrier surface branch, assuming its differentiated
small-q expansion as supplied by the dispersion equation,

    ε_group = (4/3)(q0 a)² + o((q0a)²).

The quadratic branch f=Ω²+Ω²a²q²/3 gives exactly
ε_group=4a²q0²/(3+a²q0²). Thus choose a small *nonzero* q0 first to meet
the desired relative tolerance on a fixed optical window, then choose
the half-target tag and only afterwards the macro band. Each resulting
construction has B>0 and j>0. The calculation does not take a singular
tag limit, claim a positive limiting inertia, or extend the bound to
beat times or propagation times diverging as d→0.

For an actual smooth branch, use its actual f0,A,B and the same exact
ratio, not substituted Rankine values. Carrier Taylor remainders must
be controlled at the fixed q0; symmetry can improve odd remainders, but
no such improvement is needed here. Core smoothing, localization and
intercore errors are separate. A fixed absolute realization error is
not automatically an O(d²) coefficient error: the next construction
needs differentiated spatial-row control, or an explicit finite macro
band with that error budget. This section licenses the actual axial
prepared second-jet approximation, not an unconstructed array.

## 6. Conditional isotropic assembly and retained physical tensor

Suppose an actual cell calculation additionally licenses columnwise
modulation d=n·K with no missing transverse Euler response. For axial
material angles, prepare q_n=n·Φ+O(K²), p_n=O(K), and use Haar directions.
The physical observation map is Φ=3 E[n q_n], since E[nnᵀ]=I/3. For
number density Ncell, j_hom=Ncell M0/3. This is not the transverse-tilt
observation map of 0126 and not the raw vector mean.

Define the retained, physically reconstructible symmetric tensor
Ψ=3 E[nnᵀp_n]. Through the prepared second jet,

    Φtt+f0Φ+B D_K Φ+iA K_j Ψ_ij=O(K³),
    Ψtt+f0Ψ=(iA/5)[I(K·Φ)+K⊗Φ+Φ⊗K]+O(K³),
    D_K=(|K|²I+2KKᵀ)/5.

Eliminating Ψ gives the same group memory with d² replaced by D_K.
The natural bound of section 5 holds separately on its positive
transverse/longitudinal eigenvalues |K|²/5 and 3|K|²/5.

The angular reconstruction is
p_n=(5/2)n·Ψn−tr(Ψ)/2, and its inherited metric is
E[p_n²]=[5||Ψ||²−tr(Ψ)²]/6>0 for Ψ≠0. This is a physical odd shape moment,
not a renamed spin or a supplied extra oscillator. The grading is
important: O(K²) higher harmonics of q and O(K³) harmonics of p do not
feed the displayed Φ jet. Arbitrary O(1) hidden initial harmonics would
require a larger state.

With H1,H2 as in section 4, the complete action at that grading is

    L/Ncell = 1/2 { M0/3 (|Φdot|²−f0|Φ|²)
      + M0/6 [5||Ψdot||²−tr(Ψdot)²
                    −f0(5||Ψ||²−tr(Ψ)²)]
      + M2/15 [|∇Φdot|²+2(div Φdot)²]
      − H2/15 [|∇Φ|²+2(div Φ)²] }
      − M1/3 Ψdot:∇Φdot + H1/3 Ψ:∇Φ.

Terms quadratic in ∇Ψ are fourth order under the prepared grading.
The actual observed spin, including its current correction, is

    S_i=Ncell [M0 Φdot_i/3
       + M2 (|K|²δij+2KiKj)Φdot_j/15
       + iM1 K_j Ψdot_ij/3] + O(K³).

It agrees with the functional canonical momentum of this action; the
spatial derivative in the kinetic action contributes to that momentum.
The leading positive-gradient optical benchmark therefore has the
derived 3:1 longitudinal/transverse angular-curvature ratio. It is not
yet an actual three-dimensional isotropic Euler constitutive theorem.

## 7. Physical mean/translation interface and route verdicts

The exact material hybrid-current identity of 0117/0125 remains in force:

    J_E,i−J_H,i = −iK_j C_ij − K_jK_l T_ijl/2 + O(K³),
    C_ij=I_dot,ij/2−ε_ijm S_m/2.

Hence J_H^(1)=ρ m^(1)+i C^(0)K and
J_H^(2)=ρ m^(2)+i C^(1)K+T^(0):KK/2, with the actual mean Euler rows
m^(j). The spin of section 6 is a genuine contribution to this identity;
it does not determine the quadrupole, shape rate, pressure-mediated
mean stress, or compensating ambient flow. A pure column harmonic may
have vanishing whole-fluid mean while its tag current is nonzero. No
independent translational spring is inserted to close that difference.

Route A, autonomous physical second-order closure by coherent carrier/TR
pairing at a generic nonzero carrier, is **refuted** by the explicit
A²d² memory term and initial fourth derivative. Route B, the actual
axial retained physical two-moment state, its positive coherent action,
fixed-tag spin matching through second order, and controlled positive
finite-optical-window approximation, is **established in the stated
column branch scope**. Its isotropic formulas are conditional assembly
identities. The stronger three-dimensional constitutive join is still
missing its actual transverse periodic Euler/Bloch and mean-current
construction, not refuted by the axial result.

A failure-derived autonomous alternative is a genuine nonzero-carrier
stationary point f′(q0)=0 with f″(q0)>0 and a nonsingular positive tag.
Real axial periodic modulation/Bragg coupling or a closed-cell optical
band can generate such a candidate. Its actual Euler pressure-cell
operator and observation/action rows must be constructed; pairing
frequencies does not establish it. The present uniform column has
A≈2B0q0≠0 near zero and cannot use its singular endpoint instead.

The next active construction is the actual periodic transverse and
intercore response with differentiated physical rows and mean stress.
The parent's Bernoulli-lift/periodic-array route and 0141 target that
missing object. No parent completion or scientific exhaustion follows
from this child.

## Exact oracle

`radial_rank_verify.py` and its frozen first output certify the four-row
Bessel rank construction (8 checks). `carrier_join_verify.py` and
`carrier-first-run.txt` certify the exposing carrier variance, physical
coherent mass/spin normalization, complete mass-stiffness product,
material odd-moment inversion, sphere-derived coefficients and metric,
and natural-scale finite-window bound (11 checks). These are exact
algebra/calculus checks, not numerical small-eigenvalue designs or
independent claim-acceptance reviews.
