# P253/0087 — physical high-harmonic radiation barrier for the Cao analyzer

Owner: `root`

Status: **centrally registered and schema-activated; one bounded
postactivation README precision repair was replayed before any claim body,
verifier, source download, or numerical execution.**

## Frozen objective and authority

The full Issue #203 electron and neutrino campaign remains active.  This
attempt advances the same-carrier P2/P4/P5 radiation rung exposed by
P253/0085: determine whether the high toroidal harmonic used by the proposed
Cao Schwinger--Hopf doublet suppresses its Maxwell radiation strongly enough
over a complete analyzer gate.

Reviewed P253/0052 supplies the physical core clock

    Omega_a=F(a)/a^2=Gamma/(2 pi a^2),
    omega_physical=Omega_a sigma.                           (1)

For the Cao ring, `zeta=omega_theta/r`, `dnu=r dr dz`, and therefore

    kappa=integral zeta dnu=integral omega_theta dr dz=Gamma. (2)

Reviewed P253/0066/0073 supplies the fixed-column high-radial-index law

    sigma_J(k)=k L_Phi/(pi J)+O(k/J^2)                     (3)

at fixed positive massive `k`.  Reviewed P253/0074/0078 supplies the
fixed-`J`, fixed-massive-`k` graph/Riesz transfer and two positive-Krein Cao
modes.  Final reviewed P253/0079/0082 supplies rational-ray column crossings,
the fixed-`ell`, high-`J` alternative, and the two-sided gate ledger, but does
not supply the joint high-`J`/thin-ring transfer or physical mixer
normalization.  Corrected P253/0085 has passed independent P253/0086 review:
`0085/derivation.md` SHA-256
`3bca39e05d66891d5392390f371a84a7c748893f2fe2c5dec915523a9c67a01e`,
`0085/result.yaml` SHA-256
`7fb6c0af9a7530506471c1f69a3fe271c09bd4792ec95c32cc8ad85eb264d299`,
`0086/review.md` SHA-256
`4c5c535e60a134ba589ac066aa906249952d47d843c96e7f6a3764686f85ea67`,
and `0086/verdicts.yaml` SHA-256
`afd7e563d782d56f8d9f383a8611caf41cf4f2f4bbeb1410dc089262348f9c3d`.
Only its exact Maxwell radiation shell, necessary transverse dark-current
condition, and reviewed full constrained-generator essential inclusion are
consumed here; it supplies no resonance theorem.  Active P253/0083 is not an
accepted input; its finite-window notation may be reconciled only after its
own review.

This is an analytic candidate-comparison attempt.  It makes no production
numerical claim and does not infer a stable particle, a Born rule, charge or
action selection, a strict propagation cone for Euler pressure, an electron,
or a neutrino.

## Exact physical shell and toroidal Fourier factors

Fix a translating, nonrotating Cao ring with mean radius `R`, core scale `a`,
geometric ratio `delta=a/R`, toroidal harmonic `n>=2`, and column wave number

    k=delta n,       a=R k/n.                              (4)

Let `omega` be the physical comoving mode frequency.  For translation speed
`c_g` and Maxwell speed `c_EM`, define the strict subluminal quantities

    Delta_c=c_EM-|c_g|>0,
    C_c=sqrt(c_EM^2-c_g^2)>0.                             (5)

The real Maxwell radiation shell from P253/0085 is

    |K|^2=[omega+c_g K_z]^2/c_EM^2.                       (6)

Writing `t=K_z/|K|`, the shell gives

    K_perp
      =|omega|sqrt(1-t^2)/[c_EM-sgn(omega)c_g t].         (7)

The exact maximum occurs at `t=sgn(omega)c_g/c_EM` and is

    K_perp,max=|omega|/C_c.                               (7a)

The conservative bound `K_perp<=|omega|/Delta_c` may also be retained, but
may not be called sharp.  Equations (7)--(7a) and the positive radial-root
denominator must be rederived with the adopted time/Fourier sign.  The current is

    delta J=g(chi_g delta u+u_g delta chi),               (8)

and its transverse shell trace is the physical quantity to estimate.

For a ring-supported scalar `n` character, the exact angular Fourier
integral contains `J_n(K_perp r)`.  Converting cylindrical vector components
to a fixed Cartesian frame shifts the characters and produces
`J_(n-1)`, `J_n`, and `J_(n+1)`, with the exact `2 pi i^m` phases.  The body
must derive these factors for all components of (8), including derivatives,
Leray/Hodge reconstruction, the tag term, and real `+/-n` partners.  The
transverse projector is applied covariantly: fixed-Cartesian intermediate
terms may show additional phase shifts, but rotational equivariance must be
used to derive, rather than assume, that the final vector character reduces
to orders `n-1,n,n+1`.  On the
support, for one fixed geometric constant `C_core`,

    r<=R+C_core a=R(1+C_core k/n).                        (9)

Thus the exposing barrier condition, with a fixed margin
`eta_rad in (0,1/2)`,
is

    K_perp,max(R+C_core a)<=(1-eta_rad)(n-1).            (10)

Under (10), uniform Debye bounds give constants depending only on the fixed
margin such that every vector-frame Bessel factor obeys

    |J_(n+s)(K_perp r)|<=C_eta exp(-c_eta n),
    s=-1,0,+1,                                            (11)

up to the explicitly controlled polynomial factors from derivatives and
normalization.  Merely assuming `K_perp r<n` without a fixed margin does not
earn (11), because the turning region is then unresolved.

## Route A — fixed-J massive rational-ray modes

For a fixed massive column mode with dimensionless frequency `sigma_J(k)`,
(1)--(4) give the exact physical ledger

    Omega_a=kappa n^2/(2 pi R^2 k^2),
    |omega|=kappa |sigma_J(k)| n^2/(2 pi R^2 k^2).       (12)

Therefore (10) is exactly

    [kappa |sigma_J(k)|/(2 pi R k^2 C_c)]
       [n(n+C_core k)/(n-1)] <=1-eta_rad.                (13)

For the rational-ray construction `n_i=N p_i`, `k_i` in fixed massive
compacts, define `N_rad` as the largest integer for which (13) holds for both
members, using their actual common physical frequency at the exact crossing
and their actual `C_c(N)`.  Route A must combine this ceiling with the
subluminal ceiling and all lower graph/IVT/response thresholds:

    N_min<N<min(N_EM,N_rad).                              (14)

Since the left side of (13) grows linearly in `n` while the Cao speed also
reduces `Delta_c`, this route cannot yield an `N->infinity` radiation barrier
at fixed `c_EM`.  Its positive target is a nonempty finite interval (14) with
fully exposed constants.  If the interval is empty, that is a route-scoped
refutation of high-order Bessel suppression for the fixed-`J` construction,
not a refutation of the doublet or the full LAP/Feshbach route.

## Route B — fixed-ell, high-J modes

Take one fixed nonzero harmonic step `ell` and the adjacent radial pair

    (J,n_1=J ell,k_1=k_*),
    (J+1,n_2=(J+1)ell,k_2=k_*(1+1/J)).                   (15)

at their eventual common physical frequency.  Both wave numbers remain in
one fixed positive massive compact.  For either member, with `H=J+O(1)` and
`n=H ell`, equations (1)--(4) and the reviewed column law (3) give

    sigma_J(k)=delta ell L_Phi/pi+O(delta^2),

    omega_physical
      =kappa ell L_Phi/(2 pi^2 R^2 delta)+O(1),

    n=k/delta.                                           (16)

Here `L_Phi` is the dimensionless coefficient in the normalized column
problem; its dimensional normalization must be checked against (1), rather
than absorbed into a new length.  For an actual transferred pair at a finite
subluminal integer `J`, put

    q_0(J)=C_0/sqrt(c_EM^2-c_g(J)^2),

    C_0=kappa ell L_Phi/(2 pi^2 R k_*).                  (17)

The needed joint high-`J`/thin-ring supplier must prove, uniformly on a
buffered speed range,

    Q_i(J)=[K_perp,max(R+C a)]/(|n_i|-1)
          =q_0(J)+R_i(J),
    |R_i(J)|<=C_asym/J,                 i=1,2.           (17a)

The support, `n-1`, common-frequency, and transferred-profile corrections may
all be absorbed in `R_i`.  Since `c_g(J)` grows and the fixed-`c_EM`
subluminal range ends at finite `J_EM`, (17) is a finite-`J` leading predictor,
not a literal `J`-to-infinity limit.  Fix `eta_rad in (0,1/2)`.  Debye applies
to the actual pair if one finds an integer satisfying

    q_0(J)<=1-2 eta_rad,
    J>=J_asym,
    J_asym>=C_asym/eta_rad.                              (18)

Then (17a) gives `Q_i(J)<=1-eta_rad` for both members.  The constant
`C_asym` and hence `J_asym` are not supplied by the present fixed-`k` input;
they are outputs of the open joint transfer.  If an actual `Q_i` reaches one,
the exposed scale is Airy `n^(-1/3)` rather than exponential; above one it is
oscillatory with `n^(-1/2)` envelope away from further caustics.

A buffered leading criterion requires first
`C_0<(1-2 eta_rad)c_EM` and then

    c_g(J)^2
      <=c_EM^2-[C_0/(1-2 eta_rad)]^2.                    (18b)

This gives `q_0(J)<=1-2 eta_rad`; it yields the exact memberwise margin only
when the same integer also satisfies the quantified remainder condition (18).

Combining (18b) with the Cao speed law
`c_g(J)=kappa[log J+C]/(4 pi R)+o(1)` defines a finite radiation ceiling
`J_rad`.  The charged positive target is the nonempty integer window

    J_lower=max(J_graph,J_Riesz,J_KKS,J_gate,J_asym,...)
      <J<min(J_EM,J_rad),                                (18c)

with the support and `n-1` margins included.  It may contain one or finitely
many exponentially protected gate scales when the foundation speed is large
enough and the open transfer supplies a small enough `J_asym`; it is not an
unconditional fixed-constant asymptotic family.

At fixed `c_EM`, the Cao speed grows logarithmically along the thin-ring
sequence, so a uniform `Delta_c` cannot persist to `J=infinity`.  Route B
therefore supplies an exact leading predictor and a conditional finite-window
criterion, rather than an established protected charged pair.  It earns a
same-carrier suppression result only after constructing the already named
joint high-`J`/thin-ring Riesz transfer, with error small enough for the
crossing and the physical mode normalization.  Fixed-column asymptotics alone
do not supply that transfer.

## Gate-time comparison

Pointwise Bessel decay neither restores an exact real-frequency `L2` Maxwell
eigenmode nor proves a detector or persistence theorem.  A merely nonzero but
exponentially small shell trace can support a resonance/quasimode or a
finite-gate leakage estimate; an exact eigenmode still requires exact
dark-current cancellation.  On the physically normalized doublet, derive
bounds of the form

    ||P_T delta J_hat||_(shell trace)
       <=C H^M exp(-c H),

    T_gate<=C H^G,                                       (19)

where `H=N` in Route A and `H=J` in Route B.

including both outbound and return kernels, the diagonal action drift, the
two mixer axes, and the complete coarea weight.  A weighted outgoing
limiting-absorption/Fermi-golden-rule calculation or a time-domain flux
identity must convert the shell trace into radiated power `P_rad` or a width.
The actual comparison is `T_gate P_rad` divided by the physically normalized
mode action.  Equations (11) and (19) make that loss exponentially small only
after this conversion and only when every polynomial lower threshold,
including `J_asym` and the constant-dependent loss threshold `H_loss`, lies
inside (18c).  The
exponents `M,G`, constants, and applicable finite interval are outputs; no
polynomial gate-time assumption may be inserted without the physical KKS
and response normalization.  An exponential or uncontrolled gate time leaves
the suppression theorem open even when (18) holds.

## Route C — full shell trace and resonance continuation

If (13) or (18) fails, or if physical normalization defeats (19), continue on
the same mode by evaluating the exact transverse current trace on (6).  A
nonzero trace excludes an `L2` Maxwell eigenfield for the prescribed real
frequency.  Resonance existence, sign, and width require a weighted outgoing
Maxwell limiting-absorption estimate followed by the coupled fluid/tag
analytic-Fredholm or Feshbach construction.  Exact dark/BIC currents remain
possible because an analytic compact-source Fourier transform may vanish on
a codimension-one shell.

Route C must compare the computed radiation rate with the complete 0079 gate
ledger.  It remains active after either Bessel route fails and cannot be
replaced by genericity language.

## Verification and verdict boundary

The strongest oracle is analytic: derive (2), the physical clock (12), the
three vector-frame Bessel orders, the shell maximum, the exact finite ceiling,
and the Debye exponent with a fixed margin.  A small symbolic helper may
check substitutions and integer ceilings, but it cannot establish shell
trace regularity, a mode normalization, the high-`J` transfer, or a gate-time
bound.  No production numerics are preregistered.

Each route receives exactly one scoped verdict and activates its next rung.
The full Issue #203 campaign remains open unless every frozen electron and
neutrino conjunct is established or a separately reviewed scientific-
exhaustion certificate covers the complete candidate universe.
