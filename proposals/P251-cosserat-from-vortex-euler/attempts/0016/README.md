# Attempt 0016 — N2 classical corroboration: m=2 coefficient CONFIRMED against Kelvin 1880; m=1 constant sharp split discovered

## Source

NSF review "Centrifugal Waves in Tornado-Like Vortices" (Monthly Weather
Review 149, 2021), archived at `sources/NSF_review_2021_centrifugal_waves_rankine.pdf`
(md5 64aa495d324740094071e457ec93d490, text extracted). It restates Kelvin
(Thomson 1880, his Eq. 50) as its Eq. (96), read visually from the rendered
page (pypdf text mangles primes/subscripts; Eq. (96) verified against the
page image):

  (1/(bR)) J'_m(bR)/J_m(bR)  -  2 Om m/(g b^2 R^2)  =  -(1/(kR)) K'_m(kR)/K_m(kR)

with b^2 = k^2(4 Om^2 - g^2)/g^2, g = the intrinsic (Doppler) frequency
= omega - m Om, and conditions Eqs. (86)+(93): kinematic (u_r matching)
and dynamic (pressure matching). Note: plain J_m (not I_m) — the earlier
transcription guess of I_m was wrong.

## Result 1 — m=2 channel: EXACT corroboration

Solving Kelvin's (96) directly for m = 2 gives c = -1/6 = -0.16666667 at
k = 1e-2, 1e-3, 1e-4 — identical to the campaign's J-species coefficient
(attempt 0015). The campaign m=2 root also satisfies (96) with relative
residual O((ka)^4-scale). **omega = Om(m-1) - Om(ka)^2/6 is Kelvin 1880's
own m=2 branch.** The tangential condition is implied in this channel
(consistent with the m=2 closure structure).

## Result 2 — m=1 channel: sharp structural split (open question, arbiter named)

Kelvin's (96) own m=1 slow branch gives

  c1^Kelvin = 1/4 - gamma  = -0.3272  (extracted: -0.32730/-0.32700/-0.32680
  at k = 1e-2/1e-3/1e-4, 1/ln drift toward the asymptote),

while the campaign's four-condition pose (kin pair + pressure continuity
with advected Bernoulli + frozen sheet strength [v_th] = 2 Om eta) gives
c1 = 1/2 - gamma. The difference is exactly the constant 1/4 inside the
log-bracket, and the K96-residual evaluated at the campaign root is the
pure constant +0.25 across all decades — i.e. the two dispersions differ
precisely by the content of the tangential (frozen sheet strength)
condition, which Kelvin's two-condition pose does not impose.

Status of the tangential condition: exact at k = 0 (validated against the
exact displaced-vortex Euler solution, attempt 0010); the delta-vorticity
content 2 Om eta of the displaced step is exact frozen-in geometry; the
O(eta^2) status of its finite-k corrections is argued in attempts
0013-0014. The question "does the linearized Euler m=1 eigenmode of the
Rankine vortex carry c1 = 1/4 - gamma (no frozen-strength constraint
independent of the kinematics) or c1 = 1/2 - gamma (with it)" is sharp,
decidable, and NOT settled by this turn.

## Route verdict and next route

- Route verdict (corroboration route): m=2 established against Kelvin
  1880; m=1 SPLIT discovered, evidence_scope REPRESENTATION_SCOPED on the
  matching-pose question.
- Route 0017 (named arbiter): direct numerical linearized-Euler
  eigen-solve on the Rankine profile (discretize the linearized operator
  about v0 on a radial grid; the step handled by the frozen-vorticity
  structure, no interface-condition bookkeeping) and read the m=1
  eigenfrequency at small k: it must land on c1 = 1/4 - gamma or
  1/2 - gamma. This decides which matching pose the true Euler mode
  realizes and feeds the frozen verifier as an additional oracle.
