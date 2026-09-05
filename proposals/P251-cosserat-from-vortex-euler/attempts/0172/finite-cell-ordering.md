# Ordered finite-cell choices on the same EPS field

Fix the actual C-CST-011 packet, tag and moment controls, carrier p, small
delta, time window, and a strict positive optical-curvature margin first.
Keep all packet-own errors relative to its phase mass and curvature.
The small curvature scales as Omega² delta/p², not as the optical gap.
The common-vector optical transfer0174 is an additional required input;
its absolute finite-time closeness alone would not suffice here.

## 1. The packet inertia stays bounded while its density is reduced

For a fixed nondegenerate local packet write
j_*= -beta_*/(gamma_* c_*²)>0. The actual local field/pressure transfer
and finite moment inverse are continuous at that packet. Choose their
neighborhood so each of beta, |gamma| and |c| lies between one half and
twice its reference value, with unchanged sign. Then

    j_*/16 <= j_packet,N <= 16 j_*.

This bound includes the actual retuned tag, not just the frozen KKS.
The large-period quadrature in0153 and0170 retains these local margins
while its normalized H² defect tends to zero. Enlarging the period does
not force the normalized local optical observation to vanish: its tag
quadrupole and moment inverse were fixed before this transfer.

For one packet per selected finite cell, whole-Haar reconstruction gives
j_N=j_packet,N/(3 Vol_N). Hence j_N remains positive at every finite N
and can satisfy j_N omega²/(4rho B_T)<eta<1. The acoustic density is rho,
not the diluted tagged mass density. The strict acoustic margin a0 of0170
is preserved by its period-independent response coefficient estimate.
No nonzero j is claimed in the infinite-volume limit.

## 2. Actual acoustic tag-shape rows need no false pointwise estimate

Reference tag moments and their time derivatives stay bounded per packet
on the fixed window by the uniform local flow bounds. Their densities,
including the actual I_tag, are O(Vol_N^-1).

For an acoustic response cell field known only in normalized L², a compact
tag integral against a fixed smooth weight w satisfies

    |Vol_N^-1 integral w chi_N|
       <= ||w||_L² ||chi_N||_L²_avg /sqrt(Vol_N).

The same estimate applies to the finite set of velocity/displacement/time
rows supplied by0170. Smooth compact material weights retain bounded L²
norm on this window. Thus their per-volume acoustic shape/current filters
tend to zero even without a period-independent pointwise Sobolev embedding.
This weaker bound is sufficient against the bulk scale rho a0. It is NOT
sufficient against diluted optical j; optical errors instead use the
packet-own estimates before dividing by volume.

## 3. Choice order and the remaining suppliers

Choose fixed microscopic geometry and accuracy on its own tiny B scale;
then the EPS arc radius/local transfer and finite moment margins; then a
large finite quadrature cell giving the bulk coefficient error, local
transfer and density inequalities simultaneously. These are open strict
conditions, not a zero-density theorem. Compare the mixed current error
with (j/rho)omega²(epsilon_R+epsilon_spin)/B_T; compare the pure acoustic
shape/current rows with rho a0, and the pure optical error with j B_T.

Finally choose macro k last, below the actual finite-cell Taylor radius
and the finite remainder bound divided by the selected j B_T. This never
uses a period-uniform radius or an unrestricted operator C² estimate.
The scalar acoustic second-jet coefficient transfer0170 is uniform in
period for its prepared response; the common-vector optical supplier0174
has its separate localization and infrared proof. Neither replaces the
other. Actual clock and coherent branch-variance remainders remain as in
`autonomy-remainder.md`; finite-cell dilution does not make them exactly
autonomous or extend the window to acoustic time 1/k.
