# Actual laboratory optical pair and its measured continuum angle

## 1. The clock is obtained before the action

In the exact comparison flow u=Omega t cross x, fix an axial profile h(z)
with h'(0)=1. In the radial plateau of0123 the material displacement is
xi=a(t)h(z), a perpendicular to t. Write J a=t cross a. The Kelvin/Lin
solution has a_dot=-Omega J a, and its Eulerian velocity is
v=a_dot h-(xi.grad)u=-2Omega J a h. Direct substitution into
v_t+(u.grad)v+(v.grad)u+grad p=0 gives zero with constant perturbation
pressure in this plateau. The core-vorticity angle is Phi=J a, so the
LABORATORY angle obeys Phi_dot=-Omega J Phi. Its nonzero frequency is
|Omega|. It is not the co-rotating frequency2|Omega|: with the actual
rotation matrix R_dot=Omega J R, R^-1 Phi obeys a rate -2Omega J.
Both observations are retained, rather than selecting a Floquet winding.

The globally compact completion is0123's actual divergence-free xi,
with its Leray pressure and return field. The plateau formula is its
comparison, not an exact compact invariant subspace. Every statement
below first concerns this exact leading profile system; its finite-time
actual compact remainder is then carried explicitly in section4.

The full KKS integral is beta J0 on the two angle columns, where
J0=-J and beta=2rho Omega C integral h^2. Its generator is -Omega J;
therefore its Hamiltonian matrix is -Omega_form*generator=
beta Omega identity. This follows from the actual generator and KKS,
not from prescribing an oscillator energy. Let M=beta/Omega>0 for this
leading system. The leading measured tag spin is j_tag Phi_dot, with
j_tag=rho D integral_tag z h.0123 first constructs an ideal profile root
with j_tag=M, then a nearby compact pressure-corrected root satisfying
the stronger actual bracket identity det([angle;spin])=beta. For that
compact root the exact scalar canonical inertia is beta/W(0), not
beta/Omega: its measured Wronskian includes the pressure correction.
Its actual spin-rate inertia equals beta/W(0) exactly and is close to
M. Merely observing that both coefficients are positive would not imply
either equality. Sections2--3 use the leading matched system; section4
keeps the complete compact-root corrections.

## 2. Time reversal joins actual histories rather than coefficients

Pair complete time-reversed backgrounds with equal weights1/2. Their
reference phase forms are respectively beta J0 and -beta J0, and their
positive Hamiltonians are both beta Omega identity. Let their actual
angles be Phi_plus and Phi_minus. They obey

    Phi_plus_dot=-Omega J Phi_plus,
    Phi_minus_dot=+Omega J Phi_minus.

Define q=(Phi_plus+Phi_minus)/2 and r=(Phi_plus-Phi_minus)/2.
Then q_dot=-Omega J r, r_dot=-Omega J q, and q_ddot+Omega^2 q=0.
Every q(0),q_dot(0) in the transverse plane is realized by independent
initial data r(0)=J q_dot(0)/Omega. No common arbitrary microscopic
history is imposed. Setting r=0 for all times instead forces q=0:
that shortcut would delete the actual reaction degree of freedom.

Pulling back the average of the two original actions gives, up to its
explicit total derivative,

    L_pair=beta r.dot(J q_dot)
                         -beta Omega (|q|^2+|r|^2)/2.

The independent variation of r gives r=J q_dot/Omega. Its exact
elimination yields L_q=M|q_dot|^2/2-M Omega^2|q|^2/2. The average
MEASURED spin is j_tag(Phi_plus_dot+Phi_minus_dot)/2=j_tag q_dot.
Thus the ideal moment match makes it the canonical momentum,
with no halving or doubling of a supplied rigid-body inertia. This is
an actual initial-data join of the leading Euler profiles, not the
unproved assertion that every reduced Cauchy--Born history is free Euler.

## 3. Isotropic reconstruction of a physical angle

Rotate the whole pair, including its actual tag and background, by Haar
SO(3). Let n be its axis and P_n=identity-n n^T. Prepare q_n=P_n Phi
and q_n_dot=P_n V at the initial time. Because all the comparison pairs
have the SAME measured laboratory Omega^2, their freely evolved histories
stay q_n(t)=P_n Phi(t), Phi_ddot+Omega^2 Phi=0. Averaging an arbitrary
distribution of different frequencies would not have this property;
0121 gives its exact fourth-derivative counterexample.

The measured macro angle is the least-squares reconstruction from the
actual transverse observations:

    A=E P_n=(2/3)identity,
    Phi=A^-1 E q_n.

This formula is fixed by the geometric observation experiment: each
cell measures the component of the same physical angle perpendicular
to its axis. The raw mean E q_n is NOT that angle. Equivalently Phi
minimizes E|P_n Phi-q_n|^2, whose normal equation is A Phi=E q_n.
For population intensity nu (counting a whole pair with its equal mass
weights once), the averaged action and actual spin are

    L_macro=(j/2)|Phi_dot|^2-(j Omega^2/2)|Phi|^2,
    S_macro=j Phi_dot,    j=(2/3)nu M>0.

The identical second moment A enters both expressions. If the raw mean
were renamed Phi without transforming its observation law, a spurious
factor3/2 between action and spin would result. This is an exposing
normalization mutation, not a new acceptance criterion.

## 4. Actual compact remainders and centroid response

For finite compact profiles, let the measured angle histories obey
Phi_pm_dot=mp Omega J Phi_pm+e_pm and let their actual spin rows be
S_pm=M Phi_pm_dot+d_pm, after the exact reference-time matching.0123
supplies these actual errors on a fixed finite number of periods; it
does not set them to zero. The joined equations are exactly

    q_dot=-Omega J r+e_s, r_dot=-Omega J q+e_a,
    q_ddot+Omega^2 q=e_s_dot-Omega J e_a,
    (S_plus+S_minus)/2=M q_dot+(d_plus+d_minus)/2,

where e_s=(e_plus+e_minus)/2 and e_a=(e_plus-e_minus)/2. Control of the
second-order residual consumes one derivative of e_s, or instead uses
the displayed first-order system and its Duhamel bound. A first-order
error alone is not silently differentiated. The leading positive pair
action has the actual small moving-pullback corrections from0123/0115;
they are not removed by eliminating the leading r. For fixed optical
time and bounded derivative norms these corrections remain controlled
under the same finite hierarchy. Transferring them to an actual global
EPS background is the separate construction0124, not a consequence of
local background approximation alone.

The full hybrid momentum J_H retains actual tag centroids AND all
ambient fluid.0117's exact material first jet is

    J_E=J_H+(i/2)k cross S-(i/2)dotI k+O(|k|^2).

For the completely co-rotated axial preparation, SO(3) covariance makes
the k=0 coherent symmetric shape-rate and symmetric Euler-stress rows
vanish: an isotropic rank3 tensor mapping an axial vector to a symmetric
rank2 tensor is zero. This does not delete their fluctuations or their
nonzero-k corrections. When the actual first-Bloch-jet transfer of0116/0117
applies, it therefore gives

    Delta J_H=-(i/2)k cross Delta S+O(|k|^2 C_T).

With the leading spin law this is nonzero for a transverse angle with
nonconstant rate. Its integrated form is

    Delta U_H=-(i j/(2rho))k cross Delta Phi
                       + explicit initial-momentum term + O(|k|^2 C_T).

Here rho is the TOTAL retained fluid mass density, not nu times tag mass.
An independent homogeneous Galilean mode fixes the initial-momentum
term; selecting that mode to zero extracts the optical transfer. The
same formula holds with the integrated measured spin-error contribution
when the finite compact remainder is retained. This is the physical
first-spatial-jet observation consequence, not positive shear or a
second-gradient Cosserat dispersion proof. Those remain0125's active
construction. In particular a nonzero initial-data transfer does not
establish an autonomous optical Bloch pole.

## Verdict

The lab-clock, time-reversal initial-data/action join, common-frequency
isotropic reconstruction and explicit error propagation are established
at their stated comparison and transfer scopes. The actual compact
normalization and EPS/spatial hypotheses remain attached to their named
constructions. No old alpha=L_v*T/6 formula is recovered, no accepted
claim is changed, and the stronger parent objective remains active.
