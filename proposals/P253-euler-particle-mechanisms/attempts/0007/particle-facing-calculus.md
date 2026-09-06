# Particle-facing calculus on the actual carrier and background

These are exact invariant/linear-symbol results and a controlled far-field expansion. They test specific direct identifications before particle language is attached to them. The object remains a classical incompressible Euler field; no electromagnetic field map is assumed.

## 1. Compact-vorticity impulse and velocity tail

Let omega be smooth, divergence-free, compactly supported in B_a(0), and let u be the unique decaying divergence-free velocity with curl u=omega. Define the kinematic impulse I=(1/2) integral y cross omega dy and physical impulse P=rho I. The vector potential and velocity are

    A(x)=(1/4pi) integral omega(y)/|x-y| dy, u=curl A.       (1)

These are the same full-space Hodge formulas stated in Dávila et al., arxiv2207.03263, equation(1.2), whose exact source applicability is recorded in0002. All moment identities below follow independently by integration by parts.

Divergence of y_i omega shows integral omega_i=0. Divergence of y_i y_j omega gives integral(y_i omega_j+y_j omega_i)=0. Thus the first-moment matrix is antisymmetric and

    integral y_j omega_i dy = epsilon_{j i k} I_k.         (2)

Taylor expansion of 1/|x-y| now yields, for |x|>2a,

    A(x)=I cross x/(4pi |x|^3)+O(M2/|x|^3),
    u(x)=[3 n(n dot I)-I]/(4pi |x|^3)+O(M2/|x|^4),       (3)

where n=x/|x| and M2=integral |y|^2 |omega(y)|dy. The constants are uniform when a/|x|<=1/2, by the differentiated Newton-kernel Taylor remainder. This is the leading exterior velocity of a compact-vorticity carrier with nonzero impulse. It has no velocity monopole. Zero impulse removes the displayed leading term; it does not eliminate all possible fields.

For two translated compact-vorticity cores separated by d>a1+a2, exact full-space integration by parts gives

    E_cross=rho integral u_1 dot u_2 dx
           =rho integral A_1 dot omega_2 dx.             (4)

The decays in(3) justify this energy boundary term. Expanding the Newton kernel in the exact double integral, the zero zeroth moments eliminate constant/linear terms and the two first moments give

    E_cross = rho/[4pi d^3] *
        [3(I_1 dot n)(I_2 dot n)-I_1 dot I_2] + R.        (5)

For d>2(a1+a2), a direct third-order kernel bound gives the conservative estimate

    |R| <= C rho/d^4 *
        integral integral |omega_1(y)| |omega_2(z)|
                            (|y|+|z|)^3 dy dz.           (6)

Equation(5) is the cross term in POSITIVE Euler kinetic energy. Its sign is not obtained by borrowing a maintained-current magnetic dipole potential. Taking minus its separation gradient becomes an actual force only after deriving the carrier's symplectic/kinetic reduction and retaining shape/ambient terms. Ring translation commonly uses impulse as a Hamiltonian coordinate; a material-tag mass cannot be substituted without that derivation.

The leading kernel is anisotropic d^-3 energy, not a signed Coulomb d^-1 potential. This excludes only direct identification of this bare far-separated kinetic cross term with the entire electromagnetic interaction. A different measured field map, nonlinear bound structure, or a mediating energetic background remains a distinct construction.

## 2. What no-swirl rings supply as spin

For a centered axisymmetric no-swirl field,

    x=r e_r+z e_z, u=u_r(r,z)e_r+u_z(r,z)e_z,
    omega=omega_theta(r,z)e_theta.                       (7)

Pointwise u dot omega=0. Also

    x cross u=(z u_r-r u_z)e_theta.                      (8)

Its angular integral vanishes on every centered axisymmetric tag/cutoff domain. The same holds for the tag's centered internal angular momentum, since its centroid velocity is axial. If the full angular-momentum integral converges, it too is zero. Finite kinetic energy alone is insufficient to assert absolute convergence of integral x cross u; the finite-tag result and pointwise helicity statement need no such assertion.

A directed axisymmetric ring with nonzero impulse has SO(2) as its connected rigid-rotation stabilizer; if there is no additional stabilizer, its physical orientation orbit is SO(3)/SO(2)=S^2. A rotation about its own axis does not change this Euler field. A full 2pi rigid spatial rotation gives a loop on S^2, whose fundamental group is trivial. Therefore the ORIENTATION-ONLY state space of this carrier does not have the nontrivial rotation loop invoked for odd Hopf solitons in0003. This is not a claim about the topology of the full fluid orbit, its internal excitations, or its quantum completion.

The axisymmetric no-swirl source theorem is valuable for classical persistence, but neither its circulation nor its centerline orientation is already spin half. The next candidate needs a physically nontrivial internal state, swirl/twist, or a proved collective quantum sector on the same field. Swirl existence is a concrete external route, not an invented requirement: Cao–Zhan, [arxiv2009.13210v2](https://arxiv.org/abs/2009.13210v2), constructs steady axisymmetric rings with and without swirl. Only its abstract has been used here; its stability, admissible invariants and exact rotational orbit require their own source-body audit before reuse.

## 3. Propagation about a quiescent Euler background

Linearizing actual constant-density Euler about u_0=0 gives

    partial_t v=-grad pi/rho, div v=0.                   (9)

For a nonzero Fourier wave k, incompressibility and divergence of(9) force pi_hat=0 under the usual decaying or periodic pressure normalization. Thus on the transverse two-dimensional velocity space the exact linear generator is zero. A uniform background U_0 merely gives

    (partial_t+U_0 dot grad)v=0                          (10)

after pressure elimination: the frequency is U_0 dot k, a convective shift removable by a Galilean change of frame, not two propagating photon helicities with omega=+/-c|k|.

This is an exact test of the bare quiescent/uniform Euler sector, not a global prohibition of emergent propagation. A structured energetic Euler background can supply nonzero gradients, vorticity and collective restoring effects. Its actual linear/nonlinear spectrum, currents and same-carrier coupling must be derived; the prepared target-frequency response of C-CST-017/018 does not supply an autonomous photon or neutrino band.

## 4. Route verdicts and active continuation

The compact-vorticity far-field expansion and same-field kinetic cross term are established as stated under the moment/separation hypotheses. Directly equating that bare interaction to a Coulomb potential is refuted by its anisotropic d^-3 kernel.

The no-swirl helicity and finite-tag angular momentum vanish as stated. An orientation-only spinorial-loop identification for the directed S^2 carrier is refuted at that orbit scope; the full orbit and internal quantum sector remain open.

The quiescent/uniform Euler transverse generator has only rest/convective frequencies as stated. Identifying it directly with finite-speed electromagnetic or relativistic particle propagation is refuted for that background.

These results activate intrinsic swirl/twist and structured-background routes. Neither particle objective is reduced or refuted; both still require an actual shared quantum/relativistic realization and measured electric/weak currents.
