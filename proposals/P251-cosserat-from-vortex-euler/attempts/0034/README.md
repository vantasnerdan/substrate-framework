# Attempt 0034: one-action long-wave correspondence

Active P251 N3–N5; original objective remains unchanged pending explicit user
scope decision. Candidate: closed six-filament cells, completed point-vortex
action in 0032, plus axial phase gradients from the local-induction self-energy.
This attempt derives the algebra of the candidate's long-wave continuum;
it does not assert that its thin-core truncation is an exact smooth-Euler
theorem. No numerical remainder is needed for the correspondence algebra.

## Frozen model and approximation boundary

The two triangular populations have phases phi (inner microstructure) and
beta (cage), r^2=S/2. The collective kinetic Hessian is
M11=M22=11 I/10, M12=-9 I/10, where I is the relative-angle inertia from
0032. The angle energy is K*(phi-beta)^2/2. The covariant field
psi=(11 phi-9 beta)/2 obeys psi->psi+theta under common rotation. Completing
the square yields Jpsi=2 I/55, Jbeta=4 I/11, Kpsi=4 K/121, with zero
kinetic cross term. This field map is derived, not a choice of two unrelated
alpha coefficients.

The finite-core self-energy coefficient T is the declared Biot-Savart annulus
integral rho Gamma^2 log(b/a_core)/(4pi). For six actual helical filament
curves, line length gives

    H_gradient = 3T [sqrt(1+r^2 phi_s^2)+sqrt(1+r^2 beta_s^2)-2].

Its quadratic gradient coefficient after the field map is
Cpsi=12Tr^2/121 in Cpsi psi_s^2/2. Terms containing beta_s are retained in
the exact quadratic expression; setting beta=curl u/2 makes them higher
spatial derivatives of displacement. The microscopic approximation here is
LOCAL INDUCTION: interaction and core corrections to the gradient energy
remain to be bounded; this is not the full three-dimensional Biot-Savart
second variation. Fixed local circulation, radii, and the affine cell frame
are explicit ensemble premises; their joint dynamical compatibility remains
part of the continuum construction.

For an isotropic ensemble of scalar axial phases psi_n=n·Psi and local
cage angle beta_n=n·curl(u)/2, define L_cell as the bundle-axis length per
volume. Each bundle contains six filaments, so the original proposal density
is L_v=6 L_cell at zero axial gradient. Second/fourth moments give

    alpha=L_cell Kpsi/12, j=L_cell Jpsi/3,
    c_tr=L_cell Cpsi/30, c_s=L_cell Cpsi/15, c_a=0.

The first script execution used the symbol L_v for cell density. The explicit
normalization correction is captured in stdout-density.txt; the first stdout
is preserved, and its algebraic density variable is interpreted as L_cell.

The translational kinetic term rho|u_dot|^2/2 must come from the same
coarse-grained Euler action, including cancellation or retention of the cell
translation/rotation symplectic terms. This script does not postulate its
microscopic derivation. With it supplied, the conditional N4/N5 balances
follow by variation. The extra positive macro gradient inertia changes the
full PDE, but the exact scalar characteristic equation shows it does not
change either transverse omega^2 branch through order k^2. That is a
specified asymptotic equivalence, not identical equations at finite k.

The remaining achievements are a controlled finite-core/interaction error,
the joint translational action and affine frame map, and a smooth stationary
or explicitly permitted relative-equilibrium tube bridge. The proposed
original exact claim stays open. The full derived quadratic action and its
retained remainder are the next objects for individual review.
