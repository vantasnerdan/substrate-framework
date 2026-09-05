# An actual transverse Euler acoustic route: compensated triangular cells

The parent validated the preregistration (263/12) before source inspection.
Fix a triangular lattice Lambda of primitive area A. A large finite
supercell contains N lattice points and has area NA. Its mean velocity
is fixed to zero. Write J(x,y)=(-y,x), gamma=Gamma/A, Gamma>0, rho>0.
The point idealization has physical total vorticity

    omega=Gamma sum_j delta_(z_j)-gamma.

The negative constant is part of Euler vorticity, not an imposed rotating
container. Under an area-preserving material flow it remains constant
exactly. Its velocity is included through the zero-mean torus Green
function Delta G=delta-1/(NA), not added or omitted by hand. The point
system is the singular vortex model; the smooth actual Euler realization
is established separately in smooth-dynamics.md.

## 1. Full Biot--Savart Hessian and the physical momentum pole

Use a nonzero supercell Bloch vector k with |k| below the first reciprocal
lattice vector. Let xi be the vortex-position displacement amplitude.
The primitive phase-averaged vorticity and velocity are exactly, to first
order in displacement in the point model,

    delta omega_mean=-i gamma (k dot xi),
    v_mean=-gamma Jn xi_L,  n=k/|k|, xi_L=n dot xi.                 (1)

This follows from the actual Fourier inverse curl: u_q=-i Jq omega_q/|q|².
There is no longitudinal strip cancellation: v_mean is transverse to k.
The slow velocity is actual fluid velocity, not a renamed vortex velocity.

Let f_eta(q)=exp(-eta² q)/q. The radially regularized full energy Hessian
per physical area is

    D_eta(k)=rho gamma² sum_(G in Lambda*) [
      (G+k)(G+k)^T f_eta(|G+k|²)
       -G G^T f_eta(|G|²)],                                    (2)

where the second summand at G=0 is zero. This is derived by expanding
the complete kinetic Fourier energy, including its delta²omega term;
the latter supplies the subtracted diagonal and cancels rigid core
translation self-energy. The Gaussian is a rotationally symmetric
regularization of the point self-energy, not a stationary Gaussian-core
ansatz. Its zero-radius limit defines the same periodic Green Hessian.

The G=0 term tends to rho gamma² n n^T. It is precisely the kinetic energy
rho|v_mean|²/2 from (1). Removing it would remove the actual fluid mass.
The rest has a smooth even expansion near k=0. After second derivatives,
the sole conditional sum is assigned by the displayed physical radial
regularization. Derivatives of order three and higher of the unregularized
nonzero-G symbol are absolutely summable in two dimensions; thus the
remainder after the constant and second jet is bounded by C_Lambda |k|^4.

## 2. Derived positive shear, with the ultraviolet cancellation exposed

Choose n=e_x, transverse displacement e_y. Sixfold lattice symmetry
gives the exact orbit averages

    <sin² theta>=1/2, <sin² theta cos² theta>=1/8.

Consequently the k² coefficient of the nonzero-G yy entry of (2) reduces
to

    rho gamma² sum_G [q² f_eta''(q)/4+q f_eta'(q)/2]
      =rho gamma² sum_G eta^4 q exp(-eta² q)/4.                  (3)

The equality is exact on each sixfold reciprocal orbit. Poisson summation,
or the Riemann sum after p=eta G, gives

    lim_eta->0 sum_G eta^4 |G|² exp(-eta² |G|²)/4=A/(16pi).

The reciprocal density is A/(2pi)² and the radial integral is pi/4.
Therefore the shear coefficient is derived, not supplied:

    mu=rho Gamma²/(16pi A)>0.                                   (4)

Equivalently for any smooth radial core regularizer h(q) with h(0)=1
and sufficient decay, put f=h/q; the continuum radial integral is
integral_0^infinity (q² f''+2q f')dq/4=1/4. It is the boundary term
[q² f']/4, not a finite cutoff artifact. Termwise omission of the
regularization would incorrectly give zero on every sixfold orbit.

The trace of (2) tends to rho gamma²: Poisson summation cancels the shifted
and unshifted Gaussian sums, leaving the omitted zero mode. Sixfold
symmetry fixes the second jet for every propagation direction. Thus

    D_LL=rho gamma²-mu|k|²+O(|k|^4),
    D_TT=mu|k|²+O(|k|^4),
    D_LT=O(|k|^4).                                               (5)

For all sufficiently small nonzero |k| the *full* two-by-two matrix is
positive definite. This is in the fixed zero-mean compensated Euler
frame; no frequency winding or moving-frame energy sign was selected.

## 3. Same-action mass and actual transverse mean displacement

The vortex KKS density on nonzero Bloch modes is rho gamma dxi_L wedge
dxi_T. The constant background contribution vanishes for periodic
Hamiltonian displacement generators, since its integral is a Poisson
bracket. Accordingly one phase convention for the quadratic action is

    L= -rho gamma xi_L partial_t xi_T - (1/2)xi^T D(k)xi.          (6)

The sign in (6) gives partial_t xi_T=-gamma xi_L at leading order,
the same measured velocity as (1). Eliminating xi_L yields

    L_red=(rho gamma)²/(2D_LL)
             [partial_t xi_T+D_LT xi_T/(rho gamma)]²
                       -D_TT xi_T²/2.                          (7)

The displayed real two-component convention can be applied separately
to cosine/sine sectors; inversion symmetry makes D real. The physical
mean velocity in (1), using the first Hamilton equation, is

    v_mean,T=(rho gamma²/D_LL)
                    [partial_t xi_T+D_LT xi_T/(rho gamma)].      (8)

Hence actual coarse material displacement U, defined by U_t=v_mean with
the matching initial mean displacement, satisfies U_T=xi_T+O(k² xi_T)
over a fixed number of acoustic periods. There is a genuine observer
bridge, not just equality of coordinate frequencies. In U the leading
action and physical dispersion are

    L_mean=rho |U_t|²/2-mu |grad U|²/2 + retained higher jets,
    omega²=mu|k|²/rho+O(|k|^4),
    c_T²=Gamma²/(16pi A).                                       (9)

The transverse polarization and finite ambient mass rho follow from
the same complete pressure/kinetic action. The O(k²) mass correction in
(7)--(8) is retained, not identified with an additional sound speed.

The k=0 harmonic torus velocity is a separate conserved datum. Only
nonzero Bloch displacements with sum_j Gamma xi_j=0 are used here.
They admit local Hamiltonian material lifts with no net harmonic impulse.
A common translation of every vortex need not be silently identified
with that same fixed-mean Kelvin leaf. Nothing in (9) relies on doing so.

## 4. Primary comparison and scope

Tkachenko1966, Stability of Vortex Lattices, equations29,38,42,46,
https://jetp.ras.ru/cgi-bin/dn/e_023_06_1049.pdf, treats a rotating ideal
point-vortex lattice. Its long-wave triangular speed agrees with (9)
after circulation normalization. The paper also exhibits a square-lattice
unstable mode; that whole-spectrum result is not needed for (2)--(9).
Triangular geometry is selected here because its exact sixfold moments
close the isotropic second jet. A square lattice lacks those moments;
the neutral multi-species alternative remains available for other goals.

The positive theorem above is an exact small-k point-lattice action with
an actual Euler-velocity observation. The accompanying smooth-dynamics
construction transfers finite histories to genuine finite-core Euler;
it does not claim an exact fixed-core Bloch eigenbranch or a complete
three-dimensional Cosserat/EPS continuum.
