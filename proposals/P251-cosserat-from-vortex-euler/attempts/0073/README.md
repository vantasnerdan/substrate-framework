# 0073 — material-centroid versus Eulerian mean momentum

This is a kinematic bridge for the active material-parcel route 0075,
motivated by the exact point-mean cancellation in 0072. It does not assume
that the source angular impulse equals a chosen parcel's complete spin.
The parent objective and all frozen scientific criteria remain unchanged.

For a material parcel D_a with constant mass m_a, centre X_a and relative
coordinate r=x-X_a, define V_a=Xdot_a and w=v-V_a. Incompressibility and
material transport give integral_D rho r=integral_D rho w=0. Define the
second mass moment I_ij=integral rho r_i r_j and the actual intrinsic spin
S_i=integral rho epsilon_ijk r_j w_k. Then

    Q_ij=integral_D rho w_i r_j
        =Idot_ij/2-epsilon_ijm S_m/2.

The symmetric identity is the material derivative of I; the antisymmetric
identity follows from the definition of S. No constitutive law, vortex
localization, or rigidity assumption enters. The exact kinetic split is
integral_D rho |v|²/2=m_a|V_a|²/2+integral_D rho |w|²/2.

For distributions tested against a smooth macro field, expand
delta(x-X_a-r) about X_a. The exact coefficient of its first derivative is

    p_Eulerian=p_centres-div(Idot_density)/2+curl(S_density)/2
                + terms of second and higher spatial order,
    p_centres=sum_a m_a V_a delta(x-X_a).

Translation contributions to the first moment vanish because r is centered.
Every parcel, including ambient fluid, must be included for the density to
equal the entire fluid density rho. The remainder has the ordinary integral
Taylor form and is controlled by the corresponding finite parcel moments;
it cannot be discarded in a claim about its second spatial coefficient.
In an isotropic reflection-paired ensemble a uniform axial rate cannot
produce a symmetric second-rank shape-rate tensor. That symmetry kills its
zeroth-order spin-to-Idot response, not all shape responses or higher jets.

Thus a physical centre momentum need not equal a pointwise Eulerian mean.
When a complete construction independently establishes S=j*angle_dot at
leading order and supplies the shape response, its spin part gives the
leading observation map U_E=U_centres+j curl(angle)/(2rho). A coordinate
normalization with this form is meaningful only after these hypotheses are
proved for the actual full material partition. Equal algebraic coefficients
alone do not establish it: the full physical spin can include background
affine motion and ambient parcel contributions that are absent from an
isolated core impulse. That is the explicit work assigned to 0075.

The verifier evaluates the exact first moment of a six-point material
quadrature with arbitrary translation, rigid rotation and symmetric strain.
It is a regression of the universal tensor identity, not a six-particle
Euler construction or a numerical approximation to the desired continuum.
The proof above applies directly to the continuous material integrals.
