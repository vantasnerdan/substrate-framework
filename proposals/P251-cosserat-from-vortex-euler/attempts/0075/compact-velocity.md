# Constant-vorticity compact-velocity candidate: exact restriction

This is a route-scoped result. A constant nonzero vorticity is a local
principal jet, not itself a nonzero-lambda Beltrami EPS field.

Let omega=W e_z, W!=0, and xi be smooth compactly supported and
divergence-free. Suppose the induced velocity v=P(xi cross omega) is
also compactly supported. Set F=xi cross omega=v+grad p. The pressure
potential has zero gradient in the connected exterior, so its additive
constant can be chosen with p compactly supported as well.

For any compact divergence-free vector field w, integration of
div(x_i x_j w) proves that M_ij(w)=integral x_j w_i is antisymmetric.
Therefore, writing N=M(xi),

    M(F)=W [[N_yx,N_yy,N_yz],[-N_xx,-N_xy,-N_xz],[0,0,0]],
    M(F)=M(v)-I integral p.

Its zz entry first gives integral p=0. Requiring the first matrix to be
antisymmetric then gives N_xy=N_xz=N_yz=0. Hence M(F)=M(v)=0, and

    integral x cross v = 0.

Thus a compact induced velocity cannot carry the desired nonzero
angular reaction in this constant-vorticity route. The conclusion is
about this local model, not general smooth Beltrami backgrounds.

There are nonetheless nontrivial exact compact velocities with physical
rotation jets. For any compact smooth g,

    xi=curl curl(g e_z)=(g_xz,g_yz,-Delta_perp g),
    v=W(g_yz,-g_xz,0).

Both are divergence-free and compact, and v=xi cross omega exactly.
Choosing g=q(x z^2/2+x^3/6) near a core gives
xi=(q z,0,-q x), a genuine core rotation about e_y. The complete compact
return carries the compensating angular moment. The KKS pairing within
this entire generator family vanishes by the integral of a transverse
Jacobian. Retaining the jet while ignoring its return would therefore
give a false angular rotor.

A materially different compact-velocity repair would use the spatial
variation of the actual EPS vorticity. Locally the exact inversion is

    omega dot grad p = -omega dot v,
    xi = omega cross (v+grad p)/|omega|^2 + a omega,
    omega dot grad a
        = -div[omega cross (v+grad p)/|omega|^2].

Compactness imposes two along-flow compatibility conditions in addition
to div v=0. No theorem that these are automatically solvable, with
nonzero spin and the required core jet, is imported here. This is a
specific next analytic construction, not a hidden condition declared
true by choosing a convenient gauge.
