# Cross-delta volume and Piola correction receipt

The frozen README equations (7)--(10) asked for a volume-preserving map from
one fixed reference torus to physical Cao cores whose major radius, core
radius, and generally volume vary with `delta`. That literal target is
geometrically inconsistent unless the physical volumes happen to agree.

The executed construction therefore uses the stronger correct separation:

1. the physical dilation/volume factor is retained as
   `Jbar_delta=2 I_a(delta)` relative to the fixed polar measure;
2. the action--angle map, or equivalently a normalized Hanzawa--Moser
   correction, makes the *relative* core Jacobian constant;
3. the raw cross-delta collar map is only orientation preserving and uses the
   full variable-Jacobian contravariant Piola transform
   `J_Phi^-1 D Phi`;
4. volume-preserving diffeomorphisms are asserted only for dynamically
   accessible displacements within one fixed physical carrier.

Thus cross-delta transport is an operator conjugacy between different carrier
spaces, whereas a same-carrier determinant-one displacement is a physical
coadjoint motion. The correction changes the domain-map interpretation but
does not change the 0048 leading coupling, spacing, or `O(1/log(1/delta))`
window ratio.
