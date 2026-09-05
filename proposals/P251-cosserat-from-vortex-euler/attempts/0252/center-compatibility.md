# A localizable nondegenerate axisymmetric center is a stagnation circle

Let u=v_r(r,z)e_r+w(r,z)e_theta+v_z(r,z)e_z be a C2 axisymmetric
steady Euler velocity in a neighborhood of the circle(r0,z0), r0>0,
with a single-valued axisymmetric pressure p (normalized by density).
Let v=(v_r,v_z), assume v(r0,z0)=0 and det Dv(r0,z0)!=0.
If u·grad p=0 throughout that neighborhood, then w(r0,z0)=0.
The result is local and needs neither analyticity nor a compact support
assumption. In particular it applies to a nondegenerate elliptic center.

## Direct Euler proof

Write A=Dv at the center, with A_ij=partial_j v_i. Axisymmetry turns
localizability into v·grad_(r,z) p=0. Differentiation at v=0 gives

    A^T grad_(r,z) p=0.

A is invertible, so p_r=p_z=0 there. The actual radial steady Euler
equation, retaining the cylindrical centripetal term, is

    v_r partial_r v_r+v_z partial_z v_r-w^2/r=-p_r.

At the same center it reads p_r=w^2/r0. Since r0>0 and w is real,
w=0. Both poloidal components already vanish, so the circle consists
of stagnation points rather than a nonzero closed streamline.

## Independent first-integral derivation

The toroidal Euler equation gives v·grad(rw)=0. Steady Bernoulli
conservation and localizability give v·grad(|u|^2)=0. The same
invertibility argument makes both gradients vanish at the center.
From grad(rw)=0 one obtains w_r=-w/r0 and w_z=0. Because v=0,

    partial_r |u|^2=2w w_r=-2w^2/r0,
    partial_z |u|^2=2w w_z=0.

Vanishing of the first gradient again gives w=0. This uses the physical
cylindrical metric factor r, not a flat product metric on a torus.

## Exposing the hypotheses

A pure toroidal field w(r)e_theta is steady with p_r=w^2/r and is
localizable. Its poloidal derivative is zero, so it does not meet the
nondegeneracy hypothesis and can have nonzero circular streamlines.
Conversely, the reviewed0211 ring has a nonzero nondegenerate poloidal
center; its steady radial pressure gradient is necessarily nonzero.
The theorem proves it cannot be pressure-localizable near that center.
It does not invalidate its actual nonlocalizable Euler construction.

Thus0248B's proposed exact localizable matching to the0211 nonzero
elliptic core is refuted at that specific shared-core requirement. Replacing
that core by a zero-speed cavity changes the supplier and does not preserve
its promised geometry. General smooth compact velocity constructions without
localizability, the common-frame periodic steady route, and fixed periodic
Beltrami compact-core backgrounds remain distinct live candidates.

The continuation is therefore to preserve the actual nonzero compact core
and solve either the global nonlocalizable matching/steady shape equation,
or the fixed-background actual response lift. No parent exhaustion follows.
