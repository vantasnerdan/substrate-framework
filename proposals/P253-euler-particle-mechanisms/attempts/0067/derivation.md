# Smooth carrier-medium recoil and scalar-charge test

## 1. One autonomous action closes reciprocity and momentum

Let `s_a` be smooth compact radial form factors with integral one, and set

    f(x;X,y)=sum_a s_a(x-X_a) F_a(y_a).                 (1)

On divergence-free `U`, take the frozen joint action

    A=integral dt [L_car(X,Xdot,y,ydot)
       +integral {rho_U |U_t|^2/2-mu |grad U|^2/2+f dot U} dx].  (2)

Here `L_car` must be supplied by the same carrier eventually; (2) is a
conditional autonomous carrier-medium extension, not a derivation from bare
Euler. Exact variation gives

    rho_U U_tt-mu Delta U=P_T f,                        (3)
    d/dt partial_(Xdot_a,k)L_car-partial_(X_a,k)L_car
       =integral s_a(x-X_a) F_a dot partial_k U dx,      (4)

and the analogous internal torque from `partial_y F_a`. Thus the field source,
carrier force, and work all come from one coupling.

If `L_car` is translation invariant, the field has sufficient decay for the
following momentum to converge, and the complete Euler--Lagrange system holds,
simultaneous translations give

    P_k=sum_a partial_(Xdot_a,k)L_car
         -rho_U integral U_t dot partial_k U dx.         (5)

Using (3), incompressibility, decay, and integration by parts shows
`dP_k/dt=0`:
the field derivative is minus the right side of (4). This is actual global
carrier-plus-medium recoil, unlike the prescribed-source comparison in 0064.
The theorem is conditional on the carrier action and source law in (2).

The static minimizer satisfies `-mu Delta U=P_T f`. With `G^T` the Oseen
kernel, its on-shell energy is

    E_on=-1/2 <f,G^T*f>,
    E_12=-<f_1,G^T*f_2>.                                (6)

Smooth form factors make every self pairing finite because `1/|x-y|` is
locally integrable in three dimensions. For separation much larger than the
supports,

    E_12=-F_(1,i) G^T_ij(d n) F_(2,j)+O(d^-2).           (7)

This supplies a reciprocal `1/d` interaction and finite recoil model exactly,
but its sign and angular tensor must still be compared with electric charge.

## 2. One vector source is not scalar charge

For `F_a=kappa q_a a_a`, `|a_a|=1`, (7) becomes

    E_12=-kappa^2 q_1 q_2
       {a_1 dot a_2+(a_1 dot n)(a_2 dot n)}/(8*pi*mu*d)
       +O(d^-2).                                        (8)

The leading interaction depends on both internal orientations and the line of
centers. Reversing `q` gives an opposite vector source, but it does not remove
the tensor dependence. With the positive elastic field energy in (2), common
parallel source signs have the scalar-mediator attraction sign; flipping the
common coupling convention cannot reverse the product sign. The accepted
transverse elastic mode therefore supplies neither spin-independent Coulomb
coupling nor its electric like-charge sign from this one-vector route.

If `a` is the Stokes vector of a Schwinger--Hopf doublet, the same conclusion
holds classically. Haar averaging makes the one-source mean zero. A quantum
singlet expectation can turn `a_(1,i)a_(2,j)` into a multiple of
`-delta_ij`, but that imports the very state, tensor-product, probability and
measurement structure P4 is required to derive. It cannot define the charge
of one isolated carrier.

**Route A verdict:** the smooth autonomous recoil/action construction is
established conditionally; its interpretation as scalar electric charge is
refuted by the exact tensor, sign, and unsupplied carrier constitutive law.

## 3. A symmetry-completed multiplet reveals the extra field content

Introduce three independent transverse fields `U^A` and a carrier body frame
`R_a`. The covariant source `F^A_(a,i)=kappa q_a (R_a)_(iA)` gives

    E_12=-kappa^2 q_1 q_2
      tr[R_1^T G^T(d n) R_2].                            (9)

When a separate mechanism locks `R_1=R_2`, cyclicity gives

    E_12=-kappa^2 q_1 q_2 tr G^T(d n)
        =-kappa^2 q_1 q_2/(2*pi*mu*d),                  (10)

which is isotropic. This is a useful exact comparator: finite internal
completion can remove the angular tensor, but only with a triplet of massless
vector sectors and an exact common-frame locking rule. Relative to the one
accepted prepared sector, this means two additional independent sectors and a
physical identification of the accepted one as the third component. For
general relative frame, (9) remains orientation dependent. The positive-energy
sign in (10) also remains scalar-attractive.

A two-complex-mode `CP1` carrier supplies one Stokes vector, not three
independent mediator flavours or a common-frame lock. Neither 0051 nor any
accepted claim supplies the field multiplet in (9).

**Route B verdict:** isotropic `1/d` response is established for the declared
three-field/common-frame extension and refuted as a consequence of the current
single transverse Euler-derived sector.

## 4. Neutral compensation preserves pair terms but removes isolated charge

If `sum_a F_a=0`, the leading far field of the complete assembly cancels:

    U(x)=sum_a G^T(x-X_a)F_a=O(r^-2).                   (11)

The individual cross terms (7) remain in the finite-separation energy, and
(5) still gives total momentum. Thus a neutral carrier-plus-background cell
can contain an exact `1/d` internal interaction while its external field is
dipolar. This is compatible with global recoil, but it does not assign a
nonzero scalar charge to the isolated composite. Separating one sign from its
compensator reintroduces the oriented source and the additional background
degree of freedom.

**Route C verdict:** established as an autonomous neutral recoil assembly,
refuted as an isolated signed scalar-charge construction.

## 5. Continuation

This attempt closes the smooth-source and Noether rows that 0064 left open.
It also shows exactly why they are insufficient: a single transverse source
remains a vector, while an isotropic completion uses a triplet, meaning two
extra mediator sectors beyond the current one, plus frame locking. The next current-substrate candidate is a nonlinear
configuration-dependent topological source whose far zero mode evades the
linear homogeneous theorem without losing finite energy or reciprocity. A
healthy propagating scalar with positive energy has the same same-sign
attraction problem as (6), so it is not by itself the electric alternative.
The minimal declared extension with the Gauss-law electric sign is a
constrained `U(1)` gauge field coupled to a conserved current. The reviewed
0042 transported compact integrable tag supplies the conditional current
template `j^0=g chi`, `j=g chi u`, but `g`, the Maxwell normalization and the
finite speed are independent foundation constants, and gauge invariance plus
the current domain must be proved on the enlarged state.

No electron, magnetic moment, action selection, Born rule, neutrino, P5, or
parent result follows.
