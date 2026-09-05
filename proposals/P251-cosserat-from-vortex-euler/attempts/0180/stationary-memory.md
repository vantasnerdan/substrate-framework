# Actual stationary stress kernel, preparation slip, and physical closure

## 1. Eliminate only the actual Euler/Lin cell

Use the actual operator from `stationary-internal-action.md` on
Z=(chi,w):

    G=[[C,Id],[0,L_E]], B_Y=(0,F_Y), B_V=(0,F_V).

Let

    D0 Y=P_kappa<a²>Y,
    O_v z=P_kappa<a z+u(kappa.z)>,
    O_p z=P_kappa<(kappa.grad p)z+grad p(kappa.z)>,
    O=(O_v C+O_p, O_v).

Every row is a literal full Euler pressure/stress integral. The actual
physical mean equation through its second spatial jet is

    X_tt=k²{D0 X+O Z}+O_T(k³),
    Z_t=G Z+B_Y X+B_V X_t,                             (1)

with initial Z0=(0,-P[aD+kappa(u.D)]). Replacing the bare affine inputs
by X,X_t affects the displayed second coefficient only beyond its order,
as proved in the action file. At the exact coefficient level one can
instead put X0=D+tV and V, with no approximation in that coefficient.

The stationary group gives

    Z(t)=exp(tG)Z0+integral_0^t exp((t-s)G)
                                      [B_Y X(s)+B_V X_t(s)]ds.

Therefore the derived stress memory is

    X_tt=k²[D0 X+S_D(t)D
         +integral_0^t {K_Y(t-s)X(s)+K_V(t-s)X_t(s)}ds]+
         higher spatial jets,
    S_D(t)D=O exp(tG)(0,-P[aD+kappa(u.D)]),
    K_Y(t)=O exp(tG)B_Y, K_V(t)=O exp(tG)B_V.           (2)

The kernels are time-translation invariant; the preparation slip S_D is
not omitted. Whole-field averaging at common K averages these actual
operator matrix elements, not independently selected frequencies. This
is the inherited autonomous internal-state representation of a generally
nonautonomous two-column observed action.

The causal kernel in(2) is an initial-value elimination of(3). Merely
inserting a retarded kernel into a single symmetric quadratic action
would vary to an advanced-plus-retarded kernel. The exact action remains
the local inherited action(3), with its initial cotangent data; alternatively
one must retain the boundary/preparation terms in a two-endpoint reduced
action. No unproved passive-memory or positive spectral measure is asserted.

## 2. A physical finite-moment closure criterion

For a finite list of actual time-independent stress/momentum rows L_i,
closure on an autonomous finite state requires its row span to be invariant:

    L G=A_fin L, L B_Y=B_fin,Y, L B_V=B_fin,V.           (3)

Initial slip must also factor through the same actual initial moments.
If only a chosen preparation/input subspace is claimed, the equality
need hold on its reachable space, not on every Euler perturbation. This
is the weaker and appropriate prepared-history condition. It is not an
all-k invariant-manifold requirement.

The registered material angle and tag spin are generally moving rows
L_tag(t), because the tag and its reference quadrupole are transported.
Their correct test is

    dL_tag/dt+L_tag G=A_fin L_tag                     (4)

with the actual reference phase/amplitude, initial G and integrated-current
rows included. A constant physical optical clock could make(4) close;
choosing a rotating scalar coordinate by name does not. The independent
0176 construction addresses precisely that physical clock. Full Noether
spin and the registered sheet angle remain different observations until
their actual moment equations match.

For a positive Cosserat realization, closure additionally preserves the
actual finite symplectic mass, positive restoring coefficients and the
hybrid physical field map0172. Equation(3) is an algebraic way to test
that construction, not a claim that these physical rows already close.

## 3. An exposing actual-memory test, not an isolated-mean gate

For any one nonzero-frequency optical pair with constant coefficients,
the order-k² displacement-column acceleration generated after elimination
is a constant plus one cosine (and a sine for unsymmetrized data).
Hence its even derivatives obey

    R_D''''(0)=-omega² R_D''(0).                        (5)

The actual bare two-wave Euler preparation0175 has R_D''(0)=0 but
R_D''''(0)=-2/25. Thus its displayed stress row and preparation slip
cannot be represented by that SINGLE optical pair with constant physical
coefficients. A positive initial gap or the first two zero time jets is
not a replacement for this test. A larger state, different actual
preparation or a moving observer may change the conclusion.

This is not a no-go for the coupled parent. Genuine Cosserat theories do
produce mean memory when Phi is eliminated; (5) checks its particular
one-pair structure, not the absence of memory. The actual same-EPS
corrected preparation0179 may change the stress slip and moments. The
next construction must test its ACTUAL rows, rather than transferring
the bare two-wave derivative obstruction to a different preparation.

The strongest result here is a derived stationary Euler stress-memory
representation with full physical currents and initial phase. The finite
physical-angle/spin closure remains an active construction, not a renamed
infinite-memory solution to the requested finite Cosserat objective.
