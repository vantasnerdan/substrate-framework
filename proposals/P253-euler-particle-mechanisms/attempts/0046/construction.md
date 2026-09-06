# Exact compact-pair KKS test for a physical quantum two-state sector

## 1. Route A: the accepted block is a classical canonical oscillator

Use coordinates `x=q Q+p S` on the actual `C-CST-008` tangent pair.  Its
forms are

    Omega=B dq wedge dp,
    H_2=(Hq q^2+2N q p+P p^2)/2,                  (1)

where `B!=0`, `Hq>0`, and `Hq P-N^2>0`, equivalently the symmetric Hessian
`H=[[Hq,N],[N,P]]` is positive definite.  With the physical convention
`i_X Omega=dH_2`, the linear generator is

    A=(1/B)[[N,P],[-Hq,-N]],
    A^2=-(Hq P-N^2)I/B^2.                         (2)

Thus the pair is one positive classical harmonic oscillator with

    omega_*^2=(Hq P-N^2)/B^2>0.                   (3)

This statement is exact for the tangent forms.  It does not make their span
an invariant nonlinear Euler submanifold.  On the compact `C-CST-018` field,
the optical target in one Cartesian component has canonical coordinates
`(Phi,Pi=j Phi_t)` and

    Omega=dPhi wedge dPi,
    H_opt=Pi^2/(2j)+j nu^2 Phi^2/2.               (4)

The physical initial forms are normalized to (4) exactly, while the prepared
histories satisfy its evolution only through the accepted `o(|K|^2)` error on
each fixed window.  Neither version is a free exact two-state dynamics.

If an action unit `hbar` and the CCR are added, quantization of (4) gives the
infinite ladder `E_n=hbar nu(n+1/2)`, `n=0,1,...`.  Keeping only
`{|0>,|1>}` is not closed under `Phi` or `Pi`: the raising part sends
`|1>` to `|2>`.  More invariantly, no positive finite dimension `d` admits
an exact CCR truncation of this noncompact phase plane with
`[Q,P]=i hbar I`, because

    tr([Q,P]-i hbar I)=-i hbar d !=0.              (5)

The three optical components give more oscillator modes, not a selected
two-level system.  A fixed-one-quantum sector would itself be a new occupation
constraint and is not preserved by the physical coordinate readouts.

**Route A verdict:** the exact positive KKS oscillator structure is
established.  Its direct truncation as an exact physical two-state CCR system
is refuted by the noncompact phase plane, infinite CCR ladder,
finite-dimensional trace obstruction, and absence of an invariant Euler
oscillator manifold.  This obstruction does not refute a separately
constructed compact `su(2)` orbit or reduction; that is the positive route
continued below.

## 2. Route B: the measured rows form translations, not a compact spin algebra

The scalar rows along the registered core axis are, directly from
`C-CST-008`,

    theta=q,   G=(B^2/P)q,   L=B p.                (6)

The Poisson bracket induced by (1) gives

    {theta,L}=1,  {G,L}=B^2/P,  {theta,G}=0.       (7)

These are constant central brackets: the affine observation algebra is a
Heisenberg translation algebra on `R^2`.  It is not
`{J_i,J_j}=epsilon_ijk J_k` and has no fixed-radius `S2` coadjoint orbit.
Quadratic functions of one canonical pair instead close to
`sp(2,R)=sl(2,R)`, whose relevant orbits are noncompact; that representation
change does not supply the missing compact spin sphere.

The larger `C-CST-017/018` angle, `G`, `S`, hybrid, and current inventory is
an invertible classical source-to-observation map.  Its rows can be read
simultaneously on one Euler history.  The accepted construction contains no
measurement interventions whose operators fail to commute, no outcome
probabilities, and no map from coherent phases to detector frequencies.
Tomographic controllability of classical initial data therefore remains
distinct from quantum incompatibility and the Born rule.

**Route B verdict:** the exact physical affine Poisson algebra is established
and refutes the proposed compact-`S2` reading of this block.  The full current
is valuable classical readout, but no accepted implication turns it into a
noncommuting measurement algebra.

## 3. Route C: ensemble symmetry is not an autonomous exchange loop

`C-CST-018` supplies a fixed positive whole-state translation/O(3)/time-
reversal law over disjoint compact rings.  The law rotates or translates the
field, lattice, tags, sources, and observer together.  This is a symmetry of
the ensemble and a device for isotropic averaging; it is not an Euler
trajectory that exchanges two localized carriers.

The accepted perturbations are prepared for a selected history and compact
time window.  No theorem gives a collision-free invariant two-ring sector,
an Euler path swapping its centers, an endpoint identification in the same
phase space, or the induced KKS holonomy.  Consequently the abstract `Z2`
characters in `0043` cannot be evaluated on this compact assembly, and the
Haar law cannot select either sign.

**Route C verdict:** the direct two-copy transfer is blocked by the named
missing invariant exchange path and phase-space inclusion.  This is not a
no-go for a future persistent carrier family.

## 4. Route D1: a stochastic variational extension states the missing law

A concrete extension can introduce a configuration probability `varrho`, a
phase `S`, a mass parameter `m`, and a positive diffusion coefficient `D`,
with

    varrho_t+div(varrho grad S/m)=0,
    S_t+|grad S|^2/(2m)+V
       -2mD^2 Delta sqrt(varrho)/sqrt(varrho)=0.    (8)

Then the algebraic substitution

    Psi=sqrt(varrho) exp(i S/(2mD)),
    hbar_ext=2mD                                      (9)

is the Schrödinger representation of (8).  Equations (8)--(9) make clear what
has been added.  `varrho` is a probability density on carrier configuration
space, not the constant Euler material density; `D` selects the action scale;
the osmotic term and probability reading are new dynamical/statistical laws.
Ordinary deterministic incompressible Euler neither contains them nor fixes
`D`.  Applying Brownian diffusion directly to a material tag would change
its exact transport equation; a volume-preserving stochastic Euler version
would additionally require specified noise fields and their covariance.

This route is a coherent conditional extension, but it does not derive P4
from the allowed Euler substrate.  Its next achievement would be a separate
foundational proposal deriving the stochastic law and `D` from resolved
Euler complement statistics on the same persistent carrier, followed by an
actual detector model.

## 5. Route D2: Kelvin circulation is conserved but continuously valued

Kelvin's theorem preserves a circulation `Gamma`, but supplies no integer
condition.  Under the exact Euler similarity

    u_AB(x,t)=A u(Bx,ABt),

the corresponding circulation is

    Gamma_AB=(A/B)Gamma.                             (10)

It therefore varies continuously without changing the carrier's
dimensionless knot or link type.  Requiring a phase
`exp(i m integral u dot dx/hbar)` to be single-valued would impose
`m Gamma/hbar in 2pi Z`; the phase, `m`, and `hbar` in that sentence are the
desired quantum structure, not consequences of Kelvin conservation.  Hopf
and linking numbers quantize dimensionless topology, while helicity weights
them by continuously variable circulation products.

**Route D2 verdict:** bare topological circulation selection is refuted by
the exact continuous similarity (10).  A minimal circulation/action quantum
can be declared as a new substrate axiom, but cannot be credited as derived.

## 6. Joined result and positive continuation

The accepted compact-pair block is more than an analogy: it contains an
actual positive KKS phase plane, exact physical action normalization, and a
complete classical current readout.  Its exact algebra identifies it as a
canonical oscillator/Heisenberg system.  That structure does not select a
two-dimensional Hilbert space, noncommuting physical measurements, a Born
rule, an exchange character, an action unit, or a finite invariant speed.

The most economical same-substrate continuation is an actual compact
rotation orbit of a persistent carrier, because `S2` rather than `R2` supplies
the correct finite spin phase space.  It must be joined to a physical
threshold detector or other state-dependent interaction that realizes
noncommuting rotations and records outcomes, and to a dynamical scale-
selection mechanism.  The active `0048` carrier construction addresses the
first geometric prerequisite.  A later measurement attempt must derive its
detector dynamics and statistics from the retained Euler state `0042`, rather
than relabeling prepared target histories as quantum amplitudes.
