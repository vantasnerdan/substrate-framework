# Mathematical steering: the actual adjoint pair and a scaled transfer estimate

Author: main-model supervisor, 2026-09-06. Continuation of the existing 0037
analytic routes. Inputs are the exact column oscillator in 0030 and the full
linearization in 0034. These calculations supply constructions to transfer;
they do not assert the unfinished solitary-wave scattering theorem.

## 1. The physical column adjoint pair is explicit

For a fixed nonzero axial frequency k, write

    eta_t = i k b chi,          chi_t = i k d K_k eta,
    b=2L/r^4,                 d=L'/r,           a=b/d.

Use the physical core inner product, with all poloidal exterior energy in K_k:

    <q1,q2>_E = <eta1,K_k eta2>_r + <chi1,a chi2>_r,
    <f,g>_r = integral conjugate(f) g r dr.

The positive critical radial eigenfunction f=f_1(k) solves

    K_k^-1 f = b d f/c_1(k)^2.

Set N=<f,K_k^-1 f>_r>0. The two exact right eigenvectors of the column
matrix H_k (G_0=i k H_k) and their dual amplitudes are

    r_+ = (K_k^-1 f, +d f/c_1),
    r_- = (K_k^-1 f, -d f/c_1),
    ell_+(eta,chi) = [<f,eta>_r + <b f,chi>_r/c_1]/(2N),
    ell_-(eta,chi) = [<f,eta>_r - <b f,chi>_r/c_1]/(2N).    (A1)

Here f,b,d are real; complex amplitudes use the Hermitian inner product.
Substitution gives H_k r_+=c_1 r_+, H_k r_-=-c_1 r_- and
ell_s(r_t)=delta_st. In particular,

    ||r_+||_E^2=||r_-||_E^2=2N,
    <r_+,r_->_E=0,
    ell_s(q)=<r_s,q>_E/(2N).                              (A2)

Thus P_s q=r_s ell_s(q) is an orthogonal rank-one projection in the
physical energy metric, with norm one. The apparently differentiated
eta-entry is bounded precisely by duality between K_k^-1/2 f and K_k^1/2 eta.
An arbitrary radial vorticity moment does not have this property.

The overall physical factor 2 pi rho multiplies N and the pairing equally
and cancels from P_s. No action or symplectic normalization is being selected.
At each nonzero k these formulas retain the full Bessel exterior. Their
measurable direct integral is energy-bounded; the set k=0 has zero measure.
Uniform differentiability/localization estimates near k=0 remain a separate
task and must keep the Bessel logarithm.

With the translating-frame convention of 0037 the eigenvalues are
i k(c+c_1) and i k(c-c_1). The latter is the slow channel. The other top
branch, lower radial branches and any stationary column kernel belong to
the fast complement. A smooth low-frequency localization is not itself an
idempotent projection: retain P_s fiberwise and use an explicit frequency
partition for estimates rather than silently calling the product a projector.

## 2. What the regular-label conservation coordinate can and cannot do

For the smooth regular-label D of supervisor-projection-analysis.md, let

    A_D(k)=<D(L),K_k^-1 f>_r,
    B_D(q)=<b D(L),chi>_r.

On span(r_+,r_-),

    M_D=A_D(ell_+ + ell_-),
    B_D=c_1 A_D(ell_+ - ell_-).                            (A3)

At small k, A_D has the positive threshold limit already derived. The
density and its column flux companion therefore form an invertible two-by-two
coordinate system on the critical pair. This is a useful physical interpretation
of (A1), not a proof that these two moments annihilate the other radial modes.
On arbitrary q, the Q contributions to both moments must be subtracted before
using (A3). Their time derivative brings back the Q equation. In particular,
conservation of M_D alone cannot establish the desired derivative factor for
the exact projected P equation without checking these extra terms.

## 3. A solvable transfer problem shows the actual velocity-ratio loss

Before asserting a cancellation of the slow residence time, solve the model
homological equation that the proposed proof must generalize. Let a,b>0,
a!=b, let B be smooth with B and B' integrable as needed below (compact
support suffices), and let L>0. Consider two co-directed scalar transports

    p_t+b p_z=0,
    q_t+a q_z=epsilon partial_z[B(z/L) p].                  (A4)

For compact test p define the outgoing transfer

    T p(z)=epsilon integral_0^infinity {
         B'((z-a t)/L) p(z-(a-b)t)/L
         +B((z-a t)/L) p'(z-(a-b)t)} dt.                 (A5)

It solves (a partial_z T-b T partial_z)p
=epsilon partial_z(B p). Put d=a-b. Integrating the second term by parts
in t gives

    T p = (epsilon/d) B(z/L)p(z)
          -(epsilon b/d) integral_0^infinity
                 B'((z-a t)/L) p(z-d t)/L dt.             (A6)

The integral in (A6), after y=z-dt, has kernel on one half-plane with
absolute value

    (1/(|d| L)) |B'((a y-b z)/(d L))|.

Its row and column integrals are at most ||B'||_1/a and ||B'||_1/b,
respectively. Schur's estimate proves the length-independent bound

    ||T||_2->2 <= |epsilon|/|a-b|
             [||B||_infinity + sqrt(b/a) ||B'||_1].        (A7)

The half-plane depends on the sign of d; the same estimate applies for either
ordering of the two positive velocities. Since (A6) is bounded, it extends
from compact test data to L2 and the homological identity holds on its core.

If the source is instead epsilon B(z/L) p_z, the same calculation gives

    ||T_input|| <= |epsilon|/|a-b|
             [||B||_infinity + sqrt(a/b) ||B'||_1].        (A8)

This is the consequence of derivative placement, rather than an appeal to
conservation alone. In the scaled slow/fast setting, a_s~mu, a_f~1 and
epsilon~mu: the slow-output derivative transfer is O(sqrt(mu)), the
fast-output derivative transfer is O(mu), and an input derivative reverses
which direction pays the square-root velocity ratio. Even the worse
O(sqrt(mu)) term is small, but it is not the claimed automatic removal of
every mu-dependent loss. These bounds are uniform in the spatial length L.

## 4. Transfer this mechanism, not its verdict

For 0037 the next calculation is to apply (A1) to the actual variable-coefficient
Euler perturbation V_mu and retain the exact kernels P V_mu Q and Q V_mu P.
Their radial derivatives can be integrated against the smooth physical dual
vectors when legitimate, but this earns explicit coefficient/domain estimates,
not an assumed bounded perturbation of the full transport generator.

Classify each resulting term as an output derivative, an input derivative,
or a genuine zeroth-order term. The last class needs its own weighted
estimate or a proved cancellation; neither (A7) nor (A8) covers it by name.
For dispersive column branches replace the scalar transport calculation by
the exact monotone frequency maps and physical group velocities. Track
derivatives of the radial eigenvectors and the k^2 log|k| exterior symbol.
The low-frequency singularity is precisely where a formal integration by
parts or a putative uniformly bounded inverse can fail.

The positive target is a small transfer/feedback estimate in the unchanged
physical norm. A stationary rank-one graph, a finite-frequency matrix or a
scalar KdV statement is not substituted for that estimate. Conversely a
failed instantaneous graph need not end the route: bounded incoming/outgoing
wave operators with the actual mode/translation rows can supply the needed
propagation control.
