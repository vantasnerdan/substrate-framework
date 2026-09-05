# Both transverse polarizations: exact isotropic initial Euler response

This continues the registered full-phase candidate C. It retains total
mass rho and uses a whole-field SO(3) law. The result is an actual
initial-response identity, not a replacement of the cell memory by a
constant elastic modulus at all times.

## 1. Two actual physical preparations

Let u be any smooth mean-zero periodic stationary incompressible Euler
field, with nonzero energy density. Rotate its whole field and pressure
by one Haar-distributed R, and use a uniform cell-origin translation
when expressing this as a stationary whole-field law. In each realization prescribe a transverse
macro direction κ and U,V with κ·U=κ·V=0. A displacement preparation
is the actual solenoidal Kelvin generator ξ=U exp(iεκ·x); its initial
Eulerian velocity is Pε(ξ×ω). A full-phase velocity preparation is
the actual initial velocity V exp(iεκ·x). Their superposition is an
actual linear Euler initial datum; finite-amplitude Euler data are
obtained by the material map plus that small solenoidal velocity.
Their initial circulation classes are explicitly allowed to differ.

The common-V kinetic term is rho |V|²/2 in every realization. It is not
an SH projection with a subsequently rescaled mass. The physical
observables here are the full Euler Fourier mean m and its time integral,
with the initial mean displacement assigned by the material generator.

## 2. Actual pressure expansion

For a nonzero microscopic Fourier wave g, write its velocity coefficient
b=u_g, so g·b=0. The pressure projector is

    P_(g+εκ)=I−(g+εκ)(g+εκ)ᵀ/|g+εκ|².

The Kelvin initial velocity is obtained by applying this projector to
i[g(U·b)−b(U·g)]. At ε=0 it is the translated stationary field
−i(U·g)b. Its first spatial coefficient is i Z_g, where

    Z_g=−κ(U·b)
       + g[(κ·g)(U·b)+(κ·b)(U·g)]/|g|².             (1)

The exact mean Euler equation is

    mdot=−iε Pκ〈(κ·u)w+u(κ·w)〉.                    (2)

The zero-order translated-field contribution in (2) vanishes by a
periodic total derivative. Inserting (1) gives its complete second
spatial coefficient, including the actual nonlocal pressure.

For the independent common-V preparation, differentiate the actual
Euler equation once. The first spatial coefficient of w_t(0) is i Y_g,

    Y_g=−(κ·b)V+2g(κ·b)(V·g)/|g|².                 (3)

Equations (1) and (3) are different; neither is a renamed canonical
momentum. Nevertheless their whole-field isotropic physical response
coefficients agree, as follows.

## 3. The derived whole-field coefficient

Set |κ|=1 and take a transverse observation along x with κ along z.
For each g let n=g/|g| and first consider a real transverse polarization
b of unit norm. A complex Fourier coefficient is handled by its real
and imaginary transverse polarizations; the ±g reality pair removes
spurious imaginary contributions. Under a Haar rotation of the whole
field, n and b form a uniformly rotated orthogonal pair. The exact
fourth moments are

    E[n_x² b_z²]=2/15,
    E[n_x n_z b_x b_z]=−1/30,
    E[b_x²]=E[b_z²]=1/3.

They follow by conditioning b in n's orthogonal plane,
E[b_i b_j|n]=(δij−n_i n_j)/2, and integrating n on the sphere.
For (1), the pressure-return term is

    E[(n_x b_z+n_z b_x)²]=1/5.

For (3), it is

    2E[n_x² b_z²+n_x n_z b_x b_z]=1/5.

The bare term is −1/3 in both cases. Summing the actual Fourier energy
coefficients therefore gives

    c_initial² = (2/15)〈|u|²〉 > 0,

    E_R[m_displacement,t(0)]
       =−ε² |κ|² c_initial² U+O(ε³),

    E_R[m_commonV,tt(0)]
       =−ε² |κ|² c_initial² V+O(ε³).                (4)

The identity holds on both transverse polarizations by rotational
covariance. The actual full density multiplying these accelerations
is rho. It applies in particular to the same Bernoulli-lifted array
after its mean Galilean velocity is removed. It does not assert that
the unrotated axial or complementary responses agree: 0144's axial
common-V second derivative vanishes, whereas its Kelvin displacement
force is positive. The whole-field rotation, not an independent-cell
substitution, is what yields (4).

## 4. Exact closed second jet retains the actual cell kernel

For complete histories, retain the operator (3) of
`full-cell-response.md`. Define its causal zero-data Green operator G(t)
and its propagated initial cell state χ_hom. All are actual periodic
Euler/Lin response objects, not supplied oscillator inputs. Define

    F0 U=P0[(κ·∇p)U+κ(U·∇p)],
    F1 V=−2P0[(u·κ)V],
    Jκ χ=Pκ〈(u·κ)χ+u(κ·χ)〉.

Then the physically observed mean displacement X satisfies the finite-T
second-jet Volterra law

    Xtt=ε² {Pκ〈(u·κ)²〉X
       +Jκ(∂t+A0)[χ_hom+G*(F0 X+F1 Xdot)]}+O(ε³).  (5)

Here A0 acts on the cell field before Jκ. Replacing the material mean U
by X on this right-hand side is licensed by their explicit O(ε²)
current difference on each fixed interval. No cell initial state is
silently discarded. For the common-V preparation χ_hom=0; a Kelvin
displacement carries its computed cell initial data.

Every rotated realization has the same leading free macro history
X0+tV and zero first physical mean spatial jet. Therefore replacing
its X_R by the ensemble mean in the O(ε²) right-hand side changes only
higher spatial order, at fixed T. Averaging (5) is a genuinely closed
isotropic transverse mean *memory* law at this order. It keeps mass
rho. An instantaneous Cosserat constitutive coefficient requires the
additional actual kernel reduction or controlled-window approximation;
(4) alone does not supply that reduction or an acoustic-time remainder.

The 0146 continuation constructs bounded physical normal-form/current
controls for that next step. The same-array optical join separately
uses the actual material tag and hybrid-current map, as in 0148. The
present coefficient does not replace those constructions.
