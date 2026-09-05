# Exact compact induced velocity with nonzero physical spin in an EPS tube

This theorem closes the compact-velocity existence slice, not yet the
positive angle/reaction-pair slice. The finite-jet oracle is arithmetic,
not a PDE discretization or a soft eigenvalue computation.

## 1. Explicit right-normal matrix

Use the differential operator D and order-six S frozen in
`operator-receipt.md`. For alpha a multi-index and each input component j,
the divergence row of D S contributes s_jalpha at output alpha+e_j.
Writing C_ij=epsilon_ijk omega_k, the second row is

    sum_i partial_i composed with multiplication by C_ij.

The exact multiplication identity

    m_C partial^alpha
      =sum_gamma<=alpha (-1)^|gamma| binom(alpha,gamma)
          partial^(alpha-gamma) composed with m_(partial^gamma C)

therefore supplies the column entries

    output beta=alpha-gamma+e_i:
    (-1)^|gamma| binom(alpha,gamma) partial^gamma C_ij.

This constructs the238 by252 matrix A(x) explicitly from the order-six
jet of omega. The unknown coefficients are on the right, so no
derivative of an unknown occurs in A(x)s(x)=0. Its entries are real
analytic for the actual analytic Beltrami field.

For K_e=e cross(x-x0), the mechanical spin functional has the adjoint row

    R_(e),jalpha=(-1)^|alpha| partial^alpha(omega cross K_e)_j.

At x0 the cross-product identity gives the explicit entries

    (-1)^|alpha| [e_j sum_m alpha_m partial^(alpha-e_m) omega_m
                 -alpha_j partial^(alpha-e_j)(omega dot e)].

Terms with zero alpha component vanish. All required derivatives of
omega are of order at most six, and the angular row needs at most five.

## 2. The universal rank upper bound is235, not merely238

At any point with omega0!=0, orient auxiliary orthonormal coordinates
(X,Y,Z) with omega0=W e_Z. The formal adjoint constant-principal system is

    grad f+omega0 cross grad g=0.

For homogeneous degree-seven polynomials it has the three independent
solutions

    g=Re(X+iY)^7, f=-W Im(X+iY)^7;
    g=Im(X+iY)^7, f= W Re(X+iY)^7;
    g=Z^7, f=0.

Use only these top-degree jets as test coefficients for the output
rows of A. In every derivative of the adjoint equation up to degree six,
a derivative falling on omega multiplies a lower derivative of f or g,
which vanishes at the centre. Only the constant omega0 principal system
remains. Thus these are three exact independent left null vectors of
A(x), for EVERY smooth omega at that point. Consequently

    rank A(x)<=238-3=235.

This is a structural rank bound, not an inference from a small singular
value. The separate `left_kernel_probe.py` found precisely these
highest-output-degree dependencies in the first arithmetic experiment;
the polynomial calculation above is the analytic explanation.

## 3. Exact exposing oracle

`jet_probe.py` constructs the frozen eight-wave rational unit-eigenvalue
Beltrami field and evaluates this exact matrix modulo101. It obtains

    rank A mod101 =235;
    rank [A; R_ei] mod101 =236 for each i=1,2,3;
    rank [A; R_e1; R_e2; R_e3] mod101 =238.

None of the rational wave denominators is divisible by101. A nonzero
minor modulo101 is a nonzero rational minor. Combined with the proved
universal upper bound235, these results establish over characteristic
zero that A has rank235 and that its three angular rows are independent
modulo the row space of A. This is a proof-bearing exact arithmetic
certificate; no tolerance is used.

As an independent arithmetic replay, `exact_probe.py` uses the SAME
prototype under x->15x and omega->225omega. Every wavevector and
amplitude then becomes an integer, and these invertible rescalings leave
all row ranks unchanged. FLINT's exact integer-matrix rank gives235,
236 for each separate angular augmentation, and238 for all three. The
first outputs of both computations are retained unchanged. The initial
modular probe correctly withheld its spin conclusion until the missing
analytic upper bound was supplied.

## 4. Construct the actual compact fields, not only a pointwise jet

Choose a nonzero235 by235 minor of A and independent angular augmented
minors at the prototype point. On a small neighborhood these minors
remain nonzero. Select the235 corresponding constraint rows and pivot
columns. For any chosen free column k, put s_k=1 and all other free
coefficients zero; solve the pivot coefficients by

    s_pivot(x)=-A_pivot(x)^-1 A_rows,k(x).

These are explicit analytic functions, equivalently ratios of the
specified finite determinants. Since rank A<=235 everywhere, the other
constraint rows also vanish on this s(x). Selecting the free columns
from the angular augmented minors gives three analytic null vectors
s^(l)(x) with an invertible3 by3 angular matrix

    J_il(x)=R_ei(x) s^(l)(x).

Multiply the null vectors on the RIGHT by J(x)^-1. This retains their
right-normal form and A s=0 while normalizing the angular adjoint rows
to J_il=delta_il pointwise on a still smaller ball.

Now choose any real C-infinity compact bump psi supported in that ball
with integral psi=1. Form

    xi_l=sum_j,alpha e_j partial^alpha(s^(l)_jalpha psi),
    v_l=xi_l cross omega.

Every field is smooth and compactly supported. The differential
identities D S=0 prove div xi_l=div v_l=0 exactly. Therefore

    P_R3(xi_l cross omega)=v_l

without an exterior pressure tail. Integration by parts using the
normalized angular adjoint gives

    integral (x-x0) cross v_l = e_l.

Multiplication by rho gives the actual mechanical spin. This is a
constructive formula using only finite jets, a nonzero finite minor,
analytic matrix inversion and a compact bump. It does not assume a
right inverse for an infinite boundary trace map.

The coordinate origin is immaterial because integral v_l=0 for every
compact divergence-free velocity. Its symmetric first velocity moment
also vanishes, by integrating div(x_i x_j v_l). Thus, when the support
is inside the actual D, the field automatically has zero tube centroid
moment, zero tube STF moment, zero corresponding global moments, and
identical nonzero tube/global mechanical spin. Its normal trace on the
physical material boundary is exactly zero. This replaces, for these
directions, the noncompact exterior response mechanism of0080.

## 5. Transfer to an actual knotted EPS field

The eight-wave prototype alone is not asserted to contain an EPS tube.
Instead start from the existing analytic EPS Beltrami field with its
robust invariant trefoil domain D, at the SAME curl eigenvalue. Choose
x0 in its interior. Add epsilon times the fixed eight-wave Beltrami
prototype translated to x0 and rescaled to that eigenvalue.

Each selected235 or238 minor is a polynomial in the finite jet of omega.
In a chosen minor, each first-divergence row is independent of omega;
every second-divergence and angular row is linear in omega. Thus along
the above addition the highest power of epsilon has the nonzero
prototype minor as its coefficient. The polynomial is not identically
zero. Only finitely many epsilon values are excluded. One can therefore
choose epsilon arbitrarily small, nonzero and outside that finite set,
while staying in the existing EPS/KAM persistence neighborhood.
The perturbed actual stationary Beltrami field retains its knotted
material tube and has the required open finite-jet rank property.
Shrink the support ball to lie inside this actual invariant domain.
The determinant formulas of section4 now construct its compact
nonzero-spin isovortical directions.

The strict finite-jet/minor and tube persistence margins can also be
included in the bounded good-patch event of the Gaussian stationary
assembly. That invokes its established local support theorem, not the
prototype's non-Gaussian global law. Uniform derivative bounds of the
finite determinant formulas are obtained by recording the corresponding
positive minor margin and sufficiently many background derivative norms.

## 6. Exact scope and next pair-level computation

Established: analytic actual EPS fields with compact same-isovortical
induced velocities, zero physical tube boundary-normal response and
three independently prescribed nonzero mechanical spin components.
The constant-vorticity no-spin restriction in0075 is therefore not a
no-go for varying Beltrami fields.

Not yet inferred: a nonzero KKS angle/reaction pair within THIS compact
velocity sector, the required physical core-angle observation map, or
a positive full energy Hessian on that pair. In particular, imposing
div(xi cross omega)=0 restricts its high-frequency polarization; the
earlier unrestricted circular negative-helicity cage dominance cannot
be copied without recomputing the exact operator sector.

Next exact computation: retain the local operator columns S and evaluate
their skew bilinear differential operator S^* [omega cross] S and the
full same-Euler Hessian

    rho integral [v_xi dot v_eta-v_xi dot curl v_eta/lambda].

If their required pair exists, attach a relabeling-invariant core angle
using its actual nonzero vorticity-jet variation and carry these compact
directions into the now boundary-compatible material construction. If
the restricted pair mechanism fails, retain this compact-spin theorem
and generate the next candidate from that concrete failure. No parent
completion or exhaustion follows solely from this theorem.
