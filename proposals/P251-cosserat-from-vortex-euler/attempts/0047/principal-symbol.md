# Principal carrier symbol and the failed minimal-return curvature route

This calculation freezes a local vorticity jet `omega=W e_z`, `W>0`,
inside the exact EPS action formula while retaining its Beltrami parameter
`lambda`. A nonzero constant vorticity is not itself a nonzero-lambda
Beltrami solution. The result is a principal-symbol diagnostic, not an
existence claim about such a constant background. Use the fixed length
units of 0045; take `lambda>0`, `K>0` for notation. Reversing both helicity
sign conventions gives the corresponding other sign.

For microscopic carrier sign `sigma=+/-1`, let

```
a_sigma=(1,-i sigma,0)/2,
p_sigma=(kappa_x,kappa_y,sigma K+kappa_z),
xi_1,sigma=-i p_sigma cross a_sigma/K,
xi_2,sigma=i sigma xi_1,sigma.
```

These are the exact curl-generated slow-envelope returns. The force is
`F_i=xi_i cross W e_z`, and the velocity is projected using the full
shifted wave `p_sigma`. The norm/helicity expression can be reduced before
expansion:

```
<P_p F_i,P_p F_j>
 =F_i^*.F_j-(p.F_i)^*(p.F_j)/|p|^2,
<P_p F_i,curl_p P_p F_j>=F_i^*.curl_p F_j.
```

This exact representation avoids a needlessly huge rational expansion; no
numerical approximation is involved.

Put `r_sigma^2=kappa_perp^2+(sigma K+kappa_z)^2` and

```
h_sigma=W^2 (sigma K+kappa_z)^2/K^2
 [(kappa_perp^2+2(sigma K+kappa_z)^2)/(4r_sigma^2)
    +sigma(sigma K+kappa_z)/(2lambda)].
```

The complete matrices are

```
H=(h_++h_-) I+i(h_+-h_-) J,
Omega=W[-2i kappa_z/K I+(1+kappa_z^2/K^2) J],
J=[[0,1],[-1,0]].
```

In particular the diagonal KKS first derivative is not zero. Define
`t=kappa_z/K` and the exact Darboux map
`T=(I-i t J)/(1-t^2)`. Then `T^* Omega T=W J`, while

```
T^* H T=A I+i D J,
A=W^2[1+K/lambda
       -kappa_perp^2(1/r_+^2+1/r_-^2)/4],
D=W^2[kappa_z/lambda
       -kappa_perp^2(1/r_+^2-1/r_-^2)/4].
```

The positive frequency branch is exactly `(A+D)/W`; its partner is
`-(A-D)/W`. Through the leading spatial-gradient order,

```
A=A0-W^2 kappa_perp^2/(2K^2)+O(|kappa|^3),
A0=W^2(1+K/lambda),
D=W^2 kappa_z/lambda+O(|kappa|^3).
```

Eliminating the conjugate shape in the normalized action therefore gives

```
L=(W^2/(2A)) (partial_t q+(W/lambda) partial_z q)^2-A q^2/2
```

at this derivative order, with the transverse operator dependence of A
retained. Normalizing its inertia to `I0=lambda/(K+lambda)` yields

```
C_perp=-W^2/K^2,
C_z,lab=-W^2/[lambda(K+lambda)].
```

Both are negative. The first-derivative drift `W/lambda` remains in the
single-orientation action. If counter-drifting/parity sectors are averaged
at the action level, their mixed time/space terms cancel but the displayed
negative axial coefficient remains. Averaging their frequency squares
would incorrectly produce the opposite sign for that contribution.

Route verdict: the minimal constant-jet carrier envelope does not supply
the requested positive principal gradient action. This does **not** refute
the full finite compact EPS construction: its physical core attachment and
background-gradient corrections can matter at the small curvature scale,
and a value-only O(1) Hessian estimate cannot determine their derivatives.
The actual-field continuation is therefore performed directly in
gradient-cage-proof.md, with every required derivative bounded and a
different explicitly declared geometric cage return. No sign is transferred
from this constant-jet approximation to the full object without such bounds.
