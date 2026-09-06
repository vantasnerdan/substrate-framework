import sympy as s

from substrate_framework.euler_neutral_cell import (
    octahedral_rank2_orbit_sum,
    octahedral_rank3_orbit_sum,
    octahedral_vector_orbit_sum,
    proper_octahedral_rotations,
)


rotations = proper_octahedral_rotations()
assert len(rotations) == 24
assert len(set(rotations)) == 24
assert all(R.det() == 1 and R.T*R == s.eye(3) for R in rotations)
print("PASS physical proper-octahedral group")

a, b, c, d, e, f = s.symbols("a b c d e f", real=True)
vector = s.Matrix([a, b, c])
matrix = s.Matrix([[a, d, e], [f, b, d], [c, e, f]])
assert octahedral_vector_orbit_sum(vector) == s.zeros(3, 1)
assert octahedral_rank2_orbit_sum(matrix) == 8*s.trace(matrix)*s.eye(3)
print("PASS vector and anisotropic rank-two cancellation")

rank3 = s.MutableDenseNDimArray.zeros(3, 3, 3)
counter = 1
for i in range(3):
    for j in range(i, 3):
        for k in range(3):
            rank3[i,j,k] = counter
            rank3[j,i,k] = counter
            counter += 1
assert octahedral_rank3_orbit_sum(rank3) == s.ImmutableDenseNDimArray.zeros(3,3,3)
print("PASS symmetric rank-three cancellation")

r = s.symbols("r", positive=True)
assert s.integrate(r**2*r**-5, (r,1,s.oo)) == s.Rational(1,2)
assert s.integrate(r**2*r*r**-5, (r,1,s.oo)) == 1
print("PASS r^-5 lattice zeroth and first moments")
