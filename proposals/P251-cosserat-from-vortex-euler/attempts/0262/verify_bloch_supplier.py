"""Exact product-rule and projector checks; analytic estimates remain in prose."""

import sympy as sp


def main():
    x, y, z = sp.symbols("x y z", real=True)
    kx, ky, kz = sp.symbols("kx ky kz", real=True)
    coords = (x, y, z)
    wave = sp.Matrix([kx, ky, kz])

    def div(vector):
        return sum(sp.diff(vector[index], coords[index]) for index in range(3))

    def curl(vector):
        return sp.Matrix([
            sp.diff(vector[2], y)-sp.diff(vector[1], z),
            sp.diff(vector[0], z)-sp.diff(vector[2], x),
            sp.diff(vector[1], x)-sp.diff(vector[0], y),
        ])

    def curl_k(vector):
        return curl(vector)+sp.I*wave.cross(vector)

    def derivative(vector, direction):
        return vector.jacobian(coords)*direction

    # Exact nonconstant steady constant-curl field: curl u=u, div u=0.
    velocity = sp.Matrix([sp.sin(z), sp.cos(z), 0])
    assert curl(velocity) == velocity
    assert div(velocity) == 0
    potential = sp.Matrix([x*y, y*z, z*x])
    displacement = curl_k(potential)
    assert sp.simplify(div(displacement)+sp.I*wave.dot(displacement)) == 0
    ordinary = derivative(displacement, velocity)-derivative(velocity, displacement)
    bloch = ordinary+sp.I*wave.dot(velocity)*displacement
    full = curl_k(displacement.cross(velocity))
    assert sp.simplify(full-bloch) == sp.zeros(3, 1)
    assert sp.simplify(full-ordinary) != sp.zeros(3, 1)

    px, py, pz = sp.symbols("px py pz", real=True)
    reciprocal = sp.Matrix([px, py, pz])
    total = reciprocal+wave
    projector = sp.eye(3)-total*total.T/total.dot(total)
    force = sp.Matrix(sp.symbols("f1 f2 f3"))
    assert sp.simplify(total.cross(projector*force)-total.cross(force)) == sp.zeros(3, 1)
    # Curl_K force is affine in K before projection; all mean rows are finite.
    raw_force = displacement.cross(velocity)
    for first in wave:
        for second in wave:
            assert raw_force.diff(first, second) == sp.zeros(3, 1)
    print("Bloch product cancellation, nonconstant Beltrami field, ordinary-bracket "
          "mutation, full Fourier projector curl, affine force constraints: PASS")


if __name__ == "__main__":
    main()
