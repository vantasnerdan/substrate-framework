"""Dependency-free Gamma curvature maps for NumPy and exact SymPy matrices."""


def gamma_tangent_increment(j_columns, curvature, displacement):
    return j_columns @ (curvature @ displacement)


def gamma_transpose_increment(j_columns, curvature, dual):
    # Bilinear transpose, deliberately not a Hilbert adjoint.
    return curvature.T @ (j_columns.T @ dual)
