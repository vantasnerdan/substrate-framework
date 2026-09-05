"""Classify the first trigonometric comparison without changing the field."""

import sympy as s

from substrate_framework.euler_acoustic import triangular_euler_array

x, y = s.symbols("x y", real=True)
field = triangular_euler_array(1, 1, 1, (x, y))
p, v = field.pressure, field.velocity
gradient = s.Matrix([s.diff(p, x), s.diff(p, y)])
chi = 1+p
modulated = chi*v
base_euler = v.jacobian((x, y))*v+gradient
residual = modulated.jacobian((x, y))*modulated+chi**2*gradient
chain_difference = (residual-chi*(v.dot(gradient))*v-chi**2*base_euler).applyfunc(s.expand)
print("Exact expanded chain-rule difference:", chain_difference)
print("Actual base Euler residual:", base_euler.applyfunc(s.trigsimp))
assert chain_difference == s.zeros(2, 1)
assert base_euler.applyfunc(s.trigsimp) == s.zeros(2, 1)
print("The failed direct trigsimp was an unresolved representation, not a changed Euler coefficient.")
