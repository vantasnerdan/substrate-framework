"""Exact finite convention blocks before the native/oracle composition."""
import json
import sys
from pathlib import Path
import sympy as s


def main(path):
    mass, k, t = s.symbols('m k t', positive=True, real=True)
    j = s.Matrix([[0,1],[-1,0]])
    h = s.diag(k,1/mass)
    a = j*h
    assert a == s.Matrix([[0,1/mass],[-k,0]])
    assert a*j+j*a.T == s.zeros(2)
    x, xp, p = s.symbols('Cxx Cxp Cpp', real=True)
    c = s.Matrix([[x,xp],[xp,p]])
    w = c+s.I*j/2
    row_t = (s.eye(2)+a*t)[0,:]
    row_0 = s.Matrix([[1,0]])
    greater = -s.I*(row_t*w*row_0.T)[0]
    lesser = -s.I*(row_t*w.T*row_0.T)[0]
    ret = s.simplify(greater-lesser)
    assert ret == -t/mass
    assert s.simplify(s.I*lesser.subs(t,0)) == x
    assert s.simplify(s.I*lesser.subs(t,0)+x) == 2*x  # wrong lesser sign is exposed
    # Eliminating a bath coordinate from z^2 M+K gives a minus Schur term;
    # D^R(z)=-(z^2 M+K)^-1 implies D=D0+D0*(V Db V^T)*D.
    aa,bb,v = s.symbols('A B v', nonzero=True)
    full = s.Matrix([[aa,v],[v,bb]])
    exact = -full.inv()[0,0]
    dyson = 1/(-aa+v*v/bb)
    assert s.simplify(exact-dyson) == 0
    # Noncommutative F and Fcc ordering cannot be collapsed.
    q = s.Matrix([[1,2],[0,1]])
    sigma = s.Matrix([[2,0],[3,1]])
    difference = -q*sigma+sigma*q
    assert difference != s.zeros(2)
    output = {'canonical_generator':str(a),'linear_retarded_position':str(ret),
              'initial_i_lesser':str(s.simplify(s.I*lesser.subs(t,0))),
              'schur_dyson_difference':str(s.simplify(exact-dyson)),
              'noncommuting_order_difference':str(difference),
              'scope':'encoded exact Cartesian/sign/block identities; not full numerical accuracy'}
    Path(path).write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output))

if __name__ == '__main__':
    main(sys.argv[1])
