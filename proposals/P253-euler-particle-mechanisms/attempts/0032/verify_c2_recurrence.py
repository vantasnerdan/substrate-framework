#!/usr/bin/env python3
import sympy as s
t,u=s.symbols('t u', real=True)
g=s.sqrt(2)/2
C0=s.Matrix([[0,0],[-g,0]])
C1=lambda x:s.Matrix([[s.Rational(3,2)*s.cos(x),s.sqrt(2)/2*s.sin(x)],[s.sin(x)/(2*s.sqrt(2)),-s.cos(x)]])
U=lambda x:s.Matrix([[1,0],[-g*x,1]])
D=s.integrate(U(2*s.pi-t)*C1(t)*U(t),(t,0,2*s.pi))
assert s.simplify(s.trace(D))==0
DD=s.integrate(s.integrate(s.trace(U(2*s.pi-t)*C1(t)*U(t-u)*C1(u)*U(u)),(u,0,t)),(t,0,2*s.pi))
print('double_C1_trace =',s.simplify(DD))
print('recurrence: F_c(g_c(sigma))=sigma; solve coefficient-by-coefficient')
print('source: gamma=c/4+g32(theta)c^(3/2)+g2(theta)c^2+...')
print('rho,zeta exact formulas determine C2 once g2 is supplied')
