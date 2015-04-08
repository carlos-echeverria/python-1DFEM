# coding: utf-8
import sympy as sy
import numpy as np.

x=symbols('x')
x=sp.symbols('x')

sy.diff(sy.exp(x))
sy.diff((1+x)*sy.exp(1-x)+x*(1-sy.exp(-x)))
u1=sy.diff((1+x)*sy.exp(1-x)+x*(1-sy.exp(-x)))
sy.diff(u1)
u2=sy.diff(u1)
u=(1+x)*sy.exp(1-x)+x*(1-sy.exp(-x))
u2+2*u1+u
sol=u2+2*u1+u
sy.expand(sol)
u1.subs(x,0)
u1.subs(x,1)
