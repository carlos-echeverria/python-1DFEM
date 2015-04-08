from dolfin import *
mesh = UnitSquareMesh( 32,32)

V  = FunctionSpace(mesh,"Lagrange", 1)
u  = TrialFunction(V)
v  = TestFunction(V)
f  = Constant(-6.0)
g  = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]")
bc = DirichletBC(V, g, DomainBoundary())
a  = inner(grad(u), grad(v))*dx
L  = f*v*dx
u  = Function(V)

solve(a == L, u, bc)
plot(u, interactive=True)
