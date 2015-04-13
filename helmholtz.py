from main import *
from matplotlib import pyplot as plt


"""
Use main.py to solve PDE u''+k*u=f, x in [0,1]
#1. u=1, f=x, Dirichlet b.c. (verification)
#2. u=x, f=1+x, Dirichlet b.c. (verification)
#3. u=x**2, f=x**2+2*x+2, Dirichlet b.c. (verification)
#4. u=(1+x)*exp(1-x)+x*(1-exp(-x)), f=x+2, Neumann conditions

C.Echeverria 13/04/15
"""

def helmholtz(title, N, rhsfn, exactfn):
    """
    generic verification runs for PDE u''+ k*u = rhs on [0,1]
    Dirichlet b.c.
    title = descriptive title
    N = number of elements
    rhsfn = function for rhs as function of x
    exact_sol = function for exact solution as function of x
    CE 13/04/15
    """
    # Create mesh with N elements in [0,1]
    mesh = Mesh(N,0.0,1.0)
    
    # Generate shape functions and function space
    sfns = Shapefns()
    V = FunctionSpace(mesh,sfns)
    
    # rhs and exact
    x = V.dofpts()
    exact = exactfn(x)

    rhs = rhsfn(x)
    b = V.int_phi(rhs)
    
    k = 1
    # assemble stiffness matrix: u''+2*u'+u
    # integration by parts on first term introduces minus sign
    A = -V.int_phi_phi(derivative=[True,True]) \
        +k*V.int_phi_phi(derivative=[False,False]) \
   

    # insert boundary conditions
    # left bndry u=1 (from exact solution)
    A[0,0] = 1.0
    A[0,1:] = 0.0
    b[0] = exact[0]
    # right bndry u=1 (from exact solution)
    A[-1,-1] =1.0
    A[-1,0:-1] = 0.0
    b[-1] = exact[-1]

    # solve
    u = la.solve(A,b)

    # check
    print title," relative error=",la.norm(u-exact)/la.norm(exact)
    plt.plot(u,x)
    plt.show()
    

#verification("Case 1",5,rhsfn=lambda(x):0.0*x+1.0),exactfn=lambda(x):0.0*x+1.0)
helmholtz("Case 1", 5, rhsfn=lambda(x):0.0*x+x, exactfn=lambda(x):0.0*x+x)
#verification("Case 2",5,rhsfn=lambda(x):x+2.0,exactfn=lambda(x):x)
#verification("Case 3", 5, rhsfn=lambda(x):x**2+4*x+2.0, exactfn=lambda(x):x**2)

