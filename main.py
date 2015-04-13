import numpy as np
import scipy.linalg as la
#if __name__ == '__main__':

class Mesh(object):
    """
    "hola"
    A mesh is a list of point global coordinates
    and a list of element definitions (by endpoint).
    This class defines a uniform mesh on an interval.

    Mesh(N,a,b) N intervals on [a,b]
    coordinates(n=None): return coord node n, or all coords
    cells(n=None): array of n-th cell's endpoint numbers, or all
    size(): returns number of elements
    """

    def __init__(self,N,a,b):
        self.__size = N
        self.__mesh = np.linspace(a,b,N+1)
        #self.variable  = 10.  

    def size(self):
        """ access mesh size"""
        return self.__size
        
    def coordinates(self,n=None):
        if n==None:
            return self.__mesh
        else:
            return self.__mesh[n]

    def cells(self,n=None):
        assert(0 <= n < self.__size or n is None), 'cell index out of range'
       # g = lambda i: self.mesh[i:i+2]
        if n==None:
           # cells = np.zeros((self.size,2))
           # for i in range(self.size):
           #     cells[i,:] = g(i) 
           # return np.asarray([g(i) for i in range(self.size)])
            return np.arange(self.__size)
        else:
            #return g()
            return np.array([n,n+1])
     



class Shapefns(object):
    """
    Define Quadratic Lagrange shape functions
    These will be defined on the (local) interval [0,1], with mid point 0.5
    Shapefns()
    eval(n,xi): phi[n](xi)
    ddx(n,xi):  dphi[n](xi)
    size(): number of nodes for these shape functions
    """

    def __init__(self):
        """ 
        an array of functions for phi and deriv phi
        """

        self.__phi = [lambda xi: 2.0 * (xi - 0.5) * (xi - 1.0) ,\
                     lambda xi: 4.0 * xi * (1.0 - xi), \
                     lambda xi: 2.0 * xi * (xi - 0.5)]
        
        self.__phi2 = [lambda xi: 1.0-xi,\
                       lambda xi: np.piecewise(xi,[xi<0.5,xi>=0.5],\
                       [lambda y: 2*y, lambda y:2.0*(1-y)]),\
                       lambda xi: xi]
        # and dphi (derivative of phi w.r.t. xi)
        # derivative of second factor * first + derivative of first factor * second
        self.__dphi = [lambda xi: 2.0 * (xi - 0.5) + 2.0*(xi - 1.0) ,\
                      lambda xi: -4.0 * xi + 4.0*(1.0 - xi), \
                      lambda xi: 2.0 * xi + 2.0*(xi - 0.5)]

        self.__dphi2 = [lambda xi: -1.0,\
               lambda xi: np.piecewise(xi,[xi<0.5,xi>=0.495],[lambda y: 2.0,lambda y:-2.0]),\
                lambda xi:  1.0]

        self.__N = 3  #number of nodes in quadratic Lagrange polynomial



    def eval(self,n,xi,linear=False):
        """
        the function phi[n](xi), for any xi
        """
        if linear:
           
            return self.__phi2[n](np.array([xi]))[0]
        else:
            return self.__phi[n](xi)

    def ddx(self,n,xi,linear=False):
        """
        the function dphi[n](xi), for any xi
        """
        if linear:
            return self.__dphi2[n](xi)
        else:
            return self.__dphi[n](xi)

    def size(self,linear=False):
        """
        the  number of points
        """
        if linear:
            return self.__N-1
        else:
            return self.__N

"""
# Snippet of code to test the construction of shapefunctions

N = 10
rightpt = 1.0
print "\n\nFEM1D.py Test case, dx=",rightpt/N

mesh = Mesh(N,0.0,rightpt)
coords = mesh.coordinates()
print "mesh.coordinates()=", coords

sfns = Shapefns()
print "sfns.size()-2=", sfns.size()-2
xi = np.linspace(0,1,11)
import matplotlib.pyplot as plt
if True:
    for n in range(3):
        plt.plot(xi,sfns.eval(n,xi,True))
        plt.axis([0,1,0,1.5])
        plt.show()

        plt.plot(xi,sfns.ddx(n,xi))
        plt.show()
"""      

class FiniteElement(object):
    def __init__(self,mesh,sfns,eltno,dofnos):
        """
        mesh is the mesh it is built on 
        sfns is a member of Shapefuns class
        eltno is this element's number
        endnos is a pair of ints giving the numbers of the endpoints 
               in the mesh
        dofnos is an array of ints giving the numbers of the dofs
        """ 
        assert(0 <= eltno < mesh.size()), 'element number is not part of domain'
        self.__eltno = eltno
        endnos = mesh.cells(eltno)
        assert(len(endnos) == 2), 'endnos does not have 2 entries'
        self.__endpts = np.array(mesh.coordinates(endnos))
        self.__numDofs = sfns.size()
        assert(len(dofnos)==sfns.size())
        self.__dofnos = dofnos
        self.__dofpts = np.linspace(self.__endpts[0],self.__endpts[1],self.__numDofs)
        self.__sfns = sfns

        # Gauss points and weights: 3-pts are enough?
        self.__gausspts = np.array( \
                (.112701665379258311482073460022, .5, .887298334620741688517926539978))
        self.__gausswts = np.array( (5.0/18.0, 8.0/18.0, 5.0/18.0) )
        # Generate an array for shape functions evaluated at Gauss points
        self.__gaussvals = np.empty([self.__numDofs,self.__gausspts.size])
        for n in range(self.__numDofs):
            self.__gaussvals[n,:]=sfns.eval(n,self.__gausspts[:])
        # can we pass the entire array to avoid the for loop?

        
    def endpts(self):
        """ access endpoints """
        return self.__endpts
    
    def dofpts(self):
        """ access dof points """
        return self.__dofpts

    def dofnos(self):
        """ access dof point numbers """
        return self.__dofnos

    def numDofs(self):
        """ access numDofs """
        return self.__numDofs

    def eval(self,n,x):
        """
        evalueate the n-th shape function on this element at the
        spatial coordinate x
        """
        # map x to xi
        xx = np.array(x)
        xi = (xx - self.__endpts[0]) / (self.__endpts[1] - self.__endpts[0])
        # evaluate
        return self.__sfns.eval(n,xi) * (xi >=0.) * (xi <= 1.)

    def ddx(self,n,x):
        """
        evaluate the n-th shape function derivative on this element at the
        spatial coordinate x
        """
        # map x to xi
        xi = (x - self.__endpts[0]) / (self.__endpts[1] - self.__endpts[0])
        # evaluate
        return self.__sfns.ddx(n,xi) * (xi>=0.) * (xi <= 1.0)

    def integral(self, f1=None, f2=None, derivative=False):
        """
        Integrate either phi[i](xi)*f1(xi)*f2(xi)  or dphi[i]*f1*f2 over
        this element, depending on if derivative is False or True
        Returns a vector of 3 results, 
        one for phi[0], one for phi[1], and one for phi[2].
        f1 and f2 are assumed to have been mapped to this elment as arrays.
        If derivative is True, phi is replaced with dphi
        """
        L = self.__endpts[1]-self.__endpts[0] # length of element
        t = self.__gausswts.copy()
        gp = self.__gausspts

        if f1 != None:
            assert(len(f1)==self.__numDofs), 'array dimension missmatch f1'
            fvals = np.zeros([self.__gausspts.size])
            for n in range(self.__numDofs):
                fvals += f1[n]*self.__gaussvals[n,:]
            t *= fvals

        if f2 != None:
            assert(len(f2)==self.__numDofs), 'array dimension missmatch f2'
            fvals = np.zeros([self.__gausspts.size])
            for n in range(self.__numDofs):
                fvals += f2[n]*self.__gaussvals[n,:]
            t *= fvals

        if derivative:
            # really: t *= L*(1/L)
            q = np.dot(np.array([self.__sfns.ddx(0,gp), \
                                 self.__sfns.ddx(1,gp), \
                                 self.__sfns.ddx(2,gp)]),t)
        else:
            t *= L # correct for affine map x->xi
            q = np.dot(np.array([self.__sfns.eval(0,gp), \
                                 self.__sfns.eval(1,gp), \
                                 self.__sfns.eval(2,gp)]),t)
        
        return q



"""
# snippet of code to test the Finite Element class
N = 5
rightpt = 5.0
print "\n\nFEM1D.py Test case, dx=",rightpt/N

mesh = Mesh(N,0.0,rightpt)
coords = mesh.coordinates()
print "mesh.coordinates()=", coords

sfns = Shapefns()
#print "sfns.size()-2=", sfns.size()-2

elt = FiniteElement(mesh,sfns,0,[0,1,2])

if N == 5 and np.abs(coords[-1]-5.0) < 1.e-10:
    #test some integrals
    print "elt integral() err=", max(abs(elt.integral()-[1./6,2./3,1./6]))
    print "integral(deriv) err=", \
            max(abs(elt.integral(derivative=True)-[-1,0,1]))

    #pick the function f(x)=x, find its expansion coefficients
    ex = np.empty([sfns.size()])
    ex[0] = elt.endpts()[0]
    ex[2] = elt.endpts()[1]
    ex[1] = 0.5*(ex[0]+ex[2])
    ex2 = ex**2
    print "integral(x) err=",    \
            max(abs(elt.integral(ex)-[0,1./3,1./6]))
    print "integral(x**2) err=", \
            max(abs(elt.integral(ex2)-[-1./60,1./5,3./20]))
    print "integral(x**2) err=", \
            max(abs(elt.integral(ex,ex)-[-1./60,1./5,3./20]))
    print "integral(x,phi') err=", \
            max(abs(elt.integral(ex,derivative=True)-[-1./6,-2./3,5./6]))
    print "integral(x**2,phi') err=", \
            max(abs(elt.integral(ex2,derivative=True)-[0,-2./3,2./3]))
"""

class FunctionSpace(object):
    """
    A FunctionSpace has a list of elements numbered and with coords according 
    to mesh FunctionSpace(mesh,sfns): constructor, sfns is ShapeFuns
    size(): number of elements
    ndofs(): number of all dofs
    dofpts(n=None): coordinates of dof[n] or all dofs
    int_phi_phi(c=None,derivative=[False,False]):
        integral(c*phi*phi) or
        integral(c*dphi*phi) or
        integral(c*dphi*dphi) or
        integral(c*phi*dphi)
    int_phi(f=None,derivative=False):
        integral(f*phi) or
        integral(f*dphi)
    """
    def __init__(self, mesh,sfns):
        """
        mesh is a Mesh class
        sfns is a Shapefns class
        """
        self.__size=mesh.size()
        # number the elements in the same was as mesh
        self.__elts = list([])
        self.__dofpts = list([])
        self.__nDOFs=0
        for n in range(self.__size):
            #ASSUMing only boundary points are number 0 and (self.__size)
            if n==0:
                self.__nDOFs += 3
                dofs=[2*n, 2*n+1,2*n+2]
                newdofs = range(3)
            else:
                self.__nDOFs +=2
                dofs=[2*n, 2*n+1, 2*n+2]
                newdofs = range(1,3)
            fe = FiniteElement(mesh,sfns,n,dofs)
            self.__elts.append(fe)
            for i in newdofs:
                self.__dofpts.append(fe.dofpts()[i])
        self.__dofpts = np.array(self.__dofpts)

    def size(self):
        return len(self.__elts)

    def Ndofs(self):
        return self.__nDOFs

    def dofpts(self,n=None):
        if n == None:
            return self.__dofpts
        else:
            return self.__dofpts[n]
    
    def int_phi_phi(self, c=None, derivative=[False,False]):
        """
        assemble $\int c(x)\phi(x) dx$ or with $d\phi/dx$
        """
        A = np.zeros([self.__nDOFs, self.__nDOFs])
        # loop over elements
        for elt in self.__elts:
            d0=elt.dofnos()
            if c != None:
                cc = c[d0]
            else:
                cc = None
            N = elt.numDofs()
            endpts=elt.endpts()
            L = endpts[1]-endpts[0] # length of element
            for j in range(N):
                if derivative[1]:
                    #chain rule: d(xi)/d(x) = 1/L
                    phi = elt.ddx(j,elt.dofpts())/L
                else:
                    phi = elt.eval(j,elt.dofpts())
                A[d0,d0[j]] += elt.integral(phi,cc,derivative=derivative[0])
        return A

    def int_phi(self, f= None, derivative=False):
        """
        assemble $\int f(x)\phi(x) dx$ or with $d\phi/dx$
        """
        b = np.zeros([self.__nDOFs,1])
        #loop over elements
        for elt in self.__elts:
            d0=elt.dofnos()
            if f != None:
                ff = f[d0]
            else:
                ff = None
            N = elt.numDofs()
            endpts=elt.endpts()
            L = endpts[1]-endpts[0] # length of element
            for j in range(N):
                if derivative:
                    #chain rule: d(xi)/d(x) = 1/L
                    phi = elt.ddx(j,elt.dofpts())/L
                else:
                    phi = elt.eval(j,elt.dofpts())
                b[d0,0] += elt.integral(phi,ff,derivative=derivative)
        return b

"""
snippet of code to test the Function Space class
"""
if __name__ == '__main__':

    N = 5
    rightpt = 5.0
    print "\n\nFEM1Dclasses.py Test case, dx=",rightpt/N
    mesh = Mesh(N,0.0,rightpt)
    coords = mesh.coordinates()
    if N==5:
        print "mesh.coordinates()=",coords

    # generate an element
    sfns = Shapefns()
    print "sfns.size()-3=",sfns.size()-3
    xi = np.linspace(0,1,100)
    import matplotlib.pyplot as plt
    if False:
        for n in range(3):
            plt.plot(xi,sfns.eval(n,xi))
            plt.show()
            plt.plot(xi,sfns.ddx(n,xi))
            plt.show()

    elt = FiniteElement(mesh,sfns,0,[0,1,2])

    if N==5 and np.abs(coords[-1]-5.0) < 1.e-10:
        # test some integrals
        print "elt integral() err=",max(abs(elt.integral()-[1./6,2./3,1./6]))
        print "integral(deriv) err=",max(abs(elt.integral(derivative=True)- [-1,0,1]))

        # test some more integrals
        # pick the function f(x)=x, find its expansion coefs
        ex = np.empty([sfns.size()])
        ex[0] = elt.endpts()[0]
        ex[2] = elt.endpts()[1]
        ex[1] = .5*(ex[0]+ex[2])
        ex2 = ex**2
        print "integral(x) err=",max(abs(elt.integral(ex)-[0,1./3,1./6]))
        print "integral(x**2) err=",max(abs(elt.integral(ex2)-[-1./60,1./5,3./20]))
        print "integral(x**2) err=",max(abs(elt.integral(ex,ex)-[-1./60,1./5,3./20]))

        print "integral(x,phi') err=",\
              max(abs(elt.integral(ex,derivative=True)-[-1./6,-2./3,5./6]))
        print "integral(x**2,phi') err=",\
              max(abs(elt.integral(ex2,derivative=True)-[0,-2./3,2./3]))
        print

    V = FunctionSpace(mesh,sfns)
    print "V.Ndofs()-correct=",V.Ndofs()-(2*N+1)
    print "V.size()-correct=",V.size()-N

    x = V.dofpts()
    f = x.copy()
    print "error in integral x over [",x[0],",",x[-1],"]=",\
        np.sum(V.int_phi(f))-x[-1]**2/2.
    f = 0.0*x+1
    print "error in integral 1 over [",x[0],",",x[-1],"]=",\
        np.sum(V.int_phi(f))-x[-1]
    f = x.copy()**2
    print "error in integral x**2 over [",x[0],",",x[-1],"]=",\
        np.sum(V.int_phi(f))-x[-1]**3/3.
    f = x.copy()**3
    print "error in integral x**3 over [",x[0],",",x[-1],"]=",\
        np.sum(V.int_phi(f))-x[-1]**4/4.
    f = x.copy()**4
    print "error in integral x**4 over [",x[0],",",x[-1],"]=",\
        np.sum(V.int_phi(f))-x[-1]**5/5.," should be nonzero."
    print "norm(V.dofpts()-correct)=",\
        la.norm(V.dofpts()-np.linspace(0,coords[-1],2*N+1))

    # check eq \int \phi \phi * u = 1 gives 1 back (homog Neumann b.c.)
    # generate matrix by assembling $A_{ij}=\int\phi_i\phi_j$
    A = V.int_phi_phi()
    if N ==5 and np.abs(coords[-1]-5.0) < 1.e-10:
        print "error A00=",A[0,0]-2./15.
        print "error A01=",A[0,1]-1./15.
        print "error A02=",A[0,2]+1./30.
        print "error A11=",A[1,1]-8./15.
        print "error A12=",A[1,2]-1./15.
        print
        print "error A22=",A[2,2]-4./15.
        print "error A23=",A[2,3]-1./15.
        print "error A24=",A[2,4]+1./30.
        print "error A33=",A[3,3]-8./15.
        print "error A34=",A[3,4]-1./15.
        print

    Ndofs = V.Ndofs()
    f = np.random.rand(Ndofs)
    b = V.int_phi(f)
    print "norm(A*f-b)=",la.norm(np.dot(A.transpose(),f)-b)

    #trivial check with coefficient
    c = np.ones([Ndofs])
    A1 = V.int_phi_phi(c)
    print "Norm difference matrices=",la.norm(A-A1)

    # try putting a coefficient in: c(x) = (1+x)
    # rhs = ones, soln = 1/(1+x)
    c = (1.0+x)
    B = V.int_phi_phi(c)
    if N == 5 and np.abs(coords[-1]-5.0) < 1.e-10:
        print "error B00=",B[0,0]-3./20.
        print "error B01=",B[0,1]-1./15.
        print "error B02=",B[0,2]+1./20.
        print "error B11=",B[1,1]-12./15.
        print "error B12=",B[1,2]-2./15.
        print
        print "error B22=",B[2,2]-8./15.
        print "error B23=",B[2,3]-2./15.
        print "error B24=",B[2,4]+1./12.
        print "error B33=",B[3,3]-4./3.
        print "error B34=",B[3,4]-3./15.

    C = V.int_phi_phi(derivative=[True,True])
    if N == 5 and np.abs(coords[-1]-5.0) < 1.e-10:
        print "\n Laplace Matrix"
        print "error C00*3=",C[0,0]-7./3.
        print "error C01*3=",C[0,1]+8./3.
        print "error C02*3=",C[0,2]-1./3.
        print "error C11*3=",C[1,1]-16./3.
        print "error C12*3=",C[1,2]+8./3.
        print
        print "error C22*3=",C[2,2]-14./3.
        print "error C23*3=",C[2,3]+8./3.
        print "error C24*3=",C[2,4]-1./3.
        print "error C33*3=",C[3,3]-16./3.
        print "error C34*3=",C[3,4]+8./3.
        print

    soln2 = np.ones([Ndofs])
    b2 = np.dot(C,soln2)
    print "const soln Laplace, norm check=",la.norm(b2)

    soln = x
    b0 = np.dot(C,soln)
    rhs0 = V.int_phi(np.zeros([Ndofs]))
    # natural b.c. not satisfied, don't check them
    rhs0[0] = -b0[0]
    rhs0[-1] = -b0[-1]
    print "soln=x Laplace, norm check=",la.norm(rhs0+b0)

    soln = x**2
    b1 = np.dot(C,soln)
    rhs1 = V.int_phi(2.0*np.ones([Ndofs]))
    # natural b.c. not satisfied on right, don't check it 
    rhs1[-1] = -b1[-1]
    print "soln=x**2 Laplace, norm check=",la.norm(rhs1+b1)

    D = V.int_phi_phi(derivative=[False,True])
    soln = np.ones([V.Ndofs()])
    b2 = np.dot(D,soln)
    print "norm check (rhs d/dx+Neumann, const soln)=",la.norm(b2)

    D[0,0] = 1.0
    D[0,1:] = 0.0
    D[-1,-1] = 1.0
    D[-1,0:-1] = 0.0
    soln = x
    b3 = np.dot(D,soln)
    rhs3 = V.int_phi(np.ones([Ndofs]))
    rhs3[0] = soln[0]
    rhs3[-1] = soln[-1]
    print "norm check (d/dx+Dirichlet soln=x)=",la.norm(rhs3-b3) 
