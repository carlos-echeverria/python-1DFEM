import numpy as np
#if __name__ == '__main__':

class Mesh(object):
    """
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
        #self.__g  =  

    def size(self):
        """ access mesh size"""
        return self.__size
        
    def coordinates(self,n=None):
        if n==None:
            mesh = self.__mesh
            return mesh
        else:
            mesh = self.__mesh[n-1]
        return mesh

    def cells(self,n=None):

       # g = lambda i: self.mesh[i:i+2]

        if n==None:
           # cells = np.zeros((self.size,2))
           # for i in range(self.size):
           #     cells[i,:] = g(i) 
           # return np.asarray([g(i) for i in range(self.size)])
            return np.arange(1,self.size+1)
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
#snippet of code to test the previous code

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
        sfns is the Shapefuns member
        eltno is this element's number
        endnos is a pair of ints giving the numbers of the endpoints 
               in the mesh
        dofnos is an array of ints giving the numbers of the dofs
        """ 
        # ojo is it mesh.size() or mesh..size
        assert(0 <= eltno < mesh.size()), 'element number is not part of domain'
        self.__eltno = eltno
        endnos = mesh.cells(eltno)
        assert(len(endnos) == 2), 'end nos does not have 2 entries'
        self.__endpts = np.array(mesh.coordinates(endnos))
        self.__numDofs = sfns.size()
        assert(len(dofnos)==sfns.size())
        self.__dofnos = dofnos
        self.__dofpts = np.linspace(self.__endpts[0],self.__endpts[1],self.__numDofs)
        self.__sfns = sfns

        # Gauss points and weights: 3-pts are enough?
        self.__gausspts = np.array(\
                (.112701665379258311482073460022, .5, .887298334620741688517926539978))
        self.__gausswts = np.array( (5.0/18.0, 8.0/18.0, 5.0/18.0) )
        # Generate an array fo shape functions evaluated at Gauss points
        self.__gaussvals = np.empty([self.__numDofs,self.__gausspts.size])
        for n in range(self.__numDofs):
            self.__gaussvals[n,:]=sfns.eval(n,self.__gausspts[:])


        
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

    def ddxl(self,n,x):
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
