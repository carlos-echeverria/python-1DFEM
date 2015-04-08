import numpy as np

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

    def __init__(self,a,b,N):
        self.size = N
        self.mesh = np.linspace(a,b,N+1)
        #self.__g  =  
        
    def coordinates(self,n=None):
        if n==None:
            mesh= self.mesh
            return mesh
        else:
            mesh = self.mesh[n]
        return mesh

    def cells(self,n=None):
        if n==None:
            for i in range(self.size - 2):
                cells=np.zeros((self.size,1,2))
                cells[i]=self.mesh[i:i+2] 


            return cells
        else:
            mesh = self.mesh
            cells = mesh[n:n+2]
            return cells




class Mesh2(object):
    """
    A mesh is a list of point global coordinates
    and a list of element definitions (by endpoint).
    This class defines a uniform mesh on an interval.

    Mesh(N,a,b) N intervals on [a,b]
    coordinates(n=None): return coord node n, or all coords
    cells(n=None): array of n-th cell's endpoint numbers, or all
    size(): returns number of elements
    """

    def __init__(self,a,b,N):
        self.size = N
        self.a = a
        self.b = b
        self.N = N

    def coordinates(self,n=None):
        if n==None:
            mesh = np.linspace(self.a,self.b,self.N)
            return mesh
        else:
            mesh =np.linspace(self.a,self.b,self.N)[n]
            return mesh

    def cells(self,n=None):
        if n==None:
            cells = np.linspace(self.a,self.b,self.N)
            return cells
        else:
            cells =np.linspace(self.a,self.b,self.N)[n:n+2]
            return cells



