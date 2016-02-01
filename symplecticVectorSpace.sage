# Until I learn how to do better, the intention is that you copy and paste this into a sage notebook.

# Choose a prime q and a dimension d.  For now these are hard-coded as 5 and 1 respectively.  
# We will work with the symplectic vector space of dimension 2d.
q = 5
d = 1

# Create a vector space of dimension 2d over the field with q elements.
V = VectorSpace(GF(q), 2*d)

# Create a matrix M to encode the standard symplectic form on V
hSpace = MatrixSpace(GF(q), d, d)
M1 = (hSpace(0).augment(hSpace(-1)))
M = M1.stack(hSpace(1).augment(hSpace(0)))

# create a function that evaluates the symplectic inner product of a pair of vectors u and v
def omega(u,v): return u.inner_product(M*v)

# create a class for elements of the Heisenberg Group on V
class heisenbergElement:
  # The first line says that when you create an element of the Heisenberg group, 
  # you are supposed to choose an element of V and an element of the center
  def __init__(self, vComponent, zComponent):
    self.v = vComponent
    self.z = zComponent
    
  # Say how to multiply two elements of the heisenberg group.  This is totally untested.
