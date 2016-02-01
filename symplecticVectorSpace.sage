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

# create a function that evaluates the symplectic inner product of a pair of vectors
