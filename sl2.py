# Code for the weil representation of SL_2.  This is an action on the space of functions on FF_q, and we express it using the basis of delta functions.
# There are some special elements of SL_2: w = [[0, 1], [-1, 0]], n(u) = [[1, u], [0, 1]] and z(a) = [[a, 0], [0, a^{-1}]].
# The strategy is to represent each of these elements, then factor any arbitrary matrix in terms of them.
# Finally we include some code to check that we've really got a homomorphism.


# work with the field with p elements, modelled as Z/p
p = 5
q = p # code for psi assumes q=p!
F = GF(p, repr='int')
FMat = MatrixSpace(F, 2, 2)
CMat = MatrixSpace(CC, q, q) 
CMatEven = MatrixSpace(CC, floor(p/2)+1, floor(p/2)+1)
CMatOdd = MatrixSpace(CC, floor(p/2), floor(p/2))
gen1 = FMat(matrix([[0, -1], [1, 0]]))
gen2 = FMat(matrix([[1,  1], [0, 1]]))


# fix an additive character psi.   
# When we upgrade to prime powers, we will need to precompose with a trace function.
def psi(a):
	return exp(sage.rings.integer.Integer(a)*2*I*pi/q)

# output rho(n_u) in the delta function basis
def rhoN(u):
	diagonalEntries = [0]*q
	for i in range(0, q):
		diagonalEntries[i] = psi(-u*i*i/2)
	return diagonal_matrix(diagonalEntries)

# output rho(w) in the delta function basis
def rhoW():
	m = matrix(q,q,I-I) # initialize zero matrix
	for i in range(0,q): 
		for j in range(0,q):
			m[i,j] = psi(-i*j)
	return I**((p-1)/2)/sqrt(q)*m

# Let's precompute
WW = rhoW()

# output rho(z_a) in the delta function basis
def rhoZ(a):
	m = matrix(q,q,I-I) # initialize zero matrix
	for i in range(0,q):
		j = a*F(i)    # "F(i)" means i is interpreted in F
		m[i,j] = 1    # m is now the permutation matrix for f(x) -> f(a^(-1)x) in terms of delta functions
	c = kronecker(a,q)    # kronecker is Legendre symbol
	m = c*m  
	return m

# given a 2x2 matrix m, output rho(m)
def rho(m):
	a = m[0,0]
	b = m[0,1]
	c = m[1,0]
	d = m[1,1]
	if c == 0:
		output = rhoZ(a)*rhoN(b/a)
	else:
		output = rhoN(a/c)*WW*rhoZ(-c)*rhoN(d/c)
	return output

# create a random matrix for testing, in the stupidest possible way
def randSL2():
	while(True):
		m = random_matrix(F,2)
		if m.determinant() == 1:
			return m

# Outputs f(A)f(B) - f(A*B) coerced to complex entries.
# f is a function to GLn
def homoTest(A, B, f, n):
	CMatn = MatrixSpace(CC, n, n)
	fA = CMatn(f(A))
	fB = CMatn(f(B))
	fAB = CMatn(f(A*B))
	return fA*fB - fAB

# Test the purported homomorpism f: SL2 --> GLn, blah times, and output the worst result
# If you get two zero matrices, you win.  
def manyHomoTest(blah, f, n):
	exampleA = matrix(F, q,q,0)
	exampleB = matrix(F, q,q,0)
	for i in range(0,blah):
		tempA = randSL2()
		tempB = randSL2()
		diff = homoTest(tempA, tempB, f, n)
		maxDiff = max(max(max(diff)))
		if maxDiff > 10**(-15):
			exampleA = tempA
			exampleB = tempB
			break
	return [exampleA, exampleB]

# evaluate the sum of f(g) as g runs over the two generators and their inverses
def genSum(f):
	return f(gen1) + f(gen2) + f(gen1.inverse()) + f(gen2.inverse())

################################################################################
# In what follows we want to restrict the Weil representation to the subspaces #
# of odd and even functions.  Our basis for the even functions consists of     #
# d_a + d_{p-a}, where d_a is the delta function at a, as a runs from 0 up to  #
# floor(p/2).  Our basis for odd function is d_a - d_{p-a} as a runs from 1 up #
# to floor(p/2).  This will need to be reworked for the prime power case.      #
################################################################################


# Operator projecting onto the even functions by f --> [f(x) + f(-x)]/2
def evenProjector():
	pHalf = floor(p/2)
	output = matrix(QQ, pHalf+1, q, 0)
	output[0,0]=1
	for i in range(1,pHalf+1):
		output[i,i]=0.5
		output[i,q-i]=0.5
	return output

# Operator projecting onto the odd functions by f --> [f(x)-f(-x)]/2
def oddProjector():
	pHalf = floor(p/2)
	output = matrix(QQ, pHalf, q, 0)
	for i in range(0, pHalf):
		output[i,i+1]=0.5
		output[i,q-i-1]=-0.5
	return output

# sections of evenProjector and oddProjector
es = 2*ep.transpose()
es[0,0]=1
os = 2*op.transpose()

def rhoEven(m):
	return ep*rho(m)*es

def rhoOdd(m):
	return op*rho(m)*os
	

