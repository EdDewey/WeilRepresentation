# Code for the weil representation of SL_2.  This is an action on the space of functions on FF_q, and we express it using the basis of delta functions.
# There are some special elements of SL_2: w = [[0, 1], [-1, 0]], n(u) = [[1, u], [0, 1]] and z(a) = [[a, 0], [0, a^{-1}]].
# The strategy is to represent each of these elements, then factor any arbitrary matrix in terms of them.
# Finally we include some code to check that we've really got a homomorphism.


# work with the field with p elements, modelled as Z/p
p = 3
q = p # code for psi assumes q=p!
F = GF(p, repr='int')
FMat = MatrixSpace(F, 2, 2)
CMat = MatrixSpace(CC, q, q)
gen1 = FMat(matrix([[0, -1], [1, 0]]))
gen2 = FMat(matrix([[1,  1], [0, 1]]))


# fix an additive character psi.   
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

# output rho(z_a) in the delta function basis
def rhoZ(a):
	m = matrix(q,q,I-I) # initialize zero matrix
	for i in range(0,q):
		j = a*F(i)    # "F(i)" means i is interpreted in F
		m[i,j] = 1    # m is now the permutation matrix for f(x) -> f(a^(-1)x) in terms of delta functions
	c = kronecker(a,q)
	m = c*m  
	return m

# given a 2x2 matrix m, output rho(m)
def rho(m):
	a = m[0,0]
	b = m[0,1]
	c = m[1,0]
	d = m[1,1]
	rw = rhoW()
	if a == 0:
		output = rhoZ(b)*rw*rhoN(-b*d)
	else:
		output = rhoZ(a)*rw.inverse()*rhoN(-a*c)*rw*rhoN(b/a)
	return output

# create a random matrix for testing, in the stupidest possible way
def randSL2():
	while(True):
		m = random_matrix(F,2)
		if m.determinant() == 1:
			return m

# Outputs rho(A)rho(B) - rho(A*B) coerced to complex entries
def homoTest(A, B):
	rhoA = CMat(rho(A))
	rhoB = CMat(rho(B))
	rhoAB = CMat(rho(A*B))
	return rhoA*rhoB - rhoAB

# Test n times, output the worst result
def manyHomoTest(blah):
	exampleA = matrix(F, q,q,0)
	exampleB = matrix(F, q,q,0)
	for i in range(0,blah):
		tempA = randSL2()
		tempB = randSL2()
		diff = homoTest(tempA, tempB)
		maxDiff = max(max(max(diff)))
		if maxDiff > 10**(-15):
			exampleA = tempA
			exampleB = tempB
			break
	return [exampleA, exampleB]

# evaluate the rum of rho(g) as g runs over the two generators and their inverses
def genSum():
	return rho(gen1) + rho(gen2) + rho(gen1.inverse()) + rho(gen2.inverse())
