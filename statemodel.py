from math import factorial

class StateModel:
	def __init__(self, n, N):
		self.n = n				# number of quasiholes
		self.N = N				# number of particles
		self.L = 0				# number of orbitals / flux quanta

		self.Fmax = 0			# Max value of partitioning parameter
		self.Fmin = 0			# Min value of partitioning parameter

		self.cell_size = 0		# Valid states have at most max_continguous particles in cell_size sites.
		self.max_contiguous = 0	# 	...As a check, max_continguous/cell_size should = filling factor nu

		self.max_occupation = 0	# Each site can have at most max_occupation particles.

		self.geometry=''		# A string holding the geometry, i.e. 'torus' or 'sphere'


	def binom(self,a,b):
		return factorial(a)/ (factorial(b)*factorial(a-b))

	def count(self,F):
		return int(self.pos_degen(F)*self.internal_degen(F))

	def pos_degen(self,F):
		return 1
	
	def pos_degen_string(self,F):
		return ''

	def internal_degen(self,F):
		return 1

	def internal_degen_string(self,F):
		return ''

class GaffnianTorus(StateModel):
	def __init__(self,n,N):
		self.n = n
		self.N = N
		self.L = (3*N+n)/2

		self.Fmax = min(N,n)
		self.Fmin = n%2

		self.cell_size = 3
		self.max_contiguous = 2

		self.max_occupation = 2

		self.geometry = 'torus'
	
	def pos_degen(self,F):
		return float(self.L)/float( (self.N-F)/2+self.n )*self.binom((self.N-F)/2+self.n,self.n)
	
	def pos_degen_string(self,F):
		return '{0}/{1}*({1} choose {2})'.format( self.L, (self.N-F)/2+self.n, self.n )

	def internal_degen(self,F):
		return float( self.n )/float( (self.n+F)/2 )*self.binom( (self.n+F)/2,F )

	def internal_degen_string(self,F):
		return '{0}/{1}*({1} choose {2})'.format( self.n, (self.n+F)/2, F )

class GaffnianSphere(StateModel):
	def __init__(self,n,N):
		self.n = n
		self.N = N
		self.L = (3*N+n)/2-2

		self.Fmax = min(N,n-2)
		self.Fmin = n%2

		self.cell_size = 3
		self.max_contiguous = 2

		self.max_occupation = 2

		self.geometry = 'sphere'
	
	def pos_degen(self,F):
		return self.binom( (self.N-F)/2+self.n, self.n )
	
	def pos_degen_string(self,F):
		return '({} choose {})'.format( (self.N-F)/2+self.n, self.n )

	def internal_degen(self,F):
		return self.binom( (self.n+F)/2-1, F )

	def internal_degen_string(self,F):
		return '({} choose {})'.format( (self.n+F)/2-1, F )