import math
import json

class QuadExt:
	# The coefficients of the minimal polynomial of w
	# So that w is a root of :	 X^2 - TR_W * X + DET_W = 0
	TR_W=None  # Sum of w and its conjugate
	DET_W=None  # Product of w and its conjugate
	
	# Integer generating ring of integers
	# That is d such that  Q[sqrt(d)] ~= Z[w]
	DVALUE=None
	
	def __init__(self, x, y):
		self.x = x
		self.y = y
	
	@staticmethod
	def setd(d):
		if d % 4 == 1:
			QuadExt.TR_W, QuadExt.DET_W = 1, (1 - d) // 4
		else:
			QuadExt.TR_W, QuadExt.DET_W = 0, -d
		
		QuadExt.DVALUE = d
	
	def __eq__(self, other):
		if isinstance(other, QuadExt):
			return self.x == other.x and self.y == other.y
		else:
		 	# Assume other is a real
		 	return self.x == other and self.y == 0
	
	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __hash__(self):
		return hash((self.x, self.y))
	
	
	
	def __neg__(self):
		return QuadExt(-self.x, -self.y)
	
	def __add__(self, other):
		if isinstance(other, QuadExt):
			return QuadExt(self.x + other.x, self.y + other.y)
		else:
			# Assume other is a real value
			return QuadExt(self.x.__add__(other), self.y)
	
	def __radd__(self, other):
		if isinstance(other, QuadExt):
			return QuadExt(other.x + self.x, other.y + self.y)
		else:
			# Assume other is a real value
			return QuadExt(self.x.__radd__(other), self.y)
	
	def __sub__(self, other):
		if isinstance(other, QuadExt):
			return QuadExt(self.x - other.x, self.y - other.y)
		else:
			# Assume other is a real value
			return QuadExt(self.x.__sub__(other), self.y)
	
	def __rsub__(self, other):
		if isinstance(other, QuadExt):
			return QuadExt(other.x - self.x, other.y - self.y)
		else:
			# Assume other is a real value
			return QuadExt(self.x.__rsub__(other), -self.y)
	
	
	def __mul__(self, other):	
		if isinstance(other, QuadExt):
			yy = self.y * other.y
			return QuadExt(self.x * other.x - yy * QuadExt.DET_W, self.x * other.y + self.y * other.x + yy * QuadExt.TR_W)
		else:
			# Assume other is a real value
			return QuadExt(self.x.__mul__(other), self.y.__mul__(other))
	
	def __rmul__(self, other):
		if isinstance(other, QuadExt):
			yy = other.y * self.y
			return QuadExt(other.x * self.x - yy * QuadExt.DET_W, other.y * self.x + other.x * self.y + yy * QuadExt.TR_W)
		else:
			# Assume other is a real value
			return QuadExt(self.x.__rmul__(other), self.y.__rmul__(other))
	
	def __div__(self, other):
		return self.__truediv__(other)
	
	def __truediv__(self, other):
		if isinstance(other, QuadExt):
			# Multiply `self` by the conjugate of `other`
			prod = self * other.conjugate()
			
			# Find the squared magnitude of `other`
			other_norm = abs(other)
			return QuadExt(prod.x / other_norm, prod.y / other_norm)
		else:
			# Assume other is a real scalar
			return QuadExt(self.x.__truediv__(other), self.y.__truediv__(other))
	
	def __floordiv__(self, other):
		if isinstance(other, QuadExt):
			# Multiply `self` by the conjugate of `other`
			prod = self * other.conjugate()
			
			# Find the square magnitude of `other`
			other_norm = abs(other)
			return QuadExt(prod.x // other_norm, prod.y // other_norm)
		else:
			# Check for divisibility
			# Assume other is a real scalar
			return QuadExt(self.x.__floordiv__(other), self.y.__floordiv__(other))
	
	def __mod__(self, other):
		if isinstance(other, QuadExt):
			# Find quotient
			quot = self // other
			# Remove quotient from `self` leaving residue
			return self - quot * other
		else:
			return QuadExt(self.x.__mod__(other), self.y.__mod__(other))
	
	
	# Return the magnitude of the complex number squared
	def __abs__(self):
		return self.x * self.x + self.x * self.y * QuadExt.TR_W + self.y * self.y * QuadExt.DET_W
	
	def conjugate(self):
		return QuadExt(self.x + self.y * QuadExt.TR_W, -self.y)
	
	def to_frac(self):
		""" Map element into the field of fractions of O_d """
		from fractions import Fraction
		return QuadExt(Fraction(self.x), Fraction(self.y))
	
	def to_point(self):
		return (float(self.x + self.y * QuadExt.TR_W / 2), float(self.y * math.sqrt(4 * QuadExt.DET_W - QuadExt.TR_W) / 2))
	
	
	def __str__(self):
		if self.y == 0:
			return str(self.x)
		elif self.y < 0:
			if self.x == 0:
				return str(self.y) + 'w'
			else:
				return str(self.x) + ' - ' + str(-self.y) + 'w'
		else:
			if self.x == 0:
				return str(self.y) + 'w'
			else:
				return str(self.x) + ' + ' + str(self.y) + 'w'
	
	def __repr__(self):
		return "QuadExt(%s, %s)" % (self.x, self.y)
	
	def to_json(self):
		return [self.x, self.y]
	
	@staticmethod
	def from_json(xy):
		return QuadExt(xy[0], xy[1])
	
	
	# Returns a set of representatives for the residue classes of O_d mod `self`
	def residues(self, inc_edge=False):
		"""
		Get the residues
		
		If `inc_edge` is true then those rationals lying
		on the edge of the fundamental domain are included
		"""
		
		# For zero don't do any calculations
		if self.x == 0 and self.y == 0:
			return set()
		
		# Get minimum and maximum bounds on y-value of residue representatives
		# Consider the vertices of the fundamental domain multiplied by self
		vertYs = [self.y, (self * QuadExt(0, 1)).y, (self * QuadExt(1, 1)).y]
		# The plus and minus one is necessary since
		# none of [self * 1, self * w, self * (1 + w)] are residues
		# Unless inc_edge is true
		if inc_edge:
			minY, maxY = min(0, min(vertYs)), max(0, max(vertYs))
		else:
			minY, maxY = min(0, min(vertYs) + 1), max(0, max(vertYs) - 1)
		
		mag2 = abs(self)
		resids = set()
		for y in range(minY, maxY + 1):
			# Get `x` values that are valid
			
			# With M being the squared magnitude of self
			# Should fulfill   self.x * y - M < self.y * x <= self.x * y
			# Should fulfill   -y * self.y * @@detW <= (self.x + self.y * @@trW) * x < M - y * self.y * @@detW
			# If inc_edge is true then the strict inequalities are weakened
			
			# Get the first pair of bounds for the first condition
			try:
				b1, b2 = (self.x * y - mag2) / self.y, (self.x * y) / self.y
				
				if inc_edge:
					# Swap if wrong order
					if b2 < b1:
						b1, b2 = b2, b1
					
					b1, b2 = math.ceil(b1), math.floor(b2)
				else:
					# Swap if wrong order
					# The use of both ceil and floor comes from the fact that b1 is exclusive and b2 is inclusive
					if b2 < b1:
						b1, b2 = math.ceil(b2), math.ceil(b1) - 1
					else:
						b1, b2 = math.floor(b1) + 1, math.floor(b2)
				
			except ZeroDivisionError:
				# When the bounds don't apply any constraints on `x`
				b1, b2 = None, None
			
			
			# Get the second pair of bounds for the second condition
			try:
				den = self.x + self.y * QuadExt.TR_W
				c1, c2 = (-y * self.y * QuadExt.DET_W) / den, (mag2 - y * self.y * QuadExt.DET_W) / den
				
				if inc_edge:
					# Swap if wrong order
					if c2 < c1:
						c1, c2 = c2, c1
					
					c1, c2 = math.ceil(c1), math.floor(c2)
				else:
					# Swap if wrong order
					# The use of both ceil and floor comes from the fact that c1 is inclusive and c2 is exclusive
					if c2 < c1:
						c1, c2 = math.floor(c2) + 1, math.floor(c1)
					else:
						c1, c2 = math.ceil(c1), math.ceil(c2) - 1
				
			except ZeroDivisionError:
				c1, c2 = None, None
			
			# Combine intervals
			if b1 is None:
				lower, upper = c1, c2
			elif c1 is None:
				lower, upper = b1, b2
			else:
				lower, upper = max(b1, c1), min(b2, c2)
			resids |= set(QuadExt(x, y) for x in range(lower, upper + 1))
		
		resids = sorted(resids, key=lambda r: abs(r))  # Order residues by norm
		return resids
	
	@staticmethod
	def within_norm(max_norm):
		"""
		Find all element of Quadratic Extension whose norm is less than or equal to max_norm
		"""
		
		# M / (1 - trW / (4 * detW)) = M * detW / (detW - trW / 4)  serves as an upper bound on x^2
		max_x2 = math.ceil(float(max_norm) * QuadExt.DET_W / (QuadExt.DET_W - QuadExt.TR_W / 4))
		# M / (detW - trW / 4)  serves as an upper bound on y^2
		max_y2 = math.ceil(float(max_norm) / (QuadExt.DET_W - QuadExt.TR_W / 4))
		
		elems = set()
		x = 0
		while x * x <= max_x2:
			y = 0
			while y * y <= max_y2:
				# Check each possible combination of signs
				q1 = QuadExt(x, y)
				if abs(q1) <= max_norm:
					elems.add(q1)
				
				# If `x` is zero then -(x + yw) = x - yw so the next one is redundant
				if x == 0:
					y += 1
					continue
				
				q2 = QuadExt(x, -y)
				if abs(q2) <= max_norm:
					elems.add(q2)
				
				y += 1
			x += 1
		
		elems = sorted(elems, key=lambda elm: abs(elm)) # Order by norm
		return elems
	
	def matrices(self):
		# Get representatives for O_d / self * O_d
		resids = self.residues(inc_edge=True)
		
		# Set of matrices
		mats = []
		
		# Iterate over all possible `a` values
		for a in resids:
			# Iterate over all possible `d` values
			bds = [((a * d - 1) // self, d) for d in resids if (a * d - 1) % self == 0]
			if len(bds) > 0:
				(b, d) = min(bds, key=lambda bd: 100 * abs(bd[0]) + abs(bd[1]))
				# Construct matrix and add to list
				mats.append(Matrix(a, b, self, d))
		
		return mats


class Matrix:
	def __init__(self, a, b, c, d):
		self.a, self.b, self.c, self.d = a, b, c, d
	
	@staticmethod
	def identity():
		return Matrix(1, 0, 0, 1)
	
	def __eq__(self, other):
		return isinstance(other, Matrix) and self.a == other.a and self.b == other.b and self.c == other.c and self.d == other.d
	
	
	
	# Arithmetic Operations
	def __add__(other):
		if isinstance(other, Matrix):
			return Matrix(self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)
		else:
			raise ValueError("Cannot add matrix to %s" % other)
	
	def __sub__(self, other):
		if isinstance(other, Matrix):
			return Matrix(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d)
		else:
			raise ValueError("Cannot subtract %s from matrix" % other)
	
	@staticmethod
	def shift(shf):
		"""
		Create Matrix transform which adds `shf` to its input
		   z  --->  z + shf
		"""
		return Matrix(1, shf, 0, 1)
	
	def to_shift(self):
		"""
		Convert shift matrix into its shift
		Raises ValueError if given matrix isn't a shift
		"""
		
		if self.a == 1 and self.d == 1 and self.c == 0:
			return self.b
		elif self.a == -1 and self.d == -1 and self.c == 0:
			return -self.b
		else:
			raise ValueError("%s Isn't a shift matrix" % self)
	
	def __mul__(self, other):
		if isinstance(other, Matrix):
			return Matrix(
				self.a * other.a + self.b * other.c,
				self.a * other.b + self.b * other.d,
				self.c * other.a + self.d * other.c,
				self.c * other.b + self.d * other.d
			)
		else:
			raise TypeError("Cannot multiply matrix by %s" % other)
	
	def special_inverse(self):
		""" Inverse of Special Linear group Matrices """
		return Matrix(self.d, -self.b, -self.c, self.a)
	
	# Determinant of Matrix
	@property
	def determinant(self):
		return self.a * self.d - self.b * self.c
	
	# Trace of Matrix
	@property
	def trace(self):
		return self.a + self.d
	
	
	def __str__(self):
		return "[ %s , %s | %s , %s ]" % (self.a, self.b, self.c, self.d)
	
	def __repr__(self):
		return "Matrix(%s, %s, %s, %s)" % (self.a, self.b, self.c, self.d)
	
	def to_json(self):
		return [self.a.to_json(), self.b.to_json(), self.c.to_json(), self.d.to_json()]
	
	@staticmethod
	def from_json(abcd):
		return Matrix(*map(lambda v: QuadExt.from_json(v), abcd))
	
	
	# Calculate the circle which is the intersection between the main horosphere
	# And its image under this transform
	def circle_by_height(self, height):
		from fractions import Fraction
		rad2 = Fraction(1, abs(self.c)) - height * height
		
		from circles import Circle
		return Circle(self.a.to_frac() / self.c.to_frac(), rad2)
	
	# Check if the horospherical images of the height `h` horosphere
	# For each of the transformations have any intersection
	def any_intersect(mat1, mat2, h):
		# Get radii of the two horospheres
		# Diameter (in Euclidean upper half space) of each horosphere will be 1 / (|c|^2 * h)
		r1 = Fraction(1, 2 * abs(mat1.c) * h)
		r2 = Fraction(1, 2 * abs(mat2.c) * h)
		
		# Calculate squared distance of tangency for the two radii
		# Spheres of radii `r` and `s` will be tangent for d = 2 sqrt(s * r)
		tang_dist2 = 4 * r1 * r2
		
		# Find the squared distance between the base points of each horosphere in the boundary
		# | a1 / c1 - a2 / c2 |^2 = | a1 * c2 - a2 * c1 |^2 / | c1 * c2 |^2
		d2 = Fraction(abs(mat1.a * mat2.c - mat2.a * mat1.c), abs(mat1.c * mat2.c))
		
		# If sphere are closer than tangency distance
		return d2 <= tang_dist2
	
	# Check if the horospherical images of the height `h` horosphere
	# For each of the transformations have an triple intersection
	@staticmethod
	def touching(mat1, mat2, h):
		from fractions import Fraction
		
		# Get radii of the two horospheres
		# Diameter (in Euclidean upper half space) of each horosphere will be 1 / (|c|^2 * h)
		r1 = Fraction(1, 2 * abs(mat1.c) * h)
		r2 = Fraction(1, 2 * abs(mat2.c) * h)
		
		# Find the squared distance between the base points of each horosphere in the boundary
		# | a1 / c1 - a2 / c2 |^2 = | a1 * c2 - a2 * c1 |^2 / | c1 * c2 |^2
		d2 = Fraction(abs(mat1.a * mat2.c - mat2.a * mat1.c), abs(mat1.c * mat2.c))
		
		# Check if the spheres are close enough for
		# the apex of the smaller to lie within the larger
		# This occurs if d2 <= 4 * s * (r - s) where `s` is the smaller radii and `r` the larger
		if d2 <= 4 * min(r1, r2) * abs(r1 - r2):
			# Check if minimum diameter places the apex of the inner sphere above the planar horosphere
			return 2 * min(r1, r2) >= h
		
		# Check if the highest point of the intersection of the sphere's boundaries lies above the planar horosphere
		# Height of highest and lowest intersections, for radii `r` and `s` and distance `d` will be
		#			   ___________
		#	r + s +- v/ 4sr - d^2
		#  ------------------------- = a_1, a_2
		#	   /	 (r - s)^2 \
		#	2 ( 1 + ----------- )
		#	   \		d^2	/
		# And the higher must be above or at `h` so a_1 >= h
		# So, rewriting
		#		  /	 (r - s)^2 \				 ___________
		#  k = 2h ( 1 + ----------- ) - (r + s) <= v/ 4sr - d^2
		#		  \		d^2	/
		# Thus, 4sr - d^2 >= 0 and (k <= 0 or k^2 <= 4sr - d^2)
		
		# Calculate and check the discriminant square disc2 = 4sr - d^2
		disc2 = 4 * r1 * r2 - d2
		if disc2 < 0:
			return False
		
		# Calculate the sum and difference of radii
		sm, df = r1 + r2, r1 - r2
		# Calculate `k`
		k = 2 * h * (1 + df * df / d2) - sm
		return k <= 0 or k * k <= disc2


