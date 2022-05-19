import math
from quad import QuadExt, Matrix

class Presentation:
	def __init__(self, gens, denoms=None, circs=None, height=None, relats=None):
		# Dictionary of matrices and names keyed by their base in the boundary as a tuple
		self.generators = gens
		self.height = height
		
		# Collect set of denominators
		if denoms is None:
			self.denominators = set(mat.c for mat in gens.values())
		else:
			self.denominators = set(denoms)
		
		# Collect circles of each transform
		self.circles = circs
		
		self.relations = [] if relats is None else relats
	
	@staticmethod
	def from_heights(heights, fund_dom_splits=100, d=None):
		from circles import Parallelogram, circles_cover
		
		if d is not None:
			QuadExt.setd(d)
		
		# Generate parallelogram cells
		fund_dom = Parallelogram(QuadExt(0, 0).to_frac(), QuadExt(1, 0).to_frac(), QuadExt(0, 1).to_frac())
		cells = fund_dom.cells(100, 100)
		
		heights = sorted(list(heights), reverse=True)  # Make sure heights are in descending order
		for h in heights:
			print("Checking Height: ", float(h))
			# Find maximum norm for c that will result in a horosphere that intersects the plane
			max_norm = 1 / (h * h)
			print("Maximum Norm: %s <= %s" % (float(max_norm), math.ceil(max_norm)))
			
			# Find all elements with a norm less than this
			# That are non-zero
			denoms = list(filter(lambda den: abs(den) > 0, QuadExt.within_norm(max_norm)))
			
			gens = {}  # List of generators
			# `circs` contains circles for all transforms based in the fundamental domain
			#    including those on the upper boundaries
			circs = []  # List of circles corresponding to the matrices above
			genid = 0  # Identifier for generatorsg
			for c in denoms:
				print("Denominator: ", c)
				# For each `c` value find the matrices for each coprime residue
				for mat in c.matrices():
					crc = mat.circle_by_height(h)
					circs.append(crc)
					
					# Check if matrix is a generator or not
					if crc.center.x < 1 and crc.center.y < 1:
						print(mat)
						gens['A' + str(genid)] = mat
						genid += 1
			
			# Check if circles cover fundamental domain
			if circles_cover(circs, cells):
				return Presentation(gens, denoms=denoms, circs=circs, height=h)
		
		raise ValueError("Given heights do not admit a presentation")
	
	def sphere_shifts(self, denom):
		"""
		For a given denominator and height, this returns the shifts of the fundamental domain
		which when applied to the generating set of matrices will produce a set of transformations
		containing every horosphere that may triple intersect an element of the generating set.
		"""
		from fractions import Fraction
		
		# Calculate the radius for denom
		# Diameter (in Euclidean upper half space) of each horosphere will be 1 / (|c|^2 * h)
		r = Fraction(1, 2 * abs(denom) * self.height)
		# Calculate the maximum possible radius corresponding to |c|^2 = 1
		rm = Fraction(1, 2 * self.height)
		
		# Calculate squared distance of tangency for the radii
		# Spheres of radii `r` and `s` will be tangent for d = 2 sqrt(s * r)
		tang_dist2 = 4 * r * rm
		
		# Magnitude of the shifts must be less than or equal to
		# the sum of the fundamental domain and the tangency distance
		# Longer diagonal of fundamental domain from 0 to 1 + w represents diameter
		diam2 = abs(QuadExt(1, 1))
		# Since shift_mag <= tang_dist + diam
		# shift_mag^2 <= tang_dist^2 + 2 * diam * tand_dist + diam^2
		shift_mag2 = tang_dist2 + math.ceil(2 * math.sqrt(tang_dist2 * diam2)) + diam2
		
		# Find shifts that have a squared magnitude / norm close enough to potential have an effect
		return QuadExt.within_norm(shift_mag2 * 10)
	
	def relations_by_denom(self, denom):
		"""
		For a given denominator find all of the relations between
		generators that lie in the fundamental domain and other transforms
		"""
		
		# Get all valid shifts and convert them to transforms
		shifts = self.sphere_shifts(denom)
		
		# Iterate over all possible shifts
		for nmg, g in self.generators.items():
			# Only consider generators with given denominator
			if g.c == denom:
				# Get inverse relation for generator
				[nm1, shf2, nm3, shf3] = self.get_inverse_relation(g, nmg)
				print("Found relation for %s^-1" % nmg)
				
				# Make relation into string
				# Write as single line equal to the identity
				rel = "%s * X^%i * Y^%i * %s * X^%i * Y^%i" % (nm1, shf2.x, shf2.y, nm3, shf3.x, shf3.y)
				self.relations.append(rel)
				
				# Iterate over all generators
				for nmb, b in self.generators.items():
					# Ignore when generators are the same
					if nmb == nmg:
						continue
					
					for s in shifts:
						if Matrix.touching(g, Matrix.shift(s) * b, self.height):						
							
							try:
								# Get relation for two horospheres
								[nm1, shf1, nm2, shf2, nm3, shf3] = self.get_relation(g, nmg, b, nmb, s)
								print("Found relation between ", nmg, "&", nmb)
								
								# Make relation into string
								# Write as single line equal to the identity
								rel = "%s^-1 * X^%i * Y^%i * %s * X^%i * Y^%i * %s * X^%i * Y^%i" % (nm2, -shf1.x, -shf1.y, nm1, shf2.x, shf2.y, nm3, shf3.x, shf3.y)
								self.relations.append(rel)
							except ZeroDivisionError:
								# When 
								pass
	
	def find_relations(self):
		self.relations = []
		# Get set of denominators of generators
		for c in self.denominators:
			self.relations_by_denom(c)
	
	
	def get_relation(self, mat1, name1, mat2, name2, shift):
		"""
		Finds the relation corresponding to the triple intersection of
		the horoball for `mat1` and the horoball for `shift * mat2`
		
		Each relation is returned as a list,
		[name1, shift, name2, shift2, name3, shift3]
		
		Corresponding to the relation
		mat1^-1 * shift * mat2 = shift2 * mat3 * shift3
		
		Args:
			mat1 : Matrix  --  Matrix which is one of the generators
			name1 : str  --  Name used to identify generator `mat1`
			mat2 : Matrix  --  Matrix which is one of the generators
			name2 : str  --  Name used to identify generator `mat2`
			shift : QuadExt  --  Shift to apply to mat2
		"""
		
		# Calculate product matrix
		prod = mat1.special_inverse() * Matrix.shift(shift) * mat2
		
		# Separate out post-shift (shift2)
		# Find residue of prod's basepoint lying in the fundamental domain
		shift2 = prod.a // prod.c
		# Get products residue in fundamental domain
		shifted_gener = Matrix.shift(-shift2) * prod
		
		# shifted_gener has the same horoball as one of the generators
		# Find generator with the corresponding horoball
		name3, gener = None, None
		basepoint = shifted_gener.a.to_frac() / shifted_gener.c.to_frac()
		for nm, mat in self.generators.items():
			if mat.a.to_frac() / mat.c.to_frac() == basepoint:
				name3, gener = nm, mat
				break
		
		# Get pre-shift (shift3)
		shift3 = (gener.special_inverse() * shifted_gener).to_shift()
		return [name1, shift, name2, shift2, name3, shift3]
	
	def get_inverse_relation(self, mat1, name1):
		"""
		Finds the relation corresponding to the triple intersection of
		the horoball for `mat1` and the main horoball twice
		
		Each relation is returned as a list,
		[name1, shift2, name3, shift3]
		
		Corresponding to the relation
		mat1^-1 = shift2 * mat3 * shift3
		
		Args:
			mat1 : Matrix  --  Matrix which is one of the generators
			name1 : str  --  Name used to identify generator `mat1`
		"""
		
		inv = mat1.special_inverse()
		# Separate out post-shift (shift2)
		shift2 = inv.a // inv.c
		# Get inv's residue in fundamental domain
		shifted_gener = Matrix.shift(-shift2) * inv
		
		# shifted_gener has the some horoball as one of the generators
		# Find generator with the corresponding horoball
		name3, gener = None, None
		basepoint = shifted_gener.a.to_frac() / shifted_gener.c.to_frac()
		for nm, mat in self.generators.items():
			if mat.a.to_frac() / mat.c.to_frac() == basepoint:
				name3, gener = nm, mat
				break
		
		# Get pre-shift (shift3)
		shift3 = (gener.special_inverse() * shifted_gener).to_shift()
		return [name1, shift2, name3, shift3]
	
	
	
	# Plot presentation at given height
	def plot(self, fund_dom_splits=100):
		"""
		Plot fundamental domain with the circles for each generating set for a given height
		"""
		from circles import Parallelogram, plot_circles
		
		# Calculate the parallelogram cells
		fund_dom = Parallelogram(QuadExt(0, 0).to_frac(), QuadExt(1, 0).to_frac(), QuadExt(0, 1).to_frac())
		cells = fund_dom.cells(fund_dom_splits, fund_dom_splits)
		
		plot_circles(self.circles, cells)
	
	
	def to_json(self):
		return {
			'height': str(self.height),
			'dvalue': QuadExt.DVALUE,
			'generators': {name: mat.to_json() for name, mat in self.generators.items()},
			'circles': [{'cnt': [str(crc.center.x), str(crc.center.y)], 'rad2': str(crc.radius2)}for crc in self.circles],
			'relations': self.relations
		}
	
	@staticmethod
	def from_json(obj):
		from fractions import Fraction
		from circles import Circle
		gens = {name: Matrix.from_json(mat) for name, mat in obj['generators'].items()}
		circs = {Circle(QuadExt(Fraction(crc['cnt'][0]), Fraction(crc['cnt'][1])), Fraction(crc['rad2'])) for crc in obj['circles']}
		
		# Set discriminant
		QuadExt.set_disc(obj['discriminant'])
		return Presentation(gens, circs=circs, height=Fraction(obj['height']), relats=obj['relations'])

