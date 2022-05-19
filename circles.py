import math
import matplotlib.pyplot as plt


class Parallelogram:
	# `base` is a rational complex number representing one of the vertices
	# `shift1` is the difference between `base` and one of the neighboring vertices
	# `shift2` is the difference between `base` and the other neighboring vertex
	# Vertices: [base, base + shift1, base + shift2, base + shift1 + shift2]
	# Conventions: `shift1` is treated as the horizontal axis of the parallelogram
	
	def __init__(self, base, shift1, shift2):
		self.base = base 
		self.shift1, self.shift2 = shift1, shift2
	
	def __eq__(self, other):
		return isinstance(other, Parallelogram) and self.base == other.base and self.shift1 == other.shift1 and self.shift2 == other.shift2
	
	def __hash__(self):
		return hash((self.base, self.shift1, self.shift2))
	
	
	@property
	def vertices(self):
		""" Get tuple of the four vertices that form this parallelogram """
		return (
			self.base, self.base + self.shift1,
			self.base + self.shift1 + self.shift2, self.base + self.shift2
		)
	
	# List of cells obtained by separating this parallelogram
	def cells(self, rows, columns):
		# Calculate cell shift from parallelogram shifts
		cell1, cell2 = self.shift1 / columns, self.shift2 / rows
		
		cells = []
		for r in range(rows):
			for c in range(columns):
				start = self.base + c * cell1 + r * cell2
				cells.append(Parallelogram(start, cell1, cell2))
		return cells
	
	def to_figure(self):
		return plt.Polygon(list(map(lambda vert: vert.to_point(), self.vertices)), alpha=0.3, fc='red', ec='red')


class Circle:
	def __init__(self, center, radius2):
		self.center = center
		self.radius2 = radius2
	
	def __str__(self):
		return "Circle(%s, %s)" % (str(self.center), str(self.radius2))
	
	def __repr__(self):
		return str(self)
	
	def __contains__(self, other):
		if isinstance(other, Parallelogram):
			# Check that each of the vertices lies within the circle
			return all(map(lambda v: abs(v - self.center) <= self.radius2, other.vertices))
		else:
			return abs(other - self.center) <= self.radius2
	
	def to_figure(self):
		return plt.Circle(self.center.to_point(), radius=math.sqrt(self.radius2), alpha=0.1, ec='blue')

def circles_cover(circles, parallels):
	# Check that every parallelogram is covered
	for para in parallels:
		# If none of the circles cover the cell then leave
		if not any(map(lambda circ: para in circ, circles)):
			return False
	return True

def plot_circles(circles, parallels):
	verts = [vert.to_point() for para in parallels for vert in para.vertices]
	plt.axis([-0.4, 1.4 * max(v[0] for v in verts), -0.4, 1.4 * max(v[1] for v in verts)])
#	print(circles)
#	plt.axis([-5, 5, -5, 5])
	
	for para in parallels:
		plt.gca().add_patch(para.to_figure())
	
	for circ in circles:
		plt.gca().add_patch(circ.to_figure())
	plt.show()

