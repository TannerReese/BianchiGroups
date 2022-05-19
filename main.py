from circles import Parallelogram, plot_circles, circles_cover
from fractions import Fraction
from relat import Presentation

def matrices_for_height(height):
	print("\n\nChecking height %s" % height)
	
	# Find maximum norm for c that will result in a horosphere that intersects the plane
	max_norm = 1 / (height * height)
	print("Maximum Norm %s" % max_norm)
	
	# Find all elements with a norm less than this
	denoms = QuadExt.within_norm(max_norm)
	
	matrices = []
	circs = []
	for c in denoms:
		print("Denominator %s" % c)
		
		# For each `c` value find the matrices for each coprime residue
		cmats = c.matrices()
		print("Matrices {\n%s\n}" % '\n'.join(map(str, cmats)))
		matrices += cmats
	return matrices



def plot_height(height, fund_dom_splits=100):
	"""
	Plot fundamental domain with the circles for each generating set for a given height
	"""
	
	# Get the circle for each matrix
	mats = matrices_for_height(height)
	print(mats)
	circs = list(mt.circle_by_height(height) for mt in mats)
		
	# Calculate the parallelogram cells
	fund_dom = Parallelogram(QuadExt(0, 0).to_frac(), QuadExt(1, 0).to_frac(), QuadExt(0, 1).to_frac())
	cells = fund_dom.cells(fund_dom_splits, fund_dom_splits)
	
	plot_circles(circs, cells)


def check_heights(heights, fund_dom_splits=100):
	# Calculate the parallelogram cells
	fund_dom = Parallelogram(QuadExt(0, 0).to_frac(), QuadExt(1, 0).to_frac(), QuadExt(0, 1).to_frac())
	cells = fund_dom.cells(fund_dom_splits, fund_dom_splits)
	
	for h in heights:
		print("\n\nChecking height %s" % h)
		
		# Find maximum norm for c that will result in a horosphere that intersects the plane
		max_norm = 1 / (h * h)
		print("Maximum Norm %s" % max_norm)
		
		# Find all elements with a norm less than this
		denoms = QuadExt.within_norm(max_norm)
		
		matrices = []
		circs = []
		for c in denoms:
			print("Denominator %s" % c)
			
			# For each `c` value find the matrices for each coprime residue
			cmats = c.matrices()
			print("Matrices {\n%s\n}" % '\n'.join(map(str, cmats)))
			matrices += cmats
			
			# Get the circle for each matrix
			circs += [mat.circle_by_height(h) for mat in cmats]
		
		# Check if circles cover the fundamental domain
		if circles_cover(circs, cells):
			return (h, matrices)
		else:
			print("Did not cover")
	
	return (None, None)


if __name__ == '__main__':
	import os.path
	import json
	
	filename = input("Read From: ")
	if os.path.isfile(filename):
		# Read generators from file as json
		with open(filename, 'r') as fl:
			print("Reading file", filename)
			prsnt = Presentation.from_json(json.load(fl))
	else:
		print("No file found at '%s'" % filename)
		
		# Select particular ring of integers
		dval = -abs(int(input("d value (2, 7, 11, 19, 43, 67, 163)? ")))
		print("Generating Presentation")
		prsnt = Presentation.from_heights([Fraction(29, 30) ** i / 2 for i in range(50)], fund_dom_splits=500, d=dval)
	
	# Plot graph
	prsnt.plot()
	
	# Calculate Relations
	if input("Calculate Relations? (y/N)").lower().strip() in ['y', 'yes']:
		filename = ""
		while len(filename) == 0:
			filename = input("Save To: ")
		
		with open(filename, 'w') as fl:
			# Calculate Relations
			print("Calculating Relations")
			prsnt.find_relations()
			print(prsnt.to_json())
			print("Total Relations: ", len(prsnt.relations))
			
			# Write to file
			print("Saving to", filename)
			json.dump(prsnt.to_json(), fl)

