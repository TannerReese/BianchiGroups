import json

from relat import Presentation

if __name__ == '__main__':
	jsonfile = input("JSON file containing Presentation: ")
	# Open json file to read presentation
	with open(jsonfile, 'r') as prsfl:
		prsnt = Presentation.from_json(json.load(prsfl))
		gens = prsnt.generators
		relats = prsnt.relations
		
		# Open magma file to write to
		magmafile = input("Magma file to write Presentaiton: ")
		with open(magmafile, 'w') as mgmfl:
			# Write line which declares generators
			mgmfl.write("F<X, Y, %s> := FreeGroup(%i);\n" % (', '.join(gens), len(gens) + 2))
			# Quotient out by the relations
			mgmfl.write("G<X, Y, %s> := quo<F | X * Y = Y * X, %s>;\n" % (', '.join(gens), ', '.join(relats)))
			# Simplify the presentation and print result
			mgmfl.write("H<[B]>, f := Simplify(G);\nH;\n")
			# Find Abelianization and print its presentation
			mgmfl.write("K<[C]>, g := AbelianQuotient(H);\nK;\n")
			# Leave
			mgmfl.write("quit;\n")

