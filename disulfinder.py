##### Disulfinder #####
# import __main__
# __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import sys
import pymol
from math import cos
from pymol import cmd, stored
import time
import networkx as nx

def compute_energy(chi1_i, chi1_j, theta_i, theta_j, chi3):
	# resnum_i_pymol_prefix = "/" + pdb + "///" + resnum_i + "/"
	# resnum_j_pymol_prefix = "/" + pdb + "///" + resnum_j + "/"

	# chi1_i = cmd.get_dihedral(resnum_i_pymol_prefix + "N", resnum_i_pymol_prefix + "CA", resnum_i_pymol_prefix + "CB", resnum_i_pymol_prefix + "SG")
	# theta_i = cmd.get_angle(resnum_i_pymol_prefix + "CA", resnum_i_pymol_prefix + "CB", resnum_i_pymol_prefix + "SG")

	# chi1_j = cmd.get_dihedral(resnum_j_pymol_prefix + "N", resnum_j_pymol_prefix + "CA", resnum_j_pymol_prefix + "CB", resnum_j_pymol_prefix + "SG")
	# theta_j = cmd.get_angle(resnum_j_pymol_prefix + "CA", resnum_j_pymol_prefix + "CB", resnum_j_pymol_prefix + "SG")

	# chi3 = cmd.get_dihedral(resnum_i_pymol_prefix + "CB", resnum_i_pymol_prefix + "SG", resnum_j_pymol_prefix + "SG", resnum_j_pymol_prefix + "CB")

	E_chi1_i = 1.4 * (1 + cos (3 * chi1_i))
	E_theta_i = 55.0 * (theta_i - 114.6)**2
	E_chi1_j = 1.4 * (1 + cos (3 * chi1_j))
	E_theta_j = 55.0 * (theta_j - 114.6)**2
	E_chi3 = 4.0 * (1 - cos((2 * chi3) + 160))

	Energy = E_chi1_i + E_theta_i + E_chi1_j + E_theta_j + E_chi3

	return Energy

# def mutate_to_cys(pdb, resnum_i):
# 	selection_i = "resi " + str(resnum_i)
# 	cmd.wizard("mutagenesis")

# 	cmd.get_wizard().set_mode("CYS")

# 	cmd.do("refresh_wizard")

# 	cmd.get_wizard().do_select(selection_i)
# 	cmd.get_wizard().do_select(selection_i)

# 	for frame_i in [1, 2, 3]:
# 		cmd.frame(frame_i)
# 		cmd.create("rotamer_" + str(resnum_i) + "_" + str(frame_i), "mutation", frame_i, 1)
# 		rotamer_i_pymol_prefix = "/" + "rotamer_" + str(resnum_i) + "_" + str(frame_i) + "///" + resnum_i + "/"

# 	cmd.set_wizard("done")

def mutate_to_cys(pdb, resnum_i):
	selection_i = "resi " + str(resnum_i)
	cmd.wizard("mutagenesis")

	cmd.get_wizard().set_mode("CYS")

	cmd.do("refresh_wizard")

	cmd.get_wizard().do_select(selection_i)
	cmd.get_wizard().do_select(selection_i)

	for frame_i in [1, 2, 3]:
		cmd.frame(frame_i)
		cmd.create("rotamer_" + str(resnum_i) + "_" + str(frame_i), "mutation", frame_i, 1)
		rotamer_i_pymol_prefix = "/" + "rotamer_" + str(resnum_i) + "_" + str(frame_i) + "///" + resnum_i + "/"

	cmd.set_wizard("done")

def check_disulphide_criteria(pdb, resnum_i, resnum_j, f):
	selection_i = "resi " + str(resnum_i)
	selection_j = "resi " + str(resnum_j)

	for frame_i in [1, 2]:
		rotamer_i_pymol_prefix = "/" + "rotamer_" + str(resnum_i) + "_" + str(frame_i) + "///" + resnum_i + "/"

		chi1_i = cmd.get_dihedral(rotamer_i_pymol_prefix + "N", rotamer_i_pymol_prefix + "CA", rotamer_i_pymol_prefix + "CB", rotamer_i_pymol_prefix + "SG")
		theta_i = cmd.get_angle(rotamer_i_pymol_prefix + "CA", rotamer_i_pymol_prefix + "CB", rotamer_i_pymol_prefix + "SG")

		for frame_j in [1, 2]:
			rotamer_j_pymol_prefix = "/" + "rotamer_" + str(resnum_j) + "_" + str(frame_j) + "///" + resnum_j + "/"
	
			chi1_j = cmd.get_dihedral(rotamer_j_pymol_prefix + "N", rotamer_j_pymol_prefix + "CA", rotamer_j_pymol_prefix + "CB", rotamer_j_pymol_prefix + "SG")
			theta_j = cmd.get_angle(rotamer_j_pymol_prefix + "CA", rotamer_j_pymol_prefix + "CB", rotamer_j_pymol_prefix + "SG")

			chi3 = cmd.get_dihedral(rotamer_i_pymol_prefix + "CB", rotamer_i_pymol_prefix + "SG", rotamer_j_pymol_prefix + "SG", rotamer_j_pymol_prefix + "CB")

			#Energy = compute_energy(chi1_i, chi1_j, theta_i, theta_j, chi3)

			s_gamma_distance = cmd.get_distance(rotamer_i_pymol_prefix + "SG", rotamer_j_pymol_prefix + "SG")	
			if (s_gamma_distance < 3 and ((chi3 > -110 and chi3 < -60) or (chi3 > 70 and chi3 < 130))):
				print rotamer_i_pymol_prefix, rotamer_j_pymol_prefix, s_gamma_distance, chi3
				f.write("%s\t%s\t%f\t%f\n" % (rotamer_i_pymol_prefix, rotamer_j_pymol_prefix, s_gamma_distance, chi3))
				rotamer_i = "rotamer_" + str(resnum_i) + "_" + str(frame_i)
				rotamer_j = "rotamer_" + str(resnum_j) + "_" + str(frame_j)
				cmd.enable(rotamer_i)
				cmd.enable(rotamer_j)
				cmd.distance(rotamer_i_pymol_prefix + "SG", rotamer_j_pymol_prefix + "SG")
				cmd.create(rotamer_i + rotamer_j, rotamer_i_pymol_prefix + "or" + rotamer_j_pymol_prefix)
				# cmd.bond(rotamer_i_pymol_prefix + "SG", rotamer_j_pymol_prefix + "SG")


def check_disulphide(pdb, resnum_i, resnum_j, f):

	selection_i = "resi " + str(resnum_i)
	selection_j = "resi " + str(resnum_j)

	# cmd.wizard("mutagenesis")

	# cmd.get_wizard().set_mode("CYS")

	cmd.do("refresh_wizard")
	# Select the rotamer
	for frame_i in [1, 2, 3]:
		# cmd.do("refresh_wizard")

		cmd.get_wizard().do_select(selection_i)
		cmd.frame(frame_i)
		cmd.create("rotamer_i_" + str(frame_i), "mutation", frame_i, 1)
		rotamer_i_pymol_prefix = "/" + "rotamer_i_" + str(frame_i) + "///" + resnum_i + "/"


	for frame_j in [1, 2, 3]:

		# cmd.do("refresh_wizard")

		cmd.get_wizard().do_select(selection_j)
		cmd.frame(frame_j)
		cmd.create("rotamer_j_" + str(frame_j), "mutation", frame_j, 1)
		rotamer_j_pymol_prefix = "/" + "rotamer_j_" + str(frame_j) + "///" + resnum_j + "/"

		# s_gamma_distance = cmd.get_distance(rotamer_i_pymol_prefix + "SG", rotamer_j_pymol_prefix + "SG")
		# print s_gamma_distance

	cmd.set_wizard("clear")

	for frame_i in [1, 2, 3]:
		rotamer_i_pymol_prefix = "/" + "rotamer_i_" + str(frame_i) + "///" + resnum_i + "/"
		for frame_j in [1, 2, 3]:
			rotamer_j_pymol_prefix = "/" + "rotamer_j_" + str(frame_j) + "///" + resnum_j + "/"

			s_gamma_distance = cmd.get_distance(rotamer_i_pymol_prefix + "SG", rotamer_j_pymol_prefix + "SG")

			print s_gamma_distance

"""For a given resdiue, get all neighbours within a given cutoff distance"""
def get_neighbours(resnum_i, cutoff_dist, pdb):
	stored.neighbours_list = []
	pymol_neighbour_selection = "name CA and (byres resi " + resnum_i + " around 4.5) and " + pdb
	cmd.iterate((pymol_neighbour_selection), "stored.neighbours_list.append(resi)")
	return stored.neighbours_list


def check_disulphide_in_pdb(pdb):

	"""load structure from file"""
	structure_file = "/Users/ajvenkatakrishnan/Downloads/" + pdb + ".pdb"
	cmd.load(structure_file, pdb)

	"""output file for the list of disulphide bridge candidates"""
	disulphides_outfile = "/Users/ajvenkatakrishnan/Downloads/" + pdb + "_disulphide.txt"
	f = open(disulphides_outfile, 'w')

	""" Create a selection of residues for disulphide scanning """
	stored.list = []
	pymol_selection = "name CA and (resi 127-167 or resi 198-245 or resi 280-360 or resi 382-412) and " + pdb
	# pymol_selection = "name CA and (resi 130-148 or resi 406-408)" + pdb
	#pymol_selection = "name CA and (resi 42-88 or resi 111-166 or resi 209-288 or resi 318-342) and " + pdb
	cmd.iterate((pymol_selection), "stored.list.append(resi)")
	selection_list = stored.list

	""" Generate Cys rotamers at all the positions in the selection """
	for resnum_i in stored.list:
		mutate_to_cys(pdb, resnum_i)

	cmd.disable("rotamer*")

	""" Shortlist Cys-pairs that satisfy the criteria for disulphide bond formation """
	disulphide_graph=nx.empty_graph()
	cutoff_dist = 4.5

	for resnum_i in selection_list:
		neighbours_list = get_neighbours(resnum_i, cutoff_dist, pdb)
		for resnum_j in neighbours_list:	
			if (int (resnum_i) - int (resnum_j) > 4 or int (resnum_i) - int (resnum_j) < -4):
				if str(resnum_j) in selection_list:
					if not disulphide_graph.has_edge(str(resnum_i), str(resnum_j)):
						check_disulphide_criteria(pdb, str(resnum_i), str(resnum_j), f)
						disulphide_graph.add_edge(str(resnum_i), str(resnum_j))

pymol.finish_launching()
start_time = time.time()

### Load protein
# pdb = "2RH1"


#check_disulphide_in_pdb("3P0G")
#cmd.reinitalize()

check_disulphide_in_pdb("2RH1")

#pymol_selection_residues = "name CA and (resi 130-148 or resi 406-408) and "

cmd.do("hide everything; show cartoon; show sticks, resn CYS and rota*; zoom polymer")

# pdb = "5u09"


# structure_file = "/Users/ajvenkatakrishnan/Downloads/" + pdb + ".pdb"
# print cmd.load(structure_file, pdb)

# disulphides_outfile = "/Users/ajvenkatakrishnan/Downloads/"+ pdb + "_disulphide.txt"
# f = open(disulphides_outfile, 'w')

# """ Create a selection of residues for disulphide scanning"""
# stored.list = []
# # pymol_selection = "name CA and (resi 127-167 or resi 198-245 or resi 280-360 or resi 382-412) and " + pdb
# pymol_selection = "name CA and (resi 130-148 or resi 406-408) and " + pdb
# #pymol_selection = "name CA and (resi 42-88 or resi 111-166 or resi 209-288 or resi 318-342) and " + pdb
# cmd.iterate((pymol_selection), "stored.list.append(resi)")
# selection_list = stored.list

# """ Generate Cys rotamers at all the positions in the selection """
# for resnum_i in stored.list:
# 	mutate_to_cys(pdb, resnum_i)

# cmd.disable("rotamer*")

# """ Shortlist Cys-pairs that best satisfy the criteria for disulphide bond formation"""
# disulphide_graph=nx.empty_graph()
# for resnum_i in selection_list:
# 	cutoff_dist = 4.5
# 	neighbours_list = get_neighbours(resnum_i, cutoff_dist, pdb)
# 	for resnum_j in neighbours_list:	
# 		if (int (resnum_i) - int (resnum_j) > 4 or int (resnum_i) - int (resnum_j) < -4):
# 			if str(resnum_j) in selection_list:
# 				if not disulphide_graph.has_edge(str(resnum_i), str(resnum_j)):
# 				#print "pair", resnum_i, resnum_j
# 					check_disulphide_criteria(pdb, str(resnum_i), str(resnum_j), f)
# 					disulphide_graph.add_edge(str(resnum_i), str(resnum_j))

# cmd.do("hide everything; show cartoon; show sticks, resn CYS and rota*; zoom polymer")