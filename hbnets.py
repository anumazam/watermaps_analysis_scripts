## Determines the hydrogen bonding network in the ligand binding pocket(s) of an n-meric protein

'''
Args:
	file = PDB file name for which to determine hbonds
	directory_files = path to the file
	
Returns:
	Visualization of H bond network in PyMol
	Dictionary with identified H bond participants and H bond lengths between heavy atoms
'''

import os.path
import os
import string

from pymol import cmd
from pymol import stored

def hbnets(file, directory_files):

	letters = string.ascii_uppercase

	file_num = file[-10] + file [-9]
	
	#water
	W = "W_{}".format(file_num)
	cmd.select(W,"solvent")
	#stored.water = []
	#cmd.iterate("W","stored.water.append(chain + resn + resi + name)")
	#ligand
	L = "L_{}".format(file_num)
	cmd.select(L, "hetatm and not " + W)
	#stored.ligand = []
	#cmd.iterate("L","stored.ligand.append(chain + resn + resi + name)")
	#protein
	P = "P_{}".format(file_num)
	cmd.select(P,"not hetatm and not " + W + " and not " + L + " and L around 6")
	#stored.protein = []
	#cmd.iterate("P","stored.protein.append(chain + resn + resi + name)")

	#donors
	cmd.select("don", ("(elem n or elem o) and (neighbor hydro)"))
	#stored.donors = []
	#cmd.iterate("don","stored.donors.append(chain + resn + resi + name)")
	#acceptors
	cmd.select("acc", ("elem o or (elem n and not (neighbor hydro))"))
	#stored.acceptors = []
	#cmd.iterate("acc","stored.acceptors.append(chain + resn + resi + name)")

	#create lists of all atom categories
	# w_acc = [w for w in stored.water if w in stored.acceptors]
	# w_don = [w for w in stored.water if w in stored.donors]
	# l_acc = [l for l in stored.ligand if l in stored.acceptors]
	# l_don = [l for l in stored.ligand if l in stored.donors]
	# p_acc = [p for p in stored.protein if p in stored.acceptors]
	# p_don = [p for p in stored.protein if p in stored.donors]

	#create selections of all atom categories
	cmd.select("l_don", L + " and don")
	cmd.select("w_acc", W + " and acc")
	cmd.select("w_don", W + " and don")
	cmd.select("l_acc", L + " and acc")
	cmd.select("p_don", P + " and don")
	cmd.select("p_acc", P + " and acc")

	stored.l_don_resi = []
	stored.l_don_name = []
	cmd.iterate("l_don","stored.l_don_resi.append(resi)")
	cmd.iterate("l_don","stored.l_don_name.append(name)")

	stored.w_don_resi = []
	stored.w_don_name = []
	cmd.iterate("w_don","stored.l_don_resi.append(resi)")
	cmd.iterate("w_don","stored.l_don_name.append(name)")

	stored.p_don_resi = []
	stored.p_don_name = []
	cmd.iterate("p_don","stored.l_don_resi.append(resi)")
	cmd.iterate("p_don","stored.l_don_name.append(name)")

	stored.w_acc_resi = []
	stored.w_acc_name = []
	cmd.iterate("w_acc","stored.w_acc_resi.append(resi)")
	cmd.iterate("w_acc","stored.w_acc_name.append(name)")

	stored.l_acc_resi = []
	stored.l_acc_name = []
	cmd.iterate("l_acc","stored.w_acc_resi.append(resi)")
	cmd.iterate("l_acc","stored.w_acc_name.append(name)")

	stored.p_acc_resi = []
	stored.p_acc_name = []
	cmd.iterate("p_acc","stored.w_acc_resi.append(resi)")
	cmd.iterate("p_acc","stored.w_acc_name.append(name)")

	hbonds = {}

	os.chdir(directory_files)
	lines = []
	with open(file) as f:
		for line in f.readlines():
			l = line.split()
			lines.append(l)

	#L---W
	for i,atm1 in enumerate(stored.l_don_resi):
		atm1_ = stored.l_don_name[i]
		for i,atm2 in enumerate(stored.w_acc_resi):
			atm2_ = stored.w_acc_name[i]
			HB_LW = "HB_LW_{}".format(file_num)
			dst = cmd.distance(HB_LW, atm1 + "/" + atm1_, atm2 + "/" + atm2_, 3.2)
			if dst < 3.2 and dst > 0:
				for i in lines:
					try:
						if i[5] == atm1 and i[2] == atm1_:
							if i[4] in letters:
								res1 = i[2] + i[3] + i[4] + i[5]
							else:
								res1 = i[2] + i[3] + i[4]
							break
						if i[4] == atm1 and i[1] == atm1_:
							if i[3] in letters:
								res1 = i[1] + i[2] + i[3] + i[4]
							else:
								res1 = i[1] + i[2] + i[3]
							break
					except:
						continue
				for i in lines:
					try:
						if i[5] == atm2 and i[2] == atm2_:
							if i[4] in letters:
								res2 = i[2] + i[3] + i[4] + i[5]
							else:
								res2 = i[2] + i[3] + i[4]
							break
						if i[4] == atm2 and i[1] == atm2_:
							if i[3] in letters:
								res2 = i[1] + i[2] + i[3] + i[4]
							else:
								res2 = i[1] + i[2] + i[3]
							break
					except:
						continue
				hbonds[(res1, res2)] = round(dst, 2)

	#W---L
	for i,atm1 in enumerate(stored.w_don_resi):
		atm1_ = stored.w_don_name[i]
		for i,atm2 in enumerate(stored.l_acc_resi):
			atm2_ = stored.l_acc_name[i]
			HB_WL = "HB_WL_{}".format(file_num)
			dst = cmd.distance(HB_LW, atm1 + "/" + atm1_, atm2 + "/" + atm2_, 3.2)
			if dst < 3.2 and dst > 0:
				for i in lines:
					try:
						if i[5] == atm1 and i[2] == atm1_:
							if i[4] in letters:
								res1 = i[2] + i[3] + i[4] + i[5]
							else:
								res1 = i[2] + i[3] + i[4]
							break
						if i[4] == atm1 and i[1] == atm1_:
							if i[3] in letters:
								res1 = i[1] + i[2] + i[3] + i[4]
							else:
								res1 = i[1] + i[2] + i[3]
							break
					except:
						continue
				for i in lines:
					try:
						if i[5] == atm2 and i[2] == atm2_:
							if i[4] in letters:
								res2 = i[2] + i[3] + i[4] + i[5]
							else:
								res2 = i[2] + i[3] + i[4]
							break
						if i[4] == atm2 and i[1] == atm2_:
							if i[3] in letters:
								res2 = i[1] + i[2] + i[3] + i[4]
							else:
								res2 = i[1] + i[2] + i[3]
							break
					except:
						continue
				hbonds[(res1, res2)] = round(dst, 2)
	
	#P---L
	# for i,atm1 in enumerate(stored.p_don_resi):
	# 	atm1_ = stored.p_don_name[i]
	# 	for i,atm2 in enumerate(stored.l_acc_resi):
	# 		atm2_ = stored.l_acc_name[i]
	# 		HB_PL = "HB_PL_{}".format(file_num)
	# 		dst = cmd.distance(HB_LW, atm1 + "/" + atm1_, atm2 + "/" + atm2_, 3.2)
	# 		if dst < 3.2 and dst > 0:
	# 			for i in lines:
	# 				try:
	# 					if i[5] == atm1 and i[2] == atm1_:
	# 						if i[4] in letters:
	# 							res1 = i[2] + i[3] + i[4] + i[5]
	# 						else:
	# 							res1 = i[2] + i[3] + i[4]
	# 						break
	# 					if i[4] == atm1 and i[1] == atm1_:
	# 						if i[3] in letters:
	# 							res1 = i[1] + i[2] + i[3] + i[4]
	# 						else:
	# 							res1 = i[1] + i[2] + i[3]
	# 						break
	# 				except:
	# 					continue
	# 			for i in lines:
	# 				try:
	# 					if i[5] == atm2 and i[2] == atm2_:
	# 						if i[4] in letters:
	# 							res2 = i[2] + i[3] + i[4] + i[5]
	# 						else:
	# 							res2 = i[2] + i[3] + i[4]
	# 						break
	# 					if i[4] == atm2 and i[1] == atm2_:
	# 						if i[3] in letters:
	# 							res2 = i[1] + i[2] + i[3] + i[4]
	# 						else:
	# 							res2 = i[1] + i[2] + i[3]
	# 						break
	# 				except:
	# 					continue
	# 			hbonds[(res1, res2)] = round(dst, 2)
	
	#L---P
	# for i,atm1 in enumerate(stored.l_don_resi):
	# 	atm1_ = stored.l_don_name[i]
	# 	for i,atm2 in enumerate(stored.p_acc_resi):
	# 		atm2_ = stored.p_acc_name[i]
	# 		HB_LP = "HB_LP_{}".format(file_num)
	# 		dst = cmd.distance(HB_LW, atm1 + "/" + atm1_, atm2 + "/" + atm2_, 3.2)
	# 		if dst < 3.2 and dst > 0:
	# 			for i in lines:
	# 				try:
	# 					if i[5] == atm1 and i[2] == atm1_:
	# 						if i[4] in letters:
	# 							res1 = i[2] + i[3] + i[4] + i[5]
	# 						else:
	# 							res1 = i[2] + i[3] + i[4]
	# 						break
	# 					if i[4] == atm1 and i[1] == atm1_:
	# 						if i[3] in letters:
	# 							res1 = i[1] + i[2] + i[3] + i[4]
	# 						else:
	# 							res1 = i[1] + i[2] + i[3]
	# 						break
	# 				except:
	# 					continue
	# 			for i in lines:
	# 				try:
	# 					if i[5] == atm2 and i[2] == atm2_:
	# 						if i[4] in letters:
	# 							res2 = i[2] + i[3] + i[4] + i[5]
	# 						else:
	# 							res2 = i[2] + i[3] + i[4]
	# 						break
	# 					if i[4] == atm2 and i[1] == atm2_:
	# 						if i[3] in letters:
	# 							res2 = i[1] + i[2] + i[3] + i[4]
	# 						else:
	# 							res2 = i[1] + i[2] + i[3]
	# 						break
	# 				except:
	# 					continue
	# 			hbonds[(res1, res2)] = round(dst, 2)
	
	#P---W
	for i,atm1 in enumerate(stored.p_don_resi):
		atm1_ = stored.p_don_name[i]
		for i,atm2 in enumerate(stored.w_acc_resi):
			atm2_ = stored.w_acc_name[i]
			HB_PW = "HB_PW_{}".format(file_num)
			dst = cmd.distance(HB_LW, atm1 + "/" + atm1_, atm2 + "/" + atm2_, 3.2)
			if dst < 3.2 and dst > 0:
				for i in lines:
					try:
						if i[5] == atm1 and i[2] == atm1_:
							if i[4] in letters:
								res1 = i[2] + i[3] + i[4] + i[5]
							else:
								res1 = i[2] + i[3] + i[4]
							break
						if i[4] == atm1 and i[1] == atm1_:
							if i[3] in letters:
								res1 = i[1] + i[2] + i[3] + i[4]
							else:
								res1 = i[1] + i[2] + i[3]
							break
					except:
						continue
				for i in lines:
					try:
						if i[5] == atm2 and i[2] == atm2_:
							if i[4] in letters:
								res2 = i[2] + i[3] + i[4] + i[5]
							else:
								res2 = i[2] + i[3] + i[4]
							break
						if i[4] == atm2 and i[1] == atm2_:
							if i[3] in letters:
								res2 = i[1] + i[2] + i[3] + i[4]
							else:
								res2 = i[1] + i[2] + i[3]
							break
					except:
						continue
				hbonds[(res1, res2)] = round(dst, 2)

	#W---P
	for i,atm1 in enumerate(stored.w_don_resi):
		atm1_ = stored.w_don_name[i]
		for i,atm2 in enumerate(stored.p_acc_resi):
			atm2_ = stored.p_acc_name[i]
			HB_WP = "HB_WP_{}".format(file_num)
			dst = cmd.distance(HB_LW, atm1 + "/" + atm1_, atm2 + "/" + atm2_, 3.2)
			if dst < 3.2 and dst > 0:
				for i in lines:
					try:
						if i[5] == atm1 and i[2] == atm1_:
							if i[4] in letters:
								res1 = i[2] + i[3] + i[4] + i[5]
							else:
								res1 = i[2] + i[3] + i[4]
							break
						if i[4] == atm1 and i[1] == atm1_:
							if i[3] in letters:
								res1 = i[1] + i[2] + i[3] + i[4]
							else:
								res1 = i[1] + i[2] + i[3]
							break
					except:
						continue
				for i in lines:
					try:
						if i[5] == atm2 and i[2] == atm2_:
							if i[4] in letters:
								res2 = i[2] + i[3] + i[4] + i[5]
							else:
								res2 = i[2] + i[3] + i[4]
							break
						if i[4] == atm2 and i[1] == atm2_:
							if i[3] in letters:
								res2 = i[1] + i[2] + i[3] + i[4]
							else:
								res2 = i[1] + i[2] + i[3]
							break
					except:
						continue
				hbonds[(res1, res2)] = round(dst, 2)

	print(hbonds)	

	#cmd.iterate("HB_LW","stored.hbonds.append([chain + resn + resi + name])")

	#naming scheme: HB_don acc
	# for num, atm1 in enumerate(w_acc):
	# 	for atm2 in l_don:
	# 		cmd.distance("HB_LW", "W"[num], atm2, 3.2)

	# cmd.distance("HB_LP", ("P and acc"), ("L and don"), 3.2)
	# cmd.distance("HB_PL", ("P and don"), ("L and acc"), 3.2)

	# cmd.distance("HB_WP", ("P and acc"), ("W and don"), 3.2)
	# cmd.distance("HB_PW", ("P and don"), ("W and acc"), 3.2)

	
	# cmd.distance("HB_WL", ("L and acc"), ("W and don"), 3.2)

	# #show sticks in area of interest
	cmd.show_as("sticks", "byres " + L + " around 6")
	cmd.hide("everything", "hydrogen")
	cmd.show_as("sticks", L)
	cmd.show_as("sticks", W)

	# stored.hbonds = []
	# cmd.iterate("all", "stored.hbonds.append(LP)")
	# print(stored.hbonds)

	#delete selections
	#cmd.delete("don")
	#cmd.delete("acc")

	# hbonds_all[i] = hbonds
	# print(hbonds_all)

	# save_path = "/Users/zion/Desktop/Rosetta/Kortemme_Lab/Hbonds/"
	# fileOut = os.path.join(save_path, "hbonds.txt")
	# f = open( fileOut, 'w' )		
	# f.write( repr(hbonds_all) + '\n' )
	# f.close()

	return hbonds

cmd.extend("hbnets", hbnets)

