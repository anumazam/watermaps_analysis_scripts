# Run hbnets.py over all the PDB files in a directory

directory_hbnets = "/Users/zion/Desktop/Rosetta/Kortemme_Lab/Hbonds"

import os
os.chdir(directory_hbnets)
exec(open("hbnets.py").read());
from pymol import cmd

def run_hbnets():

	##########Make sure path to directories are correct!!!!##########
	directory_files = "/Users/zion/Desktop/Rosetta/Kortemme_Lab/_3.Design/CoupledMoves19.07.18loop_output/1efaIPTGALLAAnearwater"
	save_directory = "/Users/zion/Desktop/Rosetta/Kortemme_Lab/_3.Design/CoupledMoves19.07.18loop_output/1efaIPTGALLAAnearwater"
	
	
	directory_hbnets = "/Users/zion/Desktop/Rosetta/Kortemme_Lab/Hbonds"


	hbonds_all = {}

	for file in os.listdir(directory_files):
		if file.endswith('.pdb'):
			os.chdir(directory_files)
			cmd.load(file)
			os.chdir(directory_hbnets)
			hbonds = hbnets(file, directory_files)

			hbonds_all[file] = hbonds

			#save session file of hbonds
			os.chdir(save_directory)
			cmd.save(file + "_hbonds.pse")
			cmd.delete("all")
	
	print(hbonds_all)

	#save dictionary of dictionaries of hbonds
	fileOut = os.path.join(save_directory, "hbonds.txt")
	f = open( fileOut, 'w' )		
	f.write( repr(hbonds_all) + '\n' )
	f.close()


cmd.extend("run_hbnets", run_hbnets)




