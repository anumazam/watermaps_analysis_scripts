# watermaps_analysis_scripts
Scripts for pruning, analyzing, and visualizing structural water placement predictions from Rosetta-ICO/ECO on protein structures.

### This pipeline uses the files from semi-explicit solvation in Rosetta as input structures.

*(1) Alignment of solvated structures in PyMOL*

*(2) Deletion of waters outside of a 6 Å shell from the ligand*

delete_waters(ligand_name, ligand1_location, ligand2_location, radius, chains)

PyMOL>delete_waters('IPT','resi 998','resi 999',6)

*(3) Identification of water-mediated hydrogen bonds*

This outputs a text file (‘hbonds.txt’)with a dictionary of dictionaries of the identified hydrogen bonds({ 'filename.pdb': { ('donor', 'acceptor'): distance, ..., } , {'filename2.pdb': ...} }) and a separate session file for each solvated structure with the identified hydrogen bonds.

PyMOL>run_hbnets()

*(4) Calculation of the probability of each residue or ligand atom forming a water-mediated hydrogen bond*

bin_waters.py ‘hbonds.txt’ 100

Example output listing each residue (Res) and the probability (Prob) of that residue forming a water-mediated hydrogen bond:

| Res | Prob |
|-----|----- |
|O3NPCC659 |1.00|
|NASPA148 |1.00|
|NTYRA125 |0.95|
|O2NPCC659| 0.60|
|OILEA123 |1.00|
|OLEUA147 |1.00|
|OTYRA125 |1.00|
|NALAA74  |0.10|
|O3NPDD660| 1.00|
|O2NPDD660 |1.00|
|OILEB451| 1.00|
|OTYRB453 |1.00|
|NASPB476| 1.00|
|NLEUB475| 0.55|

The Res column names are organized as follows: 'name' + 'resn' + 'chain' + 'resi', so O3NPCC659 is atom O3 on residue NPC on chain C, residue number 659. In preparation for (6), save the dataframes into separate text files, one for the residues and one for the ligands ('res_hbonds.txt' and 'lig_hbonds.txt'), and add a space between each descriptive identifier in the Res column (O3NPCC659 --> O3 NPC C 659).

*(5) Comparison of hydrogen bonding probabilities between chains and among liganded states can be done with a simple bar plot of the probabilities.*

*(6) Generation of WaterMaps in PyMOL to visualize hydrogen bond propensities.*

In PyMOL, color the input structure gray80 and run data2bfactor.py and color_b.py (adapted from R. Campbell’s online PyMOL script repository, http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/), optionally modifying the color gradient and number of bins. To color the residues by probability, run data2b_res('object','filename_hbonds.txt’) and color_b('res_to_color',minimum='0.1',maximum='1').

To color atoms by probability in the ligand, run data2b_atom('object','lig_hbonds.txt') and color_b('lig_color',minimum='0.1',maximum='1'). Residues or atoms with a probability of 0 are ignored in the WaterMap.

To better visualize the residues, copy 'res_to_color' to a separate object and execute the following in the PyMOL command line:

hide everything, res_to_color_obj
show surface, res_to_color_obj
set transparency, 0.5
show sticks, res_to_color_obj
set cartoon_transparency, 0.7
