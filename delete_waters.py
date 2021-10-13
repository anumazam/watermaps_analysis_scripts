#!/usr/bin/python
# find and delete water molecules that are not within a certain radius of a ligand atom
# assumes ligand is the only residue in chain X.

from pymol import cmd

# pymol function to find and delete the water molecules
def delete_waters(ligand_name, ligand1_location, ligand2_location, radius, chains):
    '''
    Inputs
    ligand_name = ex:'IPT'
    ligand1_location = ex:'resi 998'
    ligand2_location = ex:'resi 999'
    radius = ex:6, in Angstroms
    chains = ex: ['A','B']

    '''

    lig1_waters = "lig1_waters"
    lig2_waters = "lig2_waters"
    waters_to_delete = "waters_far_" + ligand_name
    
    # select the correct water atoms
    cmd.select(lig1_waters, ligand1_location + " around %d and hetatm and chain %s" %(radius, chains[0]))
    cmd.select(lig2_waters, ligand2_location + " around %d and hetatm and chain %s" %(radius, chains[1]))
    cmd.select(waters_to_delete, "hetatm and not " + lig1_waters + " and not " + lig2_waters + " and not " + ligand1_location + " and not " + ligand2_location)
    
    # delete the water atoms that are far from the ligand(s)
    cmd.remove(waters_to_delete)
    
    # show remaining waters as spheres and ligands as sticks
    cmd.show_as("spheres", lig1_waters)
    cmd.show_as("spheres", lig2_waters)
    cmd.set ("sphere_scale", 0.2)
    cmd.show_as("sticks", ligand1_location)
    cmd.show_as("sticks", ligand2_location)
    cmd.hide("everything", "hydrogen")
    
    # remove selection
    cmd.delete(waters_to_delete)

# let pymol know about this function
cmd.extend("delete_waters", delete_waters)



