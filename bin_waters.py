import pandas as pd
import sys

# Calculate the probability (from 1/n to 1) of each residue forming a water-mediated hydrogen bond. Ignores any residues not making a water-mediated hydrogen bond within the 6 A shell.


'''
Inputs:

filename = string filename of the calculated hydrogen bonds from run_hbnets(), ex: 'hbonds.txt'
nstruct = integer of the number of structures for which you are binning the waters

'''
filename = sys.argv[1]
nstruct = sys.argv[2]


from collections import defaultdict

# Create a pandas dataframe of hydrogen bonds
# originally the function get_structures(filename)
hbonds = eval(open(filename).read())
all_structs = {}
for structure in hbonds:
    struct_dict = hbonds[structure]
    df = pd.Series(struct_dict).reset_index()
    df.columns = ['Res1','Res2','Bond Length']
    all_structs[structure] = df

water_map_A = defaultdict(int)
water_map_B = defaultdict(int)
hbonds_A = defaultdict(list)
hbonds_B = defaultdict(list)

residues = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
ligands = ['IPT','NPF','IPC','IPD','NPC','NPD'] 
waters = ['HOH','WAT','TP3']   

for struct in all_structs:
    struct_df = all_structs[struct]
    found_atoms = []
    for res in range(0, len(struct_df['Res1'])):
        res1 = struct_df['Res1'][res]
        res2 = struct_df['Res2'][res]
        bl = struct_df['Bond Length'][res]
        for wat in waters:
            if wat in res2:
                count = 0
                found_atoms.append(res1)
                for atm in found_atoms:
                    if atm == res1:
                        count += 1
                if count == 1:
            #only count a residue once per structure
                    for aa in residues:
                        if aa in res1:
                            chain_letter = res1[res1.find(aa) + 3 ]
                            break
                    for lig in ligands:
                        if lig in res1:
                            chain_letter = res1[res1.find(lig) + 3 ]
                            break
                    if chain_letter == 'A' or chain_letter == 'C':
                        water_map_A[res1] += 1
                        hbonds_A[res1].append((res2,bl))
                    if chain_letter == 'B' or chain_letter == 'D':
                        water_map_B[res1] += 1
                        hbonds_B[res1].append((res2,bl))
        for wat in waters:       
            if wat in res1:
                count = 0
                found_atoms.append(res2)
                for atm in found_atoms:
                    if atm == res2:
                        count += 1
                if count == 1:
            #only count a residue once per structure
                    for aa in residues:
                        if aa in res2:
                            chain_letter = res2[res2.find(aa) + 3 ]
                            break
                    for lig in ligands:
                        if lig in res2:
                            chain_letter = res2[res2.find(lig) + 3 ]
                            break
                    if chain_letter == 'A' or chain_letter == 'C':
                        water_map_A[res2] += 1
                        hbonds_A[res2].append((res1,bl))
                    if chain_letter == 'B' or chain_letter == 'D':
                        water_map_B[res2] += 1
                        hbonds_B[res2].append((res1,bl))
                        
water_map_A = pd.DataFrame.from_dict(water_map_A, orient='index').reset_index()
water_map_B = pd.DataFrame.from_dict(water_map_B, orient='index').reset_index()
water_map_A.columns = ['Res','Prob']
water_map_A['Prob'] = water_map_A['Prob']/100
water_map_B.columns = ['Res','Prob']
water_map_B['Prob'] = water_map_B['Prob']/100
print(water_map_A)
print(water_map_B)


