#!/bin/bash
#
# Protein-ligand MD setup and run with GROMACS and CHARMM36/CGenFF.
# Run from a folder containing protein.pdb (ligand residue named UNK), the .mdp
# files from this folder, and the third-party files listed in the README.

missing=0
for f in cgenff_charmm2gmx_py3_nx2.py sort_mol2_bonds.pl charmm36-jul2022.ff/forcefield.itp protein.pdb; do
    if [ ! -e "$f" ]; then
        echo "Missing: $f"
        missing=1
    fi
done
if [ "$missing" -ne 0 ]; then
    echo "See the README (Third-party files) for where to download the missing files."
    exit 1
fi

# Step 1: Extract and clean the protein
grep UNK protein.pdb > unk.pdb
grep -v "UNK" protein.pdb > clean.pdb

# Step 2: Process the clean protein structure
gmx pdb2gmx -f clean.pdb -o processed.gro -ter -ignh << EOF
1
1
0
0
0
0
EOF

# Step 3: Convert UNK to mol2 format and fix
obabel unk.pdb -O unk.mol2 --addh
sed -i '2s/.*/UNK/; s/UNK0/UNK/g' unk.mol2
perl sort_mol2_bonds.pl unk.mol2 unk_fix.mol2

# Step 4: Process with CGenFF server
# Manually upload unk_fix.mol2 to the CGenFF server and download unk_fix.str
if [ ! -f unk_fix.str ]; then
    echo "Upload unk_fix.mol2 to the CGenFF server, save the result as unk_fix.str"
    echo "in this folder, then run this script again."
    exit 1
fi

# Step 5: Convert to GROMACS-compatible CHARMM format
python cgenff_charmm2gmx_py3_nx2.py UNK unk_fix.mol2 unk_fix.str charmm36-jul2022.ff

# Step 6: Convert UNK to GROMACS format
gmx editconf -f unk_ini.pdb -o unk.gro

# Step 7: Merge the protein and ligand

# Step 2: Create a copy of 3HTB_processed.gro as complex.gro
cp processed.gro complex.gro

# Step 3: Extract coordinates from unk.gro (skipping header and atom count)
tail -n +3 unk.gro | head -n -1 > unk_coords.gro

# Step 4: Extract the atom count from the second line of complex.gro
CURRENT_COUNT=$(sed -n '2p' complex.gro | tr -d ' ')

# Step 5: Count the number of atoms in unk.gro
NEW_ATOMS=$(wc -l < unk_coords.gro)

# Step 6: Calculate the new atom count
UPDATED_COUNT=$((CURRENT_COUNT + NEW_ATOMS))

# Step 7: Insert ligand coordinates into complex.gro
{
  # Copy all lines except the box vectors from complex.gro
  head -n -1 complex.gro

  # Append ligand coordinates
  cat unk_coords.gro

  # Add the box vectors from complex.gro
  tail -n 1 complex.gro
} > temp_complex.gro

# Step 8: Update the atom count in the second line
sed -i "2s/.*/   $UPDATED_COUNT/" temp_complex.gro

# Step 9: Replace complex.gro with the updated file
mv temp_complex.gro complex.gro

# Cleanup
rm unk_coords.gro

echo "complex.gro successfully updated. Total atoms: $UPDATED_COUNT."


# Step 8: Update topol.top

sed -i '/; Include Position restraint file/{:a;N;/#endif/!ba;s/\n#endif\n/\n#endif/;s/#endif/#endif\n\n; Include ligand topology\n#include "unk.itp"/}' topol.top

sed -i '/; Include forcefield parameters/{:a;N;/#include ".\/charmm36-jul2022.ff\/forcefield.itp"/!ba;s/#include ".\/charmm36-jul2022.ff\/forcefield.itp"/&\n\n; Include ligand parameters\n#include "unk.prm"/}' topol.top

sed -i '/\[ molecules \]/,/Protein_chain_A/ {/Protein_chain_A/ a\
UNK                 1
}' topol.top

# Step 9: Set up the simulation box and solvate
gmx editconf -f complex.gro -o newbox.gro -bt cubic -c -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

# Step 10: Add ions
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 <<EOF
15
EOF

# Step 11: Energy minimization
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# Step 12: Generate position restraints for the ligand
gmx make_ndx -f unk.gro -o index_unk.ndx << EOF
0 & ! a H*
q
EOF

gmx genrestr -f unk.gro -n index_unk.ndx -o posre_unk.itp -fc 1000 1000 1000 <<EOF
3
EOF

sed -i '/; Include water topology/{x;s/.*/\n; Ligand position restraints\n#ifdef POSRES\n#include "posre_unk.itp"\n#endif\n/;G}' topol.top


# Step 13: Create index file for the complex
gmx make_ndx -f em.gro -o index.ndx << EOF
1 | 13
q
EOF

# Step 14: NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt

# Step 15: NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt

# Step 16: Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_100.tpr
gmx mdrun -deffnm md_0_100
