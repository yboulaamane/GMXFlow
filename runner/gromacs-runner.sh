#!/bin/bash
#
# Protein-ligand MD with GROMACS, CHARMM36 and CGenFF: from a PDB of the
# complex to a production trajectory.
#
# Run it in a folder that contains the input PDB, the .mdp files from this
# folder and the third-party files listed in the README. Settings come from
# environment variables, for example:
#
#     LIG=LIG PROTEIN=my_complex.pdb MD_NS=100 bash gromacs-runner.sh
#
# The ligand needs a parameter file from the CGenFF server, so the script stops
# once for you to fetch it, then carries on when you run it again.

set -euo pipefail

PROTEIN=${PROTEIN:-protein.pdb}     # protein-ligand complex
LIG=${LIG:-UNK}                     # ligand residue name in that PDB
FF=${FF:-charmm36-jul2022}          # force field folder, without ".ff"
WATER=${WATER:-tip3p}               # water model defined by that force field
BOX=${BOX:-cubic}                   # box type (gmx editconf -bt)
DIST=${DIST:-1.0}                   # nm between the solute and the box edge
CONC=${CONC:-0.15}                  # mol/L NaCl on top of neutralisation
DEFFNM=${DEFFNM:-md_0_100}          # file name for the production run
MD_NS=${MD_NS:-}                    # production length in ns (empty: as set in md.mdp)
GMX=${GMX:-gmx}                     # GROMACS binary, e.g. gmx_mpi
GMX_MODULE=${GMX_MODULE:-}          # optional environment module to load on a cluster

# The CGenFF conversion script names its output files in lower case.
lig=$(echo "$LIG" | tr '[:upper:]' '[:lower:]')

die() { echo "Error: $*" >&2; exit 1; }

if [ -n "$GMX_MODULE" ]; then
    module load "$GMX_MODULE"
fi

missing=0
for f in "$PROTEIN" cgenff_charmm2gmx_py3_nx2.py sort_mol2_bonds.pl "$FF.ff/forcefield.itp" \
         ions.mdp em.mdp nvt.mdp npt.mdp md.mdp; do
    if [ ! -e "$f" ]; then
        echo "Missing: $f" >&2
        missing=1
    fi
done
for tool in "$GMX" obabel perl python; do
    command -v "$tool" >/dev/null || { echo "Missing program: $tool" >&2; missing=1; }
done
[ "$missing" -eq 0 ] || die "see the README (Third-party files) for where to get what is missing."

# Step 1: split the complex. The ligand is every ATOM/HETATM record whose residue
# name (PDB columns 18-21) is $LIG; the protein is every other ATOM record.
# Other HETATM records (crystal waters, ions, additives) are left out.
resname='rn = substr($0, 18, 4); gsub(/ /, "", rn)'
awk -v r="$LIG" "/^(ATOM|HETATM)/ { $resname; if (rn == r) print }" "$PROTEIN" > "${lig}.pdb"
awk -v r="$LIG" "/^ATOM/ { $resname; if (rn != r) print } /^(TER|END)/" "$PROTEIN" > clean.pdb
[ -s "${lig}.pdb" ] || die "no atoms with residue name $LIG in $PROTEIN. Set LIG to the ligand's residue name."
grep -q '^ATOM' clean.pdb || die "no protein ATOM records in $PROTEIN."

# Step 2: protein topology, with the force field and water model named rather
# than picked by their position in a menu, and the default termini.
"$GMX" pdb2gmx -f clean.pdb -o processed.gro -ff "$FF" -water "$WATER" -ignh

# Step 3: ligand to mol2 with hydrogens, one residue name throughout, bonds sorted for CGenFF
obabel "${lig}.pdb" -O "${lig}.mol2" --addh
sed -i "2s/.*/${LIG}/; s/${LIG}[0-9]\{1,\}/${LIG}/g" "${lig}.mol2"
perl sort_mol2_bonds.pl "${lig}.mol2" "${lig}_fix.mol2"

# Step 4: ligand parameters from the CGenFF server (manual step)
if [ ! -f "${lig}_fix.str" ]; then
    echo "Upload ${lig}_fix.mol2 to the CGenFF server, save the result as ${lig}_fix.str"
    echo "in this folder, then run this script again."
    exit 0
fi

# Step 5: CGenFF stream file to GROMACS: writes ${lig}.itp, ${lig}.prm and ${lig}_ini.pdb
python cgenff_charmm2gmx_py3_nx2.py "$LIG" "${lig}_fix.mol2" "${lig}_fix.str" "$FF.ff"
"$GMX" editconf -f "${lig}_ini.pdb" -o "${lig}.gro"

# Step 6: complex.gro = protein atoms, then ligand atoms, then the protein's box line
protein_atoms=$(sed -n '2p' processed.gro | tr -d ' ')
ligand_atoms=$(sed -n '2p' "${lig}.gro" | tr -d ' ')
{
    sed -n '1p' processed.gro
    echo "$((protein_atoms + ligand_atoms))"
    sed -n "3,$((protein_atoms + 2))p" processed.gro
    sed -n "3,$((ligand_atoms + 2))p" "${lig}.gro"
    tail -n 1 processed.gro
} > complex.gro
echo "complex.gro: $protein_atoms protein atoms + $ligand_atoms ligand atoms"

# Step 7: topology. Ligand parameters go after the force field include; the
# ligand topology and its position restraints go before the water include; the
# ligand molecule is listed last. This works for one or many protein chains.
# (The restraint file is written in step 11; it is only read when POSRES is defined.)
grep -q 'forcefield.itp"' topol.top || die "topol.top has no force field include."
grep -q '; Include water topology' topol.top || die "topol.top has no water topology include."
sed -i "/forcefield.itp\"/a\\
\\
; Include ligand parameters\\
#include \"${lig}.prm\"" topol.top
sed -i "/; Include water topology/i\\
; Include ligand topology\\
#include \"${lig}.itp\"\\
\\
; Ligand position restraints\\
#ifdef POSRES\\
#include \"posre_${lig}.itp\"\\
#endif\\
" topol.top
printf '%-20s1\n' "$LIG" >> topol.top

# Step 8: box and water
"$GMX" editconf -f complex.gro -o newbox.gro -bt "$BOX" -c -d "$DIST"
"$GMX" solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

# Step 9: neutralise and add salt, replacing water molecules (group chosen by name)
"$GMX" grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL | "$GMX" genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc "$CONC"

# Step 10: energy minimisation
"$GMX" grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
"$GMX" mdrun -v -deffnm em

# Step 11: position restraints on the ligand's heavy atoms. In a ligand-only
# structure the default groups are System, Other and the residue, so the new
# heavy-atom group is number 3; it is renamed and then selected by name.
printf '0 & ! a H*\nname 3 LIG_heavy\nq\n' | "$GMX" make_ndx -f "${lig}.gro" -o "index_${lig}.ndx"
echo LIG_heavy | "$GMX" genrestr -f "${lig}.gro" -n "index_${lig}.ndx" -o "posre_${lig}.itp" -fc 1000 1000 1000

# Step 12: index group Protein_<LIG> (protein plus ligand) for temperature coupling
printf '"Protein" | "%s"\nq\n' "$LIG" | "$GMX" make_ndx -f em.gro -o index.ndx

# The .mdp files couple "Protein_UNK"; write copies for this ligand's name and,
# if MD_NS is set, the requested production length.
for m in nvt npt md; do
    sed "s/Protein_UNK/Protein_${LIG}/" "$m.mdp" > "$m.run.mdp"
done
if [ -n "$MD_NS" ]; then
    dt=$(awk -F'[=;]' '$1 ~ /^ *dt *$/ { gsub(/ /, "", $2); print $2 }' md.mdp)
    nsteps=$(awk -v ns="$MD_NS" -v dt="$dt" 'BEGIN { printf "%d", ns * 1000 / dt }')
    sed -i "s/^\( *nsteps *=\).*/\1 $nsteps ; $MD_NS ns at $dt ps per step (set by MD_NS)/" md.run.mdp
    echo "Production run: $MD_NS ns = $nsteps steps of $dt ps"
fi

# Step 13: NVT equilibration
"$GMX" grompp -f nvt.run.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
"$GMX" mdrun -deffnm nvt

# Step 14: NPT equilibration
"$GMX" grompp -f npt.run.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
"$GMX" mdrun -deffnm npt

# Step 15: production MD
"$GMX" grompp -f md.run.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o "$DEFFNM.tpr"
"$GMX" mdrun -deffnm "$DEFFNM"
