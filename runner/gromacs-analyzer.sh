#!/bin/bash
#
# Standard analyses of a protein-ligand trajectory from gromacs-runner.sh.
# Run it in the same folder. Settings come from environment variables, e.g.
#
#     LIG=LIG DEFFNM=md_0_100 bash gromacs-analyzer.sh
#
# Groups are selected by name, so the script works whatever the index numbers
# of your system are. Output .xvg files are what GMXPlotter.ipynb plots.

set -euo pipefail

DEFFNM=${DEFFNM:-md_0_100}          # production run name (the .tpr and .xtc)
LIG=${LIG:-UNK}                     # ligand residue name
INDEX=${INDEX:-index.ndx}           # index file with the Protein_<LIG> group
GMX=${GMX:-gmx}                     # GROMACS binary, e.g. gmx_mpi
GMX_MODULE=${GMX_MODULE:-}          # optional environment module to load on a cluster

die() { echo "Error: $*" >&2; exit 1; }

if [ -n "$GMX_MODULE" ]; then
    module load "$GMX_MODULE"
fi
for f in "$DEFFNM.tpr" "$DEFFNM.xtc"; do
    [ -f "$f" ] || die "missing $f; set DEFFNM to your production run name."
done

# The runner writes the index file; recreate it if it is missing.
if [ ! -f "$INDEX" ]; then
    printf '"Protein" | "%s"\nq\n' "$LIG" | "$GMX" make_ndx -f "$DEFFNM.tpr" -o "$INDEX"
fi
sel() { printf '%s\n' "$@"; }      # answers to GROMACS group prompts, one per line
n=(-n "$INDEX")

# Step 1: make molecules whole and centre the protein (whole system kept, so the ligand stays in)
sel Protein System | "$GMX" trjconv -s "$DEFFNM.tpr" -f "$DEFFNM.xtc" "${n[@]}" -o md_pbc.xtc -center -pbc mol -ur compact

# Step 2: remove rotation and translation by fitting to the protein backbone
sel Backbone System | "$GMX" trjconv -s "$DEFFNM.tpr" -f md_pbc.xtc "${n[@]}" -o "${DEFFNM}_fit.xtc" -fit rot+trans

# Steps 3-4: first and last frames of protein plus ligand
last=$("$GMX" check -f "${DEFFNM}_fit.xtc" 2>&1 | awk '/^Last frame/ { print $NF }')
[ -n "$last" ] || die "could not read the trajectory's last frame time."
sel "Protein_${LIG}" | "$GMX" trjconv -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o start.pdb -dump 0
sel "Protein_${LIG}" | "$GMX" trjconv -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o end.pdb -dump "$last"

# Step 5: backbone RMSD
sel Backbone Backbone | "$GMX" rms -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o rmsd_backbone.xvg -fit rot+trans -tu ns

# Step 6: ligand RMSD (ligand fitted on itself)
sel "$LIG" "$LIG" | "$GMX" rms -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o rmsd_ligand.xvg -fit rot+trans -tu ns

# Step 7: backbone RMSF
sel Backbone | "$GMX" rmsf -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o rmsf_residues.xvg -ox bfactor.pdb

# Step 8: radius of gyration
sel Protein | "$GMX" gyrate -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o gyrate.xvg

# Step 9: solvent-accessible surface area
sel Protein | "$GMX" sasa -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o sasa.xvg

# Step 10: covariance matrix of the C-alpha atoms (PCA)
sel C-alpha C-alpha | "$GMX" covar -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm

# Step 11: projection on the first two eigenvectors
sel C-alpha C-alpha | "$GMX" anaeig -s "$DEFFNM.tpr" -f "${DEFFNM}_fit.xtc" "${n[@]}" -v eigenvectors.trr -first 1 -last 2 -proj pc1_pc2.xvg

# Step 12: time, PC1, PC2 columns for the free energy surface. anaeig writes one
# data set per eigenvector, separated by "&"; the header length is not fixed.
if grep -q '^&' pc1_pc2.xvg; then
    paste <(awk '/^[@#]/ { next } /^&/ { exit } { print $1, $2 }' pc1_pc2.xvg) \
          <(awk '/^[@#]/ { next } /^&/ { s = 1; next } s { print $2 }' pc1_pc2.xvg) > PC1PC2.xvg
else
    awk '!/^[@#]/ { print $1, $2, $3 }' pc1_pc2.xvg > PC1PC2.xvg
fi

# Step 13: free energy surface
"$GMX" sham -f PC1PC2.xvg -ls FES.xpm

# Step 14: FES to text (needs xpm2txt.py, not included; see the README)
if [ -f xpm2txt.py ]; then
    python xpm2txt.py -f FES.xpm -o fel.dat
else
    echo "xpm2txt.py not found: skipped converting FES.xpm to fel.dat"
fi
