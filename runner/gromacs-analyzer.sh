#!/bin/bash

# Load GROMACS module if using HPC (replace with your Gromacs version)
module load GROMACS/2021.3-foss-2021a

# Step 1: Remove PBC and center
gmx trjconv -s md_0_100.tpr -f md_0_100.xtc -o md_pbc.xtc -center -pbc mol -ur compact << EOF
1
1
EOF

# Step 2: Fit to backbone
gmx trjconv -s md_0_100.tpr -f md_pbc.xtc -o md_0_100_fit.xtc -fit rot+trans << EOF
4
1
EOF

# Step 3: Extract start structure (t = 0)
gmx trjconv -s md_0_100.tpr -f md_0_100_fit.xtc -o start.pdb -dump 0 << EOF
14
EOF

# Step 4: Extract end structure (t = 100000 ps)
gmx trjconv -s md_0_100.tpr -f md_0_100_fit.xtc -o end.pdb -dump 100000 << EOF
14
EOF

# Step 5: RMSD - backbone
gmx rms -s md_0_100.tpr -f md_0_100_fit.xtc -o rmsd_backbone.xvg -fit rot+trans -tu ns << EOF
4
4
EOF

# Step 6: RMSD - ligand
gmx rms -s md_0_100.tpr -f md_0_100_fit.xtc -o rmsd_ligand.xvg -fit rot+trans -tu ns << EOF
13
13
EOF

# Step 7: RMSF - backbone
gmx rmsf -s md_0_100.tpr -f md_0_100_fit.xtc -o rmsf_residues.xvg -ox bfactor.pdb << EOF
4
EOF

# Step 8: Radius of gyration
gmx gyrate -s md_0_100.tpr -f md_0_100_fit.xtc -o gyrate.xvg << EOF
1
EOF

# Step 9: SASA
gmx sasa -s md_0_100.tpr -f md_0_100_fit.xtc -o sasa.xvg << EOF
1
EOF

# Step 10: Covariance matrix (PCA)
gmx covar -s md_0_100.tpr -f md_0_100_fit.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm << EOF
3
3
EOF

# Step 11: Project onto top 2 eigenvectors
gmx anaeig -s md_0_100.tpr -f md_0_100_fit.xtc -v eigenvectors.trr -first 1 -last 2 -proj pc1_pc2.xvg << EOF
3
3
EOF

# Step 12: Create input for FES
paste <(awk 'NR>17 {print $1, $2}' pc1_pc2.xvg) <(awk 'NR>17 {print $2}' pc1_pc2.xvg) > PC1PC2.xvg

# Step 13: Free Energy Surface
gmx sham -f PC1PC2.xvg -ls FES.xpm

# Step 14: Convert FES to text (needs xpm2txt.py, not included; see the README)
if [ -f xpm2txt.py ]; then
    python xpm2txt.py -f FES.xpm -o fel.dat
else
    echo "xpm2txt.py not found: skipped converting FES.xpm to fel.dat"
fi
