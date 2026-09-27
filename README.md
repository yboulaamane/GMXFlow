# GMXPlotter

A lightweight, reproducible workflow for protein-ligand molecular dynamics with GROMACS: scripts to set up and run the simulation, a script to run the standard analyses, and a Google Colab notebook that plots the results.

![Picture2](https://github.com/user-attachments/assets/7454382d-9afc-47e3-b9bb-8c0b9441e0af)

| Step | Where | What it does |
|---|---|---|
| 1. Set up and run | `runner/gromacs-runner.sh` | Protein topology (CHARMM36), ligand parameters (CGenFF), complex assembly, solvation, ions, minimisation, NVT and NPT equilibration, production MD |
| | `runner/submit_gromacs_gpu.sh` | Example SLURM job for running the production step on a GPU node |
| 2. Analyse | `runner/gromacs-analyzer.sh` | PBC correction and fitting, RMSD (backbone and ligand), RMSF, radius of gyration, SASA, PCA and the free energy landscape |
| 3. Plot | `GMXPlotter.ipynb` [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/yboulaamane/GMXPlotter/blob/main/GMXPlotter.ipynb) | Publication-style plots of the analysis output. The folders at the top level hold example data. |

## Running a simulation

Requirements: GROMACS, Open Babel, Perl and Python 3.

1. Put `protein.pdb` in a working folder, with the ligand as residue `UNK`.
2. Copy the `.mdp` files and scripts from `runner/`, plus the third-party files below, into the same folder.
3. Run `bash gromacs-runner.sh`. It stops after writing `unk_fix.mol2`: upload that file to the CGenFF server, save the result as `unk_fix.str` in the folder, and run the script again.
4. Run `bash gromacs-analyzer.sh` on the finished trajectory, then open the notebook to plot the `.xvg` output.

The scripts use fixed GROMACS group numbers (for example `1 | 13` for protein plus ligand, `15` for the solvent) and a 100 ns trajectory. Check them against your own system's index groups before running. The `module load` lines are for one specific HPC cluster; change or remove them.

## Third-party files

These are needed but not included in this repository; download them from their authors:

| File | Source |
|---|---|
| `charmm36-jul2022.ff` (unpacked folder) | CHARMM36 force field for GROMACS, MacKerell lab: https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs |
| `cgenff_charmm2gmx_py3_nx2.py` | CGenFF to GROMACS conversion script, same MacKerell lab page (GNU AGPL-3.0) |
| `sort_mol2_bonds.pl` | Justin Lemkul's GROMACS protein-ligand tutorial: http://www.mdtutorials.com/gmx/complex/ |
| `xpm2txt.py` (optional) | Converts the free energy landscape (`FES.xpm`) to text; the analysis script skips that step without it |

The `.mdp` parameter files are adapted from Justin Lemkul's protein-ligand complex tutorial.
